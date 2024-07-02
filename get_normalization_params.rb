# for MOTIF_LENGTH in $(seq 3 42); do
#   echo ruby get_normalization_params.rb --mutational-signature 'TOPMed_10kb_spectra.tsv:*' \
#        --genomic-context genomic_contexts/hg38_contexts.tsv --motif-lengths ${MOTIF_LENGTH} \
#        --motifs ./motifs_freeze --iterations 10000 \
#        \> background_metrics/bg_metrics_len_${MOTIF_LENGTH}.json;
# done | /usr/bin/time -v parallel
# ruby -rjson -e 'puts Dir.glob("background_metrics/*").map{|fn| JSON.parse(File.read(fn)) }.inject(&:merge).sort.to_h.to_json' > background_metrics.json

require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/context'
require_relative 'lib/mutation_process'
# require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'
require_relative 'lib/score_to_pvalue'
require_relative 'basic_stats'
require_relative 'utils'
require 'optparse'
require 'json'

def information_content(pos)
  2 + pos.map{|freq| (freq == 0) ? 0 : Math.log2(freq) * freq }.sum
end

num_iterations = 100

contexts_wg = nil
motifs_folder = nil

signatures = {}
motif_lengths = (6..15).to_a

option_parser = OptionParser.new{|opts|
  CLI.mutational_context_option(opts)
  CLI.genomic_context_option(opts)
  CLI.mutational_signature_option(opts)

  opts.on('--motif-lengths LENGTH', 'Specify motif lengths. E.g. 5,9-12. Default: 6-15.'){|str|
    motif_lengths = str.split(',').flat_map{|range|
      vals = range.split('-').map{|v| Integer(v) }
      if vals.size == 1
        [vals[0]]
      elsif vals.size == 2
        (vals[0] .. vals[1]).to_a
      else
        raise
      end
    }
  }

  opts.on('--iterations N', 'Specify number of iterations for bootstraping (default: #{num_iterations})'){|value|
    num_iterations = Integer(value)
  }

  opts.on('--motifs FOLDER', 'Specify folder with motifs to sample positions from them') {|fn|
    motifs_folder = fn
  }
}

option_parser.parse!(ARGV)

raise 'Specify motifs folder'  unless motifs_folder
raise 'Specify genomic context'  unless contexts_wg
raise 'Specify mutation process(es)'  unless signatures.size > 0

mutational_processes = signatures.transform_values{|signature|
  MutationProcess.new(signature, contexts_wg).symmetrized
}

original_motifs = Dir.glob("#{motifs_folder}/*").map{|fn|
  File.readlines(fn).drop(1).map{|l|
    counts = l.chomp.split("\t").map{|x| Float(x) }
    counts.map{|cnt| cnt.to_f / counts.sum }
  }
}; nil
positions = original_motifs.flatten(1)
ic_threshold = positions.map{|pos| information_content(pos) }.quantile(0.25)
positions = positions.select{|pos| information_content(pos) <= ic_threshold }.shuffle; nil


result = motif_lengths.map{|motif_length|
  $stderr.puts "Motif length: #{motif_length}"
  simulated_motif_infos = num_iterations.times.map{
    $stderr.print '.'
    # positions = original_motifs.select{|matrix| matrix.length == motif_length }.flatten(1)
    # pcm = PCM.new(original_motifs.sample.shuffle).to_ppm.to_pcm(100)
    # pcm = PCM.new(original_motifs.sample).to_ppm.to_pcm(100)
    pcm = PCM.new(positions.shuffle.first(motif_length)).to_ppm.to_pcm(100)
    ppm = pcm.to_ppm
    pwm = pcm.to_pwm(pseudocount: :log)
    pwm_fn = write_temp_matrix(pwm)

    thresholds_fn ||= make_thresholds(pwm_fn, thresholds_fn: nil)
    bsearch_table = read_bsearch_table(thresholds_fn)
    {pcm: pcm, ppm: ppm, pwm: pwm, pwm_fn: pwm_fn, thresholds_fn: thresholds_fn, bsearch_table: bsearch_table}
  }

  background_by_signature = mutational_processes.map{|signature_name, mutational_process|
    $stderr.print '+'
    background_infos = simulated_motif_infos.map{|motif_info|
      mutation_action_infos(mutational_process, motif_info[:ppm], motif_info[:pwm], motif_info[:bsearch_table])
    }

    metrics = background_infos.first.keys
    background_stats = metrics.map{|metric_name|
      metric_values = background_infos.map{|info| info[metric_name] }
      [metric_name, {mean: metric_values.mean, stddev: metric_values.stddev}]
    }.to_h

    [signature_name, background_stats]
  }.to_h

  $stderr.puts
  [motif_length, background_by_signature]
}.to_h

puts result.to_h.to_json
