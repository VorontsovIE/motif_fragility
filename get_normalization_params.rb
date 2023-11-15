require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/context'
require_relative 'lib/mutation_process'
# require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'
require_relative 'lib/score_to_pvalue'
require_relative 'basic_stats'
require 'optparse'
require 'json'

def mutation_action_infos(mutational_process, ppm, pwm, bsearch_table)
  attractiveness = mutational_process.motif_attractiveness(ppm)
  mean_weight_change = mutational_process.mean_weight_change(ppm, pwm)
  weight_stddev = mutational_process.weight_stddev(ppm, pwm)

  mean_score = pwm.mean_score(ppm)
  mean_pv = -Math.log10(pvalue_by_score(mean_score, bsearch_table))
  pv_after = -Math.log10(pvalue_by_score(mean_score + mean_weight_change, bsearch_table))
  pv_after_corrected = -Math.log10(pvalue_by_score(mean_score + attractiveness * mean_weight_change, bsearch_table))
  info = {
    mean_score: mean_score,
    mean_pv: mean_pv,
    attractiveness: attractiveness, #.round(2),
    mean_weight_change: mean_weight_change, #.round(2),
    stddev_weight_change: weight_stddev, #.round(2),
    pvalue_change: (pv_after - mean_pv), #.round(2),
    corrected_pvalue_change: (pv_after_corrected - mean_pv), #.round(2),
  }
  info
end


def information_content(pos)
  2 + pos.map{|freq| (freq == 0) ? 0 : Math.log(freq) * freq }.sum
end

num_iterations = 100
mutation_counts = nil

contexts_wg = nil

motif_length = nil
motifs_folder = nil

option_parser = OptionParser.new{|opts|
  opts.on('--mutational-context FILE', 'Specify mutational process contexts'){|fn|
    # mutation_counts = val
    mutation_counts = ContextDependentCounts.from_file(fn){|ctx|
      DirectedContext.from_string(ctx)
    }
  }
  opts.on('--mutational-signature FILE:SIGNATURE', 'Specify mutational process file and signature name'){|value|
    fn, signature = value.split(':', 2) # 'COSMIC_v3.4_SBS_GRCh38.txt'
    mutation_counts = ContextDependentCounts.from_multisignature_file(fn, signature){|ctx|
      DirectedContext.from_string(ctx)
    }
  }
  opts.on('--genomic-context FILE', 'Specify genomic contexts'){|fn|
    contexts_wg = ContextDistribution.from_file(fn).without_unknown
  }
  opts.on('--motif-length LENGTH', 'Specify motif length'){|value|
    motif_length = Integer(value)
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
raise 'Specify motif length'  unless motif_length
raise 'Specify genomic context'  unless contexts_wg
raise 'Specify mutation process'  unless mutation_counts

# raise 'Specify region specific contexts'  unless region_specific_contexts_fn = ARGV[2]
# contexts_rs_frequencies = ContextDistribution.from_file(region_specific_contexts_fn).without_unknown.symmetrized.frequencies


# mutation_counts = ContextDependentCounts.from_multisignature_file('COSMIC_v3.4_SBS_GRCh38.txt', 'SBS1'){|ctx|
#   DirectedContext.from_string(ctx)
# }

mutational_process = MutationProcess.new(mutation_counts, contexts_wg).symmetrized


motifs = Dir.glob("#{motifs_folder}/*").map{|fn|
  File.readlines(fn).drop(1).map{|l|
    counts = l.chomp.split("\t").map{|x| Float(x) }
    counts.map{|cnt| cnt.to_f / counts.sum }
  }
}; nil
positions = motifs.flatten(1)
ic_threshold = positions.map{|pos| information_content(pos) }.quantile(0.25)
positions = positions.select{|pos| information_content(pos) <= ic_threshold }.shuffle; nil

# motifs.select!{|matrix| matrix.length == motif_length }


background_infos = num_iterations.times.map {
  $stderr.print '.'
  pcm = PCM.new(positions.shuffle.first(motif_length)).to_ppm.to_pcm(100)
  # pcm = PCM.new(motifs.sample.shuffle).to_ppm.to_pcm(100)
  # pcm = PCM.new(motifs.sample).to_ppm.to_pcm(100)
  ppm = pcm.to_ppm
  pwm = pcm.to_pwm(pseudocount: :log)
  pwm_fn = write_temp_matrix(pwm)

  thresholds_fn ||= make_thresholds(pwm_fn, thresholds_fn: nil)
  bsearch_table = read_bsearch_table(thresholds_fn)

  mutation_action_infos(mutational_process, ppm, pwm, bsearch_table)
}
$stderr.puts

metrics = background_infos.first.keys
result = metrics.map{|metric_name|
  metric_values = background_infos.map{|info| info[metric_name] }
  [metric_name, {mean: metric_values.mean, stddev: metric_values.stddev}]
}.to_h

puts result.to_json
