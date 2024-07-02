require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/context'
require_relative 'lib/mutation_process'
# require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'
require_relative 'lib/score_to_pvalue'
require_relative 'basic_stats'
require_relative 'utils'
require 'json'
require 'optparse'

signatures = {}

contexts_wg = nil

original_pcm = nil
original_ppm = nil
original_pwm = nil
original_bsearch_table = nil
pcm_count = nil
pseudocount = :log
background_metrics_fn = nil

option_parser = OptionParser.new{|opts|
  CLI.mutational_context_option(opts)
  CLI.genomic_context_option(opts)
  CLI.mutational_signature_option(opts)
  
  opts.on('--pcm-motif FILE', 'Specify PFM/PCM file'){|fn|
    original_pcm = PCM.from_file(fn)
  }
  opts.on('--pfm-motif FILE', 'Specify PFM/PCM file'){|fn|
    original_ppm = PCM.from_file(fn).to_ppm
  }
  opts.on('--pcm-count FILE', 'Specify PCM count to multiply PFM'){|value|
    pcm_count = Float(value)
  }
  opts.on('--pseudocount VALUE', "Specify PCM pseudocount (certain number, `log` or `sqrt`). Default: #{pseudocount}") {|value|
    pseudocount = Float(value) rescue value.to_sym
  }
  opts.on('--pwm-motif FILE', 'Specify PWM file'){|fn|
    original_pwm = PWM.from_file(fn)
  }
  opts.on('--thresholds FILE', 'Specify thresholds grid file'){|fn|
    original_bsearch_table = read_bsearch_table(fn)
  }
  opts.on('--background-metrics FILE', 'Specify background metrics file. One can use pattern with `{length}` in it'){|fn|
    background_metrics_fn = fn
  }
}

option_parser.parse!(ARGV)


if !original_pcm && !original_ppm
  raise 'Specify either PCM or PFM. PCM is prefferable'
elsif !original_ppm
  original_ppm = original_pcm.to_ppm
elsif !original_pcm
  if pcm_count
    original_pcm = original_ppm.to_pcm(pcm_count)
  else
    raise "Specify PCM count"
  end
end

original_pwm = original_pcm.to_pwm(pseudocount: pseudocount)  if !original_pwm
if !original_bsearch_table
  pwm_fn = write_temp_matrix(original_pwm)
  thresholds_fn ||= make_thresholds(pwm_fn, thresholds_fn: nil)
  original_bsearch_table = read_bsearch_table(thresholds_fn)
end

raise 'Specify background metrics file'  if !background_metrics_fn
background_metrics_fn = background_metrics_fn.gsub('{length}', original_pcm.length.to_s)

mutational_processes = signatures.transform_values{|signature|
  MutationProcess.new(signature, contexts_wg).symmetrized
}

background_metrics_full = JSON.parse(File.read(background_metrics_fn))
background_metrics = background_metrics_full[original_pcm.length.to_s]

all_motif_metrics = mutational_processes.map{|signature, mutational_process|
  signature_background = background_metrics[signature]

  original_motif_info = mutation_action_infos(mutational_process, original_ppm, original_pwm, original_bsearch_table)
  motif_metrics = ['attractiveness', 'pvalue_change', 'corrected_pvalue_change', 'mean_score', 'mean_pv', 'mean_weight_change', 'stddev_weight_change'].map{|key|
    value = original_motif_info[key.to_sym]
    value_mean, value_stddev = signature_background[key].values_at('mean', 'stddev')
    [key, {
        value: value,
        zscore: zscore(value, value_mean, value_stddev),
        distribution: {mean: value_mean, stddev: value_stddev},
    }]
  }.to_h
  [signature, motif_metrics]
}.to_h

puts replace_nan(all_motif_metrics, replacement: nil).to_json
