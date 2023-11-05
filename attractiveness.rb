# require 'optparse'
require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/context'
require_relative 'lib/mutation_process'
# require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'
require_relative 'lib/score_to_pvalue'
require 'optparse'
require 'gruff'

module Enumerable
  def mean
    (length == 0) ? nil : sum.to_f / length
  end

  def variance
    m = mean
    (length < 2) ? nil : map{|x| (x - m) ** 2 }.sum.to_f / (length - 1)
  end

  def stddev
    variance ** 0.5
  end
end

def zscore(val, m, std)
  (val - m) / std
end

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


mutation_counts = nil

contexts_wg = nil

original_ppm = nil
original_pwm = nil
original_bsearch_table = nil

option_parser = OptionParser.new{|opts|
  opts.on('--mutational-context FILE', 'Specify mutational process contexts'){|fn|
    # mutation_counts = val
    mutation_counts = ContextDependentCounts.from_file(fn){|ctx|
      DirectedContext.from_string(ctx)
    }
  }
  opts.on('--mutational-signature FILE:SIGNATURE', 'Specify mutational process file and signature name'){|value|
    cosmic_fn, signature = value.split(':', 2) # 'COSMIC_v3.4_SBS_GRCh38.txt'
    mutation_counts = ContextDependentCounts.from_cosmic_file(cosmic_fn, signature){|ctx|
      DirectedContext.from_string(ctx)
    }
  }
  opts.on('--genomic-context FILE', 'Specify genomic contexts'){|fn|
    contexts_wg = ContextDistribution.from_file(fn).without_unknown
  }
  opts.on('--pfm-motif FILE', 'Specify PFM/PCM file'){|fn|
    original_ppm = PCM.from_file(fn).to_ppm
  }
  opts.on('--pwm-motif FILE', 'Specify PWM file'){|fn|
    original_pwm = PWM.from_file(fn)
  }
  opts.on('--thresholds FILE', 'Specify thresholds grid file'){|fn|
    original_bsearch_table = read_bsearch_table(fn)
  }
}

option_parser.parse!(ARGV)

# raise 'Specify region specific contexts'  unless region_specific_contexts_fn = ARGV[2]
# contexts_rs_frequencies = ContextDistribution.from_file(region_specific_contexts_fn).without_unknown.symmetrized.frequencies


# mutation_counts = ContextDependentCounts.from_cosmic_file('COSMIC_v3.4_SBS_GRCh38.txt', 'SBS1'){|ctx|
#   DirectedContext.from_string(ctx)
# }

mutational_process = MutationProcess.new(mutation_counts, contexts_wg).symmetrized



motifs = Dir.glob('motifs_freeze_best_approved/*').map{|fn| File.readlines(fn).drop(1).map{|l| l.chomp.split("\t").map{|x| Float(x) } } }; nil
# positions = motifs.flatten(1).shuffle; nil

motifs.select!{|matrix| matrix.length == original_ppm.length }


background_infos = 100.times.map {
  $stderr.print '.'
  # pcm = PCM.new(positions.shuffle.first(original_ppm.matrix.length)).to_ppm.to_pcm(100)
  # pcm = PCM.new(motifs.sample.shuffle).to_ppm.to_pcm(100)
  pcm = PCM.new(motifs.sample).to_ppm.to_pcm(100)
  ppm = pcm.to_ppm
  pwm = pcm.to_pwm(pseudocount: :log)
  pwm_fn = write_temp_matrix(pwm)

  thresholds_fn ||= make_thresholds(pwm_fn, thresholds_fn: nil)
  bsearch_table = read_bsearch_table(thresholds_fn)

  mutation_action_infos(mutational_process, ppm, pwm, bsearch_table)
}
puts

attractivenesses = background_infos.map{|info| info[:attractiveness] }
pvalue_changes = background_infos.map{|info| info[:pvalue_change] }
corrected_pvalue_changes = background_infos.map{|info| info[:corrected_pvalue_change] }


attractiveness_mean = 6.561522856702855
attractiveness_std = 2.7399736654803424
pvalue_change_mean = -0.8852969382787378
pvalue_change_std = 0.28796623799385
corrected_pvalue_change_mean = -4.096566980674194
corrected_pvalue_change_std = 1.879946505156352

attractiveness_mean = attractivenesses.mean
attractiveness_std = attractivenesses.stddev

pvalue_change_mean = pvalue_changes.mean
pvalue_change_std = pvalue_changes.stddev

corrected_pvalue_change_mean = corrected_pvalue_changes.mean
corrected_pvalue_change_std = corrected_pvalue_changes.stddev




# motif_name = File.basename(pcm_or_ppm_fn, File.extname(pcm_or_ppm_fn))
# motif_name = 'Motif'

# attractivenesses = []
# pvalue_changes = []
# corrected_pvalue_changes = []


# 1.times do
#   # order = (0 ... original_ppm.matrix.length).to_a.shuffle
#   # ppm = PPM.new(original_ppm.matrix.values_at(*order))
#   # pwm = PWM.new(original_pwm.matrix.values_at(*order))



  
#   # pwm_fn = make_pwm_from_pfm(pfm_fn, pwm_fn: nil)
#   # thresholds_fn ||= make_thresholds(pwm_fn, thresholds_fn: nil)
#   # bsearch_table = read_bsearch_table(thresholds_fn)


#   ## contexts_rs_frequencies = ppm.mean_context_distribution
#   ## puts mutational_process.attractiveness(contexts_rs_frequencies)

#   ## pp DirectedContext.each.map{|dirctx|
#   ##   [dirctx, mutational_process.unnormed_positional_frequency(ppm, 14, dirctx)]
#   ## }.sort_by{|dirctx, val| val }.to_h.transform_values{|v| v.round(3) }

#   ## pp Context.each.map{|ctx|
#   ##   [ctx, ppm.context_frequency_at_pos(ctx, 14)]
#   ## }.sort_by{|ctx, val| val }.to_h.transform_values{|v| v.round(3) }

#   ## pp mutational_process.normalized_densities.sort_by{|dirctx, val| val }.to_h.transform_values{|v| v.round(3) }

#   # positional_profile = mutational_process.positional_profile(ppm)
#   # consensus = ppm.matrix.map{|pos|
#   #   max = pos.max
#   #   pos.zip('ACGT'.chars).select{|v,n| v == max }.map(&:last).join
#   # }
#   # positional_profile.zip(consensus).each_with_index.map{|(val, consensus_nuc), pos|
#   #   puts "#{pos+1}\t#{consensus_nuc}\t#{val.round(3)}"
#   # }
#   # puts '========='

#   # g = Gruff::Line.new
#   # g.data(motif_name, positional_profile)
#   # g.write("#{motif_name}@profile.png")

#   info = mutation_action_infos(mutational_process, ppm, pwm, bsearch_table)
#   info.merge!(motif: motif_name)

#   # pp info
#   # puts [motif, attractiveness.round(2), mean_weight_change.round(2), weight_stddev.round(2)].join("\t")
#   attractivenesses << info[:attractiveness]
#   pvalue_changes << info[:pvalue_change]
#   corrected_pvalue_changes << info[:corrected_pvalue_change]
# end

original_motif_info = mutation_action_infos(mutational_process, original_ppm, original_pwm, original_bsearch_table)

attractivenesses_zscore = zscore(original_motif_info[:attractiveness], attractiveness_mean, attractiveness_std)
pvalue_change_zscore = zscore(original_motif_info[:pvalue_change], pvalue_change_mean, pvalue_change_std)
corrected_pvalue_change_zscore = zscore(original_motif_info[:corrected_pvalue_change], corrected_pvalue_change_mean, corrected_pvalue_change_std)


puts "attractiveness: #{original_motif_info[:attractiveness].round(2)}; z-score: #{attractivenesses_zscore.round(2)}; #{attractiveness_mean.round(2)} ± #{attractiveness_std.round(2)}"
puts "pvalue_change: #{original_motif_info[:pvalue_change].round(2)}; z-score: #{pvalue_change_zscore.round(2)}; #{pvalue_change_mean.round(2)} ± #{pvalue_change_std.round(2)}"
puts "corrected_pvalue_change: #{original_motif_info[:corrected_pvalue_change].round(2)}; z-score: #{corrected_pvalue_change_zscore.round(2)}; #{corrected_pvalue_change_mean.round(2)} ± #{corrected_pvalue_change_std.round(2)}"
