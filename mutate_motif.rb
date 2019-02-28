require 'optparse'
require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/mutation_process'
require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'
require_relative 'lib/score_to_pvalue'

report_site_exposure = false
report_mean_score_change = false
pwm_filename = nil
pvalues_filename = nil
flank_expansion = 1
pseudocount = false

OptionParser.new{|opts|
  opts.banner = "Usage:\n\truby mutated_motif_renorm.rb  <motif PCM>  <mutational signature>  <genomic context distribution>  [options]"
  opts.separator 'Options:'
  opts.on('--report-site-exposure', "Report site exposure to a mutational process (instead of mutated motif)") {
    report_site_exposure = true
  }
  opts.on('--report-mean-score-change [PWM_FILE]', "Report mean PWM score of motif given original word distribution and word distribution of mutated word set") {|val|
    report_mean_score_change = true
    pwm_filename = val
  }
  opts.on('--add-pseudocount', "add pseudocount (logarithmic) to PCM") {
    pseudocount = true
  }
  opts.on('--flank-length LENGTH', "Motif flanks length (default: 1nt so that boundary nucleotides could be mutated)") {|val|
    flank_expansion = Integer(val)
  }

  opts.on('--report-as-pvalues FILE', 'specify PWM score <--> P-value conversion') {|val|
    pvalues_filename = val
  }
}.parse!(ARGV)

raise 'Specify PCM file'  unless pcm_fn = ARGV[0]
raise 'Specify mutational context counts'  unless mutational_ctx_fn = ARGV[1]
raise 'Specify genomic contexts'  unless genomic_contexts_fn = ARGV[2]

raise 'Either site exposure of mean score change can be reported, not both'  if report_site_exposure && report_mean_score_change

mutation_counts = MutationProcess.from_file(mutational_ctx_fn) # already symmetrized
pcm = PCM.from_file(pcm_fn)
count = pcm.count  # we store count before adding pseudocount

pcm = pcm.with_pseudocount  if pseudocount
ppm = pcm.to_ppm
ppm_expanded = ppm.expand_flanks(flank_expansion)

if genomic_contexts_fn.downcase == 'uniform'
  genomic_contexts = ContextDistribution.uniform
else
  genomic_contexts = ContextDistribution.from_file(genomic_contexts_fn).without_unknown.symmetrized
end
mutation_rates = mutation_counts.normalized_by(genomic_contexts)

motif_context_frequencies = ppm_expanded.mean_context_probabilities




site_exposure = mutation_contexts.site_exposure(genomic_context_frequencies, motif_context_frequencies)
if report_site_exposure
  puts site_exposure
  exit
end

# renormed_mutation_contexts = mutation_contexts.renormalize_at_reduced_set(genomic_context_frequencies, motif_context_frequencies)

motif_mutation_process = MotifMutation.new(ppm_expanded, renormed_mutation_contexts)

# we calculate fragility related to whole-genome process, not to site-specific. 
# So we put prob_total not equal to 1 but to site attractiveness (or site exposure) 
mutated_expanded_ppm = motif_mutation_process.mutated(prob_total: site_exposure)

mutated_ppm = mutated_expanded_ppm.drop_flanks(flank_expansion)
if report_mean_score_change
  if pwm_filename
    pwm = PWM.from_file(pwm_filename)
  else
    pwm = pcm.to_pwm
  end
  score_1 = pwm.mean_score(ppm)
  score_2 = pwm.mean_score(mutated_ppm)
  if pvalues_filename
    threshold_pvalue_table = read_bsearch_table(pvalues_filename)
    pvalue_1 = pvalue_by_score(score_1, threshold_pvalue_table)
    pvalue_2 = pvalue_by_score(score_2, threshold_pvalue_table)
    puts [pvalue_1, pvalue_2, '-->', Math.log10(pvalue_1 / pvalue_2)].join("\t")
  else
    puts [score_1, score_2].join("\t")
  end
  exit
end

mutated_pcm = mutated_ppm.to_pcm(count)

puts mutated_pcm
