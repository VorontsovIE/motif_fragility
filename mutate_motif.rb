require 'optparse'
require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/mutation_process'
require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'

report_site_exposure = false
flank_expansion = 1
OptionParser.new{|opts|
  opts.banner = "Usage:\n\truby mutated_motif_renorm.rb  <motif PCM>  <mutational signature>  <genomic context distribution>  [options]"
  opts.separator 'Options:'
  opts.on('--report-site-exposure', "Report site exposure to a mutational process (instead of mutated motif)") {
    report_site_exposure = true
  }
  opts.on('--flank-length LENGTH', "Motif flanks length (default: 1nt so that boundary nucleotides could be mutated)") {|val|
    flank_expansion = Integer(val)
  }
}.parse!(ARGV)

raise 'Specify PCM file'  unless pcm_fn = ARGV[0]
raise 'Specify mutational context counts'  unless mutational_ctx_fn = ARGV[1]
raise 'Specify genomic contexts'  unless genomic_contexts_fn = ARGV[2]

mutation_contexts = MutationProcess.from_file(mutational_ctx_fn) # already symmetrized
pcm = PCM.from_file(pcm_fn)
ppm = pcm.to_ppm
ppm_expanded = ppm.expand_flanks(flank_expansion)

genomic_contexts = ContextDistribution.from_file(genomic_contexts_fn).without_unknown.symmetrized
genomic_context_frequencies = ContextDistribution.as_nested_indexed_hash(genomic_contexts.frequencies)

motif_context_frequencies = ppm_expanded.mean_context_probabilities

if report_site_exposure
  puts mutation_contexts.site_exposure(genomic_context_frequencies, motif_context_frequencies)
  exit
end

renormed_mutation_contexts = mutation_contexts.renormalize_at_reduced_set(genomic_context_frequencies, motif_context_frequencies)

motif_mutation_process = MotifMutation.new(ppm_expanded, renormed_mutation_contexts)
mutated_expanded_ppm = motif_mutation_process.mutated(prob_total: 1)
mutated_ppm = mutated_expanded_ppm.drop_flanks(flank_expansion)
mutated_pcm = mutated_ppm.to_pcm(pcm.count)

puts mutated_pcm
