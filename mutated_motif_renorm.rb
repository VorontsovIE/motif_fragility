require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/mutation_process'
require_relative 'lib/motif_mutation'
require_relative 'lib/context_distribution'

raise 'Specify PCM file'  unless pcm_fn = ARGV[0]
raise 'Specify mutational context counts'  unless mutational_ctx_fn = ARGV[1]
raise 'Specify genomic contexts'  unless genomic_contexts_fn = ARGV[2]

mutation_contexts = MutationProcess.from_file(mutational_ctx_fn) # already symmetrized
pcm = PCM.from_file(pcm_fn)
ppm = pcm.to_ppm

genomic_contexts = ContextDistribution.from_file(genomic_contexts_fn).without_unknown.symmetrized
genomic_context_frequencies = ContextDistribution.as_nested_indexed_hash(genomic_contexts.frequencies)

motif_context_frequencies = ppm.context_probabilities

renormed_mutation_contexts = mutation_contexts.renormalize_at_reduced_set(genomic_context_frequencies, motif_context_frequencies)

motif_mutation_process = MotifMutation.new(ppm, renormed_mutation_contexts)
mutated_ppm = motif_mutation_process.mutated(prob_total: 1)
mutated_pcm = mutated_ppm.to_pcm(pcm.count)

puts mutated_pcm
