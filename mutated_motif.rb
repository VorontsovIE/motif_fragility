require_relative 'lib/pcm'
require_relative 'lib/ppm'
require_relative 'lib/mutation_process'
require_relative 'lib/motif_mutation'

Nucleotides = %w[A C G T]

raise 'Specify PCM file'  unless pcm_fn = ARGV[0]
raise 'Specify mutational context counts'  unless mutational_ctx_fn = ARGV[1]
raise 'Specify flank expansion length'  unless flank_expansion = ARGV[2]
flank_expansion = Integer(flank_expansion)

mutation_contexts = MutationProcess.from_file(mutational_ctx_fn)

pcm = PCM.from_file(pcm_fn)
ppm = pcm.to_ppm

ppm_expanded = ppm.expand_flanks(flank_expansion)
motif_mutation_process = MotifMutation.new(ppm_expanded, mutation_contexts)
mutated_expanded_ppm = motif_mutation_process.mutated(prob_total: 1)
mutated_ppm = mutated_expanded_ppm.drop_flanks(flank_expansion)
mutated_pcm = mutated_ppm.to_pcm(pcm.count)

puts mutated_pcm
