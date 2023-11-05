require_relative 'lib/fasta_reader'
require_relative 'lib/context'

raise 'Specify FASTA file with sequences' unless fasta_fn = ARGV[0]
raise 'Specify flank length' unless flank_length = Integer(ARGV[1])

raise 'Flank should exist' if flank_length < 1

context_counts = Hash.new(0)
seq_reader = FastaReader.new(fasta_fn)
seq_reader.each do |hdr, seq|
  raise 'Flank length incompatible with sequence size'  unless seq.length == 2*flank_length + 1
  match = hdr.match(/@\[(?<ref>[ACGT])\/(?<alt>[ACGT])\]::/)
  raise "Reference nucleotide in substitution doesn't match genomic sequence"  unless seq[flank_length].upcase == match[:ref].upcase
  dirctx_str = "#{seq[flank_length - 1]}[#{match[:ref]}/#{match[:alt]}]#{seq[flank_length + 1]}".upcase
  next  if dirctx_str.match /N/
  dirctx = DirectedContext.from_string(dirctx_str)

  context_counts[dirctx] += 1
end

DirectedContext.each do |dirctx|
  puts [dirctx, context_counts[dirctx]].join("\t")
end
