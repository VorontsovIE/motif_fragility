Nucleotides = %w[A C G T]

def complement(nuc_idx)
  3 - nuc_idx
end

SNVContext = Struct.new(:flank_5, :mid_before, :flank_3, :mid_after) do
  def self.from_string(str)
    flank_5, mid_before, mid_after, flank_3 = str.upcase.match(/^([ACGT])\[([ACGT])\/([ACGT])\]([ACGT])$/i).values_at(1,2,3,4).map{|nuc| Nucleotides.index(nuc) }
    self.new(flank_5, mid_before, flank_3, mid_after)
  end

  def revcomp
    SNVContext.new(complement(flank_3), complement(mid_before), complement(flank_5), complement(mid_after))
  end

  def to_s
    "#{Nucleotides[flank_5]}[#{Nucleotides[mid_before]}/#{Nucleotides[mid_after]}]#{Nucleotides[flank_3]}"
  end
  alias_method :inspect, :to_s
end
