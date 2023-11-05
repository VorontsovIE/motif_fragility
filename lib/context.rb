NUCLEOTIDES = %w[A C G T N]

def complement_idx(nuc_idx)
  (nuc_idx < 4) ? 3 - nuc_idx : 4
end

# (alpha, beta, gamma)
class Context
  def initialize(alpha, beta, gamma)
    raise 'Only indices for ACGTN allowed'  unless [alpha, beta, gamma].all?{|idx|  0 <= idx && idx <= 4}
    @alpha, @beta, @gamma = alpha, beta, gamma
  end

  def indices
    [@alpha, @beta, @gamma]
  end

  def nucleotides
    [@alpha, @beta, @gamma].map{|idx| NUCLEOTIDES[idx] }
  end

  def to_s
    nucleotides.join
  end
  alias_method :inspect, :to_s

  include Comparable
  def <=>(other)
    self.to_s <=> other.to_s
  end

  def eql?(other); indices == other.indices; end
  def hash; indices.hash; end

  def revcomp
    Context.new(complement_idx(@gamma), complement_idx(@beta), complement_idx(@alpha))
  end

  def fully_defined?
    [@alpha, @beta, @gamma].all?{|idx| idx < 4 }
  end

  def directed(delta)
    DirectedContext.new(@alpha, @beta, @gamma, delta)
  end

  def each_directed
    return enum_for(:each_directed)  unless block_given?
    (0..3).each{|delta|
      yield directed(delta)  if delta != @beta
    }
  end

  def self.from_string(str)
    alpha, beta, gamma = str.upcase.each_char.map{|letter| NUCLEOTIDES.index(letter) }
    self.new(alpha, beta, gamma)
  end

  def self.each(with_n: false)
    return enum_for(:each, with_n: with_n)  unless block_given?
    nuc_indices = [0,1,2,3]
    nuc_indices << 4  if with_n
    nuc_indices.repeated_permutation(3).each{|alpha, beta, gamma|
      yield Context.new(alpha, beta, gamma)
    }
  end
end

# (alpha, beta, gamma) --> (alpha, delta, gamma)
class DirectedContext
  attr_reader :alpha, :beta, :gamma, :delta
  SNV_REGEXP = /^(?<flank_5>[ACGT])\[(?<before>[ACGT])[\/>](?<after>[ACGT])\](?<flank_3>[ACGT])$/i
  def initialize(alpha, beta, gamma, delta)
    raise 'Only indices for ACGT allowed'  unless [alpha, beta, gamma, delta].all?{|idx|  0 <= idx && idx <= 3}
    raise 'Central nucleotide not altered'  if beta == delta
    @alpha, @beta, @gamma, @delta = alpha, beta, gamma, delta
  end

  def indices
    [@alpha, @beta, @gamma, @delta]
  end

  def nucleotides
    [@alpha, @beta, @gamma, @delta].map{|idx| NUCLEOTIDES[idx] }
  end

  def to_s
    alpha, beta, gamma, delta = *nucleotides
    "#{alpha}[#{beta}/#{delta}]#{gamma}"
  end
  alias_method :inspect, :to_s

  include Comparable
  def <=>(other)
    self.to_s <=> other.to_s
  end

  def eql?(other); indices == other.indices; end
  def hash; indices.hash; end

  def revcomp
    DirectedContext.new(complement_idx(@gamma), complement_idx(@beta), complement_idx(@alpha), complement_idx(@delta))
  end

  def fully_defined?
    [@alpha, @beta, @gamma, @delta].all?{|idx| idx < 4 }
  end

  def original_context
    Context.new(@alpha, @beta, @gamma)
  end

  def mutated_context
    Context.new(@alpha, @delta, @gamma)
  end

  def self.from_string(str)
    match = str.upcase.match(SNV_REGEXP)
    alpha, beta, gamma, delta = match.values_at(:flank_5, :before, :flank_3, :after).map{|nuc|
      NUCLEOTIDES.index(nuc)
    }
    self.new(alpha, beta, gamma, delta)
  end

  def self.each
    return enum_for(:each)  unless block_given?
    nuc_indices = [0,1,2,3]
    nuc_indices.repeated_permutation(4).each{|alpha, beta, gamma, delta|
      yield DirectedContext.new(alpha, beta, gamma, delta)  if delta != beta
    }
  end
end
