class PPM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def length
    matrix.length
  end

  def context_probability_at_pos(ctx_nuc_indices, position)
    raise 'Position out of range'  unless (position - 1 >= 0) && (position + 1 < length)
    raise 'Wrong context'  unless ctx_nuc_indices.length == 3
    alpha, beta, gamma = *ctx_nuc_indices
    matrix[position - 1][alpha] * matrix[position][beta] * matrix[position + 1][gamma]
  end

  def context_probability(ctx_nuc_indices)
    @context_probability_cache ||= {}
    @context_probability_cache[ctx_nuc_indices] ||= begin
      (1 ... (length - 1)).map{|pos|
        context_probability_at_pos(ctx_nuc_indices, pos)
      }.inject(0.0, &:+)
    end
    @context_probability_cache[ctx_nuc_indices]
  end

  def expand_flanks(flank_length)
    new_matrix = matrix.map(&:dup)
    flank_length.times{
      new_matrix.push(Array.new(4, 0.25))
      new_matrix.unshift(Array.new(4, 0.25))
    }
    PPM.new(new_matrix)
  end

  def drop_flanks(flank_length)
    if flank_length == 0
      self
    else
      new_matrix = matrix[flank_length ... (-flank_length)]
      PPM.new(new_matrix)
    end
  end

  def to_pcm(count)
    new_matrix = matrix.map{|probs|
      probs.map{|prob| prob * count }
    }
    PCM.new(new_matrix)
  end

  def to_s
    matrix.map{|pos| pos.join("\t") }.join("\n")
  end
  alias_method :inspect, :to_s
end
