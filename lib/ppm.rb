class PPM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def length
    matrix.length
  end

  def context_probability_at_pos(alpha, beta, gamma, position)
    @context_probability_at_pos_cache ||= Hash.new{|ha, a|
      ha[a] = Hash.new{|hb, b|
        hb[b] = Hash.new{|hg, g|
          hg[g] = Hash.new{|hpos, pos|
            hpos[pos] = context_probability_at_pos_calculate(a, b, g, pos)
          }
        }
      }
    }
    @context_probability_at_pos_cache[alpha][beta][gamma][position]
  end

  def context_probability_at_pos_calculate(alpha, beta, gamma, position)
    raise 'Position out of range'  unless (position - 1 >= 0) && (position + 1 < length)
    matrix[position - 1][alpha] * matrix[position][beta] * matrix[position + 1][gamma]
  end

  def context_probability(alpha, beta, gamma)
    @context_probability_cache ||= Hash.new{|ha, a|
      ha[a] = Hash.new{|hb, b|
        hb[b] = Hash.new{|hg, g|
          hg[g] = context_probability_calculate(a, b, g)
        }
      }
    }
    @context_probability_cache[alpha][beta][gamma]
  end

  def context_probability_calculate(alpha, beta, gamma)
    (1 ... (length - 1)).map{|pos|
      context_probability_at_pos(alpha, beta, gamma, pos)
    }.inject(0.0, &:+)
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
