def numeric?(str)
  !!Float(str) rescue false
end

def read_matrix(filename)
  File.readlines(filename).map(&:chomp).map(&:split).drop_while{|row|
    ! numeric?( row.first )
  }.map{|row|
    row.map(&:to_f)
  }
end

class PCM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def count
    matrix.map{|pos| pos.inject(0.0, &:+) }.max
  end

  def to_ppm
    new_matrix = matrix.map{|counts|
      sum = counts.inject(0.0, &:+)
      counts.map{|cnt| cnt / sum }
    }
    PPM.new(new_matrix)
  end

  def to_pwm
    new_matrix = matrix.map do |pos|
      count = pos.inject(0.0, &:+)
      actual_pseudocount = Math.log([count, 2].max)
      pos.map do |el|
        Math.log((el + 0.25 * actual_pseudocount).to_f / (0.25*(count + actual_pseudocount)) )
      end
    end
    PWM.new(new_matrix)
  end

  def self.from_file(filename)
    matrix = read_matrix(filename)
    PCM.new(matrix)
  end

  def to_s
    matrix.map{|pos| pos.join("\t") }.join("\n")
  end
  alias_method :inspect, :to_s

  def with_pseudocount
    new_matrix = matrix.map{|pos|
      count = pos.inject(0.0, &:+)
      actual_pseudocount = Math.log([count, 2].max)
      pos.map{|el| el + 0.25 * actual_pseudocount}
    }
    PCM.new(new_matrix)
  end
end

###############

class PWM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def self.from_file(filename)
    matrix = read_matrix(filename)
    PWM.new(matrix)
  end

  def to_s
    matrix.map{|pos| pos.join("\t") }.join("\n")
  end

  def mean_score(ppm)
    matrix.zip(ppm.matrix).map{|pwm_pos, ppm_pos|
      pwm_pos.zip(ppm_pos).map{|pwm_el, ppm_el|
        pwm_el * ppm_el
      }.inject(0.0, &:+)
    }.inject(0.0, &:+)
  end

  def weight_delta(pos, before_idx, after_idx)
    matrix[pos][after_idx] - matrix[pos][before_idx]
  end

  def expand_flanks(flank_length)
    new_matrix = matrix.map(&:dup)
    flank_length.times{
      new_matrix.push(Array.new(4, 0.0))
      new_matrix.unshift(Array.new(4, 0.0))
    }
    PWM.new(new_matrix)
  end


  alias_method :inspect, :to_s
end
