require_relative 'matrix_operations'

class PCM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def count
    matrix.map{|pos| pos.inject(0.0, &:+) }.max
  end

  def to_ppm
    new_matrix = pcm2pfm({matrix: self.matrix, name: nil})[:matrix]
    PPM.new(new_matrix)
  end

  def to_pwm(pseudocount: :log)
    new_matrix = pcm2pwm({matrix: self.matrix, name: nil}, pseudocount: pseudocount)[:matrix]
    PWM.new(new_matrix)
  end

  def self.from_file(filename)
    matrix = read_matrix(filename)[:matrix]
    PCM.new(matrix)
  end

  def to_s
    matrix.map{|pos| pos.join("\t") }.join("\n")
  end
  alias_method :inspect, :to_s

  def with_pseudocount
    new_matrix = matrix.map{|pos|
      pos_count = pos.inject(0.0, &:+)
      actual_pseudocount = Math.log([pos_count, 2].max)
      pos.map{|el| el + 0.25 * actual_pseudocount}
    }
    PCM.new(new_matrix)
  end

  def length; matrix.length; end
  def size; length; end
end

###############

class PWM
  attr_reader :matrix
  def initialize(matrix)
    @matrix = matrix
  end

  def self.from_file(filename)
    matrix = read_matrix(filename)[:matrix]
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

  def length; matrix.length; end
  def size; length; end

  alias_method :inspect, :to_s
end
