def numeric?(str)
  !!Float(str) rescue false
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
    matrix = File.readlines(filename).map(&:chomp).drop_while{|l|
      ! numeric?( l.split.first )
    }.map{|l|
      l.split.map(&:to_f)
    }
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
    matrix = File.readlines(filename).map(&:chomp).drop_while{|l|
      ! numeric?( l.split.first )
    }.map{|l|
      l.split.map(&:to_f)
    }
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

  alias_method :inspect, :to_s
end
