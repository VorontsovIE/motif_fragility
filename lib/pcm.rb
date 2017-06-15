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
end
