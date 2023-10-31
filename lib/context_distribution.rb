require_relative 'context'
require_relative 'context_dependent_counts'

class ContextDistribution
  def initialize(counts)
    @counts = counts
    raise 'Empty counts distribution'  unless counts.total_counts > 0
  end

  def self.from_file(filename)
    counts = ContextDependentCounts.from_file(filename, parse_value: ->(str){ str.to_i }){|ctx|
      Context.from_string(ctx)
    }
    self.new(counts)
  end

  def symmetrized
    ContextDistribution.new(@counts.symmetrized)
  end

  def without_unknown
    new_counts = @counts.select{|ctx, cnt| ctx.fully_defined? }.to_h
    ContextDistribution.new(ContextDependentCounts.new(new_counts))
  end

  def self.uniform
    hsh = Context.each.map{|ctx|
      [ctx, 1.0 / (4**3)]
    }.to_h
    self.new(ContextDependentCounts.new(hsh))
  end

  def total_counts; @counts.total_counts; end
  def count(ctx); @counts[ctx]; end
  # def frequency(ctx); @counts.frequencies[ctx]; end
  # def frequencies; @counts.frequencies; end

  def to_s; @counts.to_s; end
  def inspect; "<ContextDistribution: #{@counts}>"; end
end
