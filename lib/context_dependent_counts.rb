class ContextDependentCounts
  attr_reader :counts
  def initialize(counts)
    @counts = counts.dup.tap{|x| x.default = 0 }
  end

  def self.from_file(filename, parse_value: ->(s){ s.to_f }, &block)
    result = Hash.new(0)
    File.readlines(filename).each{|l|
      ctx, count = l.chomp.split("\t")
      ctx = block.call(ctx)  if block_given?
      result[ctx] = parse_value.call(count)
    }
    self.new(result)
  end

  def self.from_cosmic_file(filename, signature_name, parse_value: ->(s){ s.to_f }, &block)
    result = Hash.new(0)
    lns = File.readlines(filename)
    header = lns.first.chomp.split("\t") # Type SBS1 SBS2 ...
    signature_index = header.index(signature_name)
    lns.drop(1).each{|l|
      ctx, count = l.chomp.split("\t").values_at(0, signature_index) # A[C>A]A 0.0008760226353560026
      ctx = block.call(ctx)  if block_given?
      result[ctx] = parse_value.call(count)
    }
    self.new(result)
  end

  def [](ctx)
    counts[ctx]
  end

  def total_counts
    @total_counts ||= counts.values.inject(0, &:+)
  end

  def frequencies
    @frequencies ||= counts.transform_values{|cnt|
      cnt.to_f / total_counts
    }.tap{|x| x.default = 0.0 }.tap{|x| pp x; pp x.values.sum }
  end

  # doubles total number
  def symmetrized
    result = Hash.new(0)
    counts.each{|ctx, cnt|
      result[ctx] += cnt
      result[ctx.revcomp] += cnt
    }
    ContextDependentCounts.new(result)
  end

  def each(&block); @counts.each(&block); end
  include Enumerable

  def to_s
    counts.map{|cnt, ctx|
      [cnt, ctx].join("\t")
    }.join("\n")
  end
  alias_method :inspect, :to_s
end
