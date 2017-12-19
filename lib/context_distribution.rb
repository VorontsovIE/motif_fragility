require_relative 'snv_context'

class ContextDistribution
  attr_reader :count_distribution
  def initialize(count_distribution)
    raise 'Empty contexts count distribution'  unless count_distribution.values.inject(0, &:+) > 0
    upcased_count_distribution = count_distribution.map{|ctx, cnt|
      [ctx.upcase, cnt]
    }.to_h
    upcased_count_distribution.default = 0
    @count_distribution = upcased_count_distribution
  end

  def self.from_file(filename)
    count_distribution = File.readlines(filename).map{|l|
      ctx, cnt = l.chomp.split("\t")
      [ctx, cnt.to_i]
    }.to_h
    self.new(count_distribution)
  end

  def symmetrized
    new_count_distribution = Hash.new(0)
    count_distribution.each{|ctx, cnt|
      ctx_revcomp = ctx.tr('ACGTN', 'TGCAN').reverse
      new_count_distribution[ctx] += cnt
      new_count_distribution[ctx_revcomp] += cnt
    }
    ContextDistribution.new(new_count_distribution)
  end

  def without_unknown
    new_count_distribution = count_distribution.reject{|ctx, cnt|
      ctx.match(/N/)
    }
    ContextDistribution.new(new_count_distribution)
  end

  def frequencies
    @frequencies ||= begin
      sum_cnt = count_distribution.values.inject(0, &:+)
      count_distribution.map{|ctx, cnt|
        [ctx, cnt.to_f / sum_cnt]
      }.to_h
    end
  end

  def to_s
    count_distribution.sort.map{|cnt, ctx|
      [cnt, ctx].join("\t")
    }.join("\n")
  end
  alias_method :inspect, :to_s

  def self.as_nested_indexed_hash(hsh)
    result = Hash.new{|ha, a|
      ha[a] = Hash.new{|hb, b|
        hb[b] = Hash.new(0)
      }
    }
    hsh.each{|ctx, cnt|
      a,b,g = ctx.chars.map{|ch| Nucleotides.index(ch) }
      result[a][b][g] = cnt
    }
    result
  end
end
