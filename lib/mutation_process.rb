require_relative 'snv_context'

class MutationProcess
  def initialize(mut_rates_by_ctx)
    @mut_rates_by_ctx = mut_rates_by_ctx
  end
  
  def mutation_rate_from_to(alpha, beta, gamma, delta)
    @mut_rates_by_ctx[alpha][beta][gamma][delta]
  end

  def mutation_rate_from_to_any(alpha, beta, gamma)
    @mut_rates_by_ctx[alpha][beta][gamma].each_value.inject(0.0, &:+)
  end

  def self.from_file(filename)
    self.new(normalize_context_counts(read_mutational_profile(filename)))
  end

  def self.read_mutational_profile(filename)
    File.readlines(filename).map{|l|
      ctx, count = l.chomp.split("\t")
      [SNVContext.from_string(ctx), count.to_f]
    }.each_with_object(Hash.new(0)){|(ctx, cnt), result|
      result[ctx] = cnt
    }
  end

  def self.normalize_context_counts(context_counts)
    sum = 2 * context_counts.values.inject(0.0, &:+)
    result = Hash.new{|ha, a|
      ha[a] = Hash.new{|hb, b|
        hb[b] = Hash.new {|hg, g|
          hg[g] = Hash.new(0)
        }
      }
    }
    context_counts.each_with_object(Hash.new(0)){|(ctx, cnt), result|
      result[ctx] += cnt
      result[ctx.revcomp] += cnt
    }.each{|ctx, cnt|
      alpha, beta, gamma, delta = *ctx.values
      result[alpha][beta][gamma][delta] = cnt / sum
    }
    result
  end
end
