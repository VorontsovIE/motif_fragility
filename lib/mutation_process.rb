require_relative 'snv_context'

class MutationProcess
  attr_reader :mut_rates_by_ctx
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

  # reduce original distribution on whole genome to a motif-specific ensemble of site with different context distribution
  def renormalize_at_reduced_set(total_ctx_freqs, reduced_ctx_freqs)
    # total_ctx_freqs = ContextDistribution.as_nested_indexed_hash( total_ctx_freqs.frequencies )
    # reduced_ctx_freqs = ContextDistribution.as_nested_indexed_hash( reduced_ctx_freqs.frequencies )

    unnormed_result = MutationProcess.empty_mutation_in_context_hsh
    sum = 0.0
    ((0...4).to_a).repeated_permutation(4){|a,b,g,d|
      unnormed_result[a][b][g][d] = mutation_rate_from_to(a, b, g, d) * reduced_ctx_freqs[a][b][g] / total_ctx_freqs[a][b][g]
      sum += unnormed_result[a][b][g][d]
    }

    result = MutationProcess.empty_mutation_in_context_hsh
    ((0...4).to_a).repeated_permutation(4){|a,b,g,d|
      result[a][b][g][d] = unnormed_result[a][b][g][d] / sum
    }
    MutationProcess.new(result)
  end

  def self.empty_mutation_in_context_hsh
    Hash.new{|ha, a|
      ha[a] = Hash.new{|hb, b|
        hb[b] = Hash.new {|hg, g|
          hg[g] = Hash.new(0)
        }
      }
    }
  end

  def self.normalize_context_counts(context_counts)
    sum = 2 * context_counts.values.inject(0.0, &:+)
    result = empty_mutation_in_context_hsh
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
