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

  # probability to hit site divided by a fraction of genome occupied by sites ($\varkappa = P_0 / J$ in terms of an old paper; Q in new terms)
  def site_exposure(total_ctx_freqs, reduced_ctx_freqs)
    ((0...4).to_a).repeated_permutation(3).map{|a,b,g|
      mutation_rate_from_to_any(a, b, g) * reduced_ctx_freqs[a][b][g] / total_ctx_freqs[a][b][g]
    }.inject(0.0, &:+)
  end

  # reduce original distribution on whole genome to a motif-specific ensemble of site with different context distribution
  def renormalize_at_reduced_set(total_ctx_freqs, reduced_ctx_freqs)
    norm = 1.0 / site_exposure(total_ctx_freqs, reduced_ctx_freqs)
    result = MutationProcess.empty_mutation_in_context_hsh
    ((0...4).to_a).repeated_permutation(4){|a,b,g,d|
      result[a][b][g][d] = norm * mutation_rate_from_to(a, b, g, d) * reduced_ctx_freqs[a][b][g] / total_ctx_freqs[a][b][g]
    }
    MutationProcess.new(result)
  end

  def self.empty_mutation_in_context_hsh
    (0...4).each_with_object(Hash.new){|a, h1|
      h1[a] = (0...4).each_with_object(Hash.new){|b, h2|
        h2[b] = (0...4).each_with_object(Hash.new){|g, h3|
          h3[g] = (0...4).each_with_object(Hash.new){|d, h4|
            h4[d] = 0
          }
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
