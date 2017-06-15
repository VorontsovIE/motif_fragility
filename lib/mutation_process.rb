require_relative 'snv_context'

class MutationProcess
  def initialize(mut_rates_by_ctx)
    @mut_rates_by_ctx = mut_rates_by_ctx
  end
  
  def mutation_rate_from_to(orig_ctx_indices, final_nuc_idx)
    substitution = SNVContext.new(*orig_ctx_indices, final_nuc_idx)
    @mut_rates_by_ctx[substitution]
  end

  def mutation_rate_from_to_any(orig_ctx_indices)
    Nucleotides.each_index.map{|ksi|
      mutation_rate_from_to(orig_ctx_indices, ksi)
    }.inject(0.0, &:+)
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
    context_counts.each_with_object(Hash.new(0)){|(ctx, cnt), result|
      result[ctx] += cnt
      result[ctx.revcomp] += cnt
    }.each_with_object(Hash.new(0)){|(ctx, cnt), result|
      result[ctx] = cnt / sum
    }
  end
end
