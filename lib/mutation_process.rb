require_relative 'context'
require_relative 'context_dependent_counts'

class MutationProcess
  def initialize(mutation_counts, seq_context)
    @mutation_counts = mutation_counts
    @seq_context = seq_context
    raise 'Empty counts distribution'  unless mutation_counts.total_counts > 0
  end

  def symmetrized
    MutationProcess.new(@mutation_counts.symmetrized, @seq_context.symmetrized)
  end

  def rate_from_to(directed_context)
    @mutation_counts.frequencies[directed_context]
  end

  def rate_from_to_any(context)
    context.each_directed.map{|directed_context| 
      @mutation_counts.frequencies[directed_context]
    }.inject(0.0, &:+)
  end

  def direction_fraction(directed_context)
    original_context = directed_context.original_context
    sum = rate_from_to_any(original_context)
    sum == 0 ? 0.25 : rate_from_to(directed_context).to_f / sum
  end

  # μ(ctx->ctx') = m(ctx->ctx') / N(ctx)
  def densities
    DirectedContext.each.map{|directed_ctx, cnt|
      orig_ctx = directed_ctx.original_context
      density = @mutation_counts[directed_ctx].to_f / @seq_context.count(orig_ctx)
      [directed_ctx, density]
    }.to_h
  end

  # μ(ctx->ctx') / (m/N)
  def normalized_densities
    @normalized_densities ||= begin
      context_independent_density = @mutation_counts.total_counts / @seq_context.total_counts
      densities.transform_values{|density|
        density / context_independent_density
      }
    end
  end

  def attractiveness(contexts_rs_frequencies)
    DirectedContext.each.map{|directed_ctx|
      orig_ctx = directed_ctx.original_context
      normalized_densities[directed_ctx] * contexts_rs_frequencies[orig_ctx]
    }.inject(0.0, &:+)
  end
  
  def motif_attractiveness(ppm)
    rs_frequencies = ppm.mean_context_distribution
    attractiveness(rs_frequencies)
  end

  def unnormed_positional_frequency(ppm, pos, directed_ctx)
    orig_ctx = directed_ctx.original_context
    normalized_densities[directed_ctx] * ppm.context_frequency_at_pos(orig_ctx, pos)
  end

  def kappa(ppm)
    inv_kappa = (0 ... ppm.length).map{|pos|
      DirectedContext.each.map{|directed_ctx|
        orig_ctx = directed_ctx.original_context
        unnormed_positional_frequency(ppm, pos, directed_ctx)
      }.sum
    }.sum
    1.0 / inv_kappa
  end

  def positional_profile(ppm)
    ppm_kappa = kappa(ppm)
    non_normalized_context_profile = (0 ... ppm.length).map{|pos|
      DirectedContext.each.map{|directed_ctx|
        unnormed_positional_frequency(ppm, pos, directed_ctx)
      }.sum * ppm_kappa
    }
  end

  def mean_weight_change(ppm, pwm)
    ppm_kappa = kappa(ppm)
    non_normalized_context_profile = (0 ... ppm.length).map{|pos|
      val = DirectedContext.each.map{|directed_ctx|
        unnormed_positional_frequency(ppm, pos, directed_ctx) * pwm.weight_delta(pos, directed_ctx.beta, directed_ctx.delta)
      }.sum
    }.sum * ppm_kappa
  end

  def weight_stddev(ppm, pwm)
    ppm_kappa = kappa(ppm)
    sqr_change = non_normalized_context_profile = (0 ... ppm.length).map{|pos|
      val = DirectedContext.each.map{|directed_ctx|
        unnormed_positional_frequency(ppm, pos, directed_ctx) * pwm.weight_delta(pos, directed_ctx.beta, directed_ctx.delta) ** 2
      }.sum
    }.sum * ppm_kappa
    sqr_change ** 0.5
  end
end
