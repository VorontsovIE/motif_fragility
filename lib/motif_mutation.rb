require_relative 'ppm'
require_relative 'mutation_process'

class MotifMutation
  attr_reader :ppm, :mut_process
  def initialize(ppm, mut_process)
    @ppm = ppm
    @mut_process = mut_process
  end

  def mutation_gradient_at(pos, delta)
    Nucleotides.each_index.map{|alpha|
      Nucleotides.each_index.map{|gamma|
        sum_mut_rate_alpha_delta_gamma = mut_process.mutation_rate_from_to_any(alpha, delta, gamma)
        alpha_delta_gamma_at_pos = ppm.context_probability_at_pos(alpha, delta, gamma, pos)
        alpha_delta_gamma = ppm.context_probability(alpha, delta, gamma)
        
        # alpha-beta-gamma for all beta --> alpha-delta-gamma
        incoming = Nucleotides.each_index.map{|beta|
          mut_rate_alpha_beta_gamma_to_delta = mut_process.mutation_rate_from_to(alpha, beta, gamma, delta)
          alpha_beta_gamma_at_pos = ppm.context_probability_at_pos(alpha, beta, gamma, pos)
          alpha_beta_gamma = ppm.context_probability(alpha, beta, gamma)

          alpha_beta_gamma_at_pos * mut_rate_alpha_beta_gamma_to_delta * (alpha_beta_gamma_at_pos / alpha_beta_gamma)
        }.inject(0.0, &:+)

        # alpha-delta-gamma for all ksi --> alpha-ksi-gamma
        outgoing = alpha_delta_gamma_at_pos * sum_mut_rate_alpha_delta_gamma * (alpha_delta_gamma_at_pos / alpha_delta_gamma)

        incoming - outgoing
        
      }.inject(0.0, &:+)
    }.inject(0.0, &:+)
  end

  def mutation_gradient
    result = Array.new(ppm.length){ Array.new(4, 0.0) }
    (1 ... (ppm.length - 1)).each do |pos|
      Nucleotides.each_index do |delta|
        result[pos][delta] = mutation_gradient_at(pos, delta)
      end
    end
    result
  end

  # prob_total - probability that any mutation occured
  def mutated(prob_total: 1)
    new_matrix = ppm.matrix.zip(mutation_gradient).map{|pos, pos_grad|
      pos.zip(pos_grad).map{|el, el_grad|
        el + prob_total * el_grad
      }
    }
    PPM.new(new_matrix)
  end
end
