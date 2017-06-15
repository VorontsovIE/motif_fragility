Nucleotides = %w[A C G T]

def numeric?(str)
  !!Float(str) rescue false
end

def complement(str)
  str.tr('ACGTacgt', 'TGCAtgca')
end
def revcomp(str)
  complement(str).reverse
end

def read_matrix(filename)
  File.readlines(filename).map(&:chomp).drop_while{|l|
    ! numeric?( l.split.first )
  }.map{|l|
    l.split.map(&:to_f)
  }
end

def print_matrix(mat)
  mat.each{|pos|
    puts pos.join("\t")
  }
end

def expand_ppm_flanks(ppm, flank_length)
  result = ppm.map(&:dup)
  flank_length.times{
    result.push(Array.new(4, 0.25))
    result.unshift(Array.new(4, 0.25))
  }
  result
end

def drop_flanks(mat, flank_length)
  flank_length == 0 ? mat : mat[flank_length ... (-flank_length)]
end


def add_pseudocount(pcm, pseudocount = 1)
  pcm.map{|counts|    
    counts.map{|cnt| cnt + pseudocount }
  }
end

def pcm_to_ppm(pcm)
  pcm.map{|counts|
    sum = counts.inject(0.0, &:+)
    counts.map{|cnt| cnt / sum }
  }
end

def ppm_to_pcm(ppm, count)
  ppm.map{|probs|
    probs.map{|prob| prob * count }
  }
end


SNVContext = Struct.new(:flank_5, :mid_before, :flank_3, :mid_after) do
  def self.from_string(str)
    flank_5, mid_before, mid_after, flank_3 = str.match(/^([ACGT])\[([ACGT])\/([ACGT])\]([ACGT])$/i).values_at(1,2,3,4)
    self.new(flank_5, mid_before, flank_3, mid_after)
  end

  def before
    "#{flank_5}#{mid_before}#{flank_3}"
  end

  def after
    "#{flank_5}#{mid_after}#{flank_3}"
  end

  def revcomp
    SNVContext.new(complement(flank_3).upcase, complement(mid_before).upcase, complement(flank_5).upcase, complement(mid_after).upcase)
  end

  def to_s
    "#{flank_5}[#{mid_before}/#{mid_after}]#{flank_3}"
  end
  alias_method :inspect, :to_s
end

def read_mutational_profile(filename)
  File.readlines(filename).map{|l|
    ctx, count = l.chomp.split("\t")
    [SNVContext.from_string(ctx), count.to_f]
  }.each_with_object(Hash.new(0)){|(ctx, cnt), result|
    result[ctx] = cnt
  }
end

def normalize_context_counts(context_counts)
  sum = 2 * context_counts.values.inject(0.0, &:+)
  context_counts.each_with_object(Hash.new(0)){|(ctx, cnt), result|
    result[ctx] += cnt
    result[ctx.revcomp] += cnt
  }.each_with_object(Hash.new(0)){|(ctx, cnt), result|
    result[ctx] = cnt / sum
  }
end

def context_probability(ppm, position, ctx_nuc_indices)
  raise 'Position out of range'  unless (position - 1 >= 0) && (position + 1 < ppm.length)
  raise 'Wrong context'  unless ctx_nuc_indices.length == 3
  ((position - 1) .. (position + 1)).zip(ctx_nuc_indices).map{|pos, nuc_at_pos|
    ppm[pos][nuc_at_pos]
  }.inject(1.0, &:*)
end

def summary_context_probability(ppm, ctx_nuc_indices)
  (1 ... (ppm.length - 1)).map{|pos|
    context_probability(ppm, pos, ctx_nuc_indices)
  }.inject(0.0, &:+)
end

def mutation_rate(mut_rates_by_ctx, orig_ctx_indices, final_nuc_idx)
  substitution = SNVContext.new( *Nucleotides.values_at(*orig_ctx_indices, final_nuc_idx) )
  mut_rates_by_ctx[substitution]
end

def summary_mutation_rate(mut_rates_by_ctx, orig_ctx_indices)
  Nucleotides.each_index.map{|ksi|
    mutation_rate(mut_rates_by_ctx, orig_ctx_indices, ksi)
  }.inject(0.0, &:+)
end

# class MutationProcessAtMotif
#   def intialize(ppm, mut_rates_by_ctx)
#     @ppm = ppm
#     @mut_rates_by_ctx = mut_rates_by_ctx
#   end
# end

def mutation_gradient(ppm, mut_rates_by_ctx)
  result = Array.new(ppm.length){ Array.new(4, 0.0) }
  (1 ... (ppm.length - 1)).each do |pos|
    Nucleotides.each_index do |delta|

      result[pos][delta] = begin
        Nucleotides.each_index.map{|alpha|
          Nucleotides.each_index.map{|gamma|
            sum_mut_rate_alpha_delta_gamma = summary_mutation_rate(mut_rates_by_ctx, [alpha, delta, gamma])
            alpha_delta_gamma_at_pos = context_probability(ppm, pos, [alpha, delta, gamma])
            alpha_delta_gamma = summary_context_probability(ppm, [alpha, delta, gamma])
            
            # alpha-beta-gamma for all beta --> alpha-delta-gamma
            incoming = Nucleotides.each_index.map{|beta|
              mut_rate_alpha_beta_gamma_to_delta = mutation_rate(mut_rates_by_ctx, [alpha, beta, gamma], delta)
              alpha_beta_gamma_at_pos = context_probability(ppm, pos, [alpha, beta, gamma])
              alpha_beta_gamma = summary_context_probability(ppm, [alpha, beta, gamma])

              alpha_beta_gamma_at_pos * mut_rate_alpha_beta_gamma_to_delta * (alpha_beta_gamma_at_pos / alpha_beta_gamma)
            }.inject(0.0, &:+)

            # alpha-delta-gamma for all ksi --> alpha-ksi-gamma
            outgoing = alpha_delta_gamma_at_pos * sum_mut_rate_alpha_delta_gamma * (alpha_delta_gamma_at_pos / alpha_delta_gamma)

            incoming - outgoing
            
          }.inject(0.0, &:+)
        }.inject(0.0, &:+)
      end
    end
  end
  result
end

def mutate(ppm, mut_rates_by_ctx, rho)
  ppm_gradient = mutation_gradient(ppm, mut_rates_by_ctx)
  ppm.zip(ppm_gradient).map{|pos, pos_grad|
    pos.zip(pos_grad).map{|el, el_grad|
      el + rho * el_grad
    }
  }
end

raise 'Specify PCM file'  unless pcm_fn = ARGV[0]
raise 'Specify mutational context counts'  unless mutational_ctx_fn = ARGV[1]
raise 'Specify flank expansion length'  unless flank_expansion = ARGV[2]
flank_expansion = Integer(flank_expansion)

pcm = read_matrix(pcm_fn)
ppm = pcm_to_ppm(pcm)
# ppm = pcm_to_ppm(add_pseudocount(read_matrix(pcm_fn), 1))
mutation_contexts = normalize_context_counts(read_mutational_profile(mutational_ctx_fn))
# print_matrix(ppm)
# puts '================'
mutated_ppm = drop_flanks(mutate(expand_ppm_flanks(ppm, flank_expansion), mutation_contexts, 1), flank_expansion)
count = pcm.map{|pos| pos.inject(0.0, &:+) }.max
mutated_pcm = ppm_to_pcm(mutated_ppm, count)
print_matrix(mutated_pcm)
