
module Enumerable
  def mean
    (length == 0) ? nil : sum.to_f / length
  end

  def variance
    m = mean
    (length < 2) ? nil : map{|x| (x - m) ** 2 }.sum.to_f / (length - 1)
  end

  def stddev
    variance ** 0.5
  end

  def quantile(rate)
    sort[ (length * rate).to_i ]
  end
end

def zscore(val, m, std)
  (val - m) / std
end
