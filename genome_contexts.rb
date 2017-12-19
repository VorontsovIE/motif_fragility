NUCS = %w[A C G T N]
NucToIdx = {'A'=>0, 'C'=>1, 'G'=>2, 'T'=>3, 'N'=>4, 'a'=>0,'c'=>1,'g'=>2,'t'=>3,'n'=>4}

hsh = Hash.new(0)
# chromosomes = ((1..22).to_a + ['X', 'Y']).map{|chr_num| "chr#{chr_num}"} ## human
chromosomes = ((1..19).to_a + ['X', 'Y']).map{|chr_num| "chr#{chr_num}"} ## mouse
chromosomes.each do |chr|
  $stderr.puts(chr)
  # seq = File.read("/home/ilya/iogen/genome/hg19/#{chr}.plain") ## human
  seq = File.read("/home/ilya/iogen/genome/mm10/#{chr}.plain")  ## mouse
  seq_idx = seq.each_char.lazy.map{|ch| NucToIdx[ch] }
  idx = 0

  n1,n2,n3 = seq_idx.first(3)
  idx = n1 * 25 + n2 * 5 + n3
  hsh[idx] += 1

  seq_idx.drop(3).each{|idx_add|
    idx = (5 * idx + idx_add) % 125
    hsh[idx] += 1
  }
end

res = NUCS.repeated_permutation(3).map{|a,b,c|
  ai,bi,ci = [a,b,c].map{|ch|
    NUCS.index(ch)
  }
  [[a,b,c].join, hsh[25*ai+5*bi+ci]]
}.to_h

res.sort.each{|ctx, cnt|
  puts [ctx, cnt].join("\t")
}
