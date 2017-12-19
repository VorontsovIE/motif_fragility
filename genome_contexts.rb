NUCS = %w[A C G T N]
NucToIdx = {'A'=>0, 'C'=>1, 'G'=>2, 'T'=>3, 'N'=>4, 'a'=>0,'c'=>1,'g'=>2,'t'=>3,'n'=>4}

class ChromosomesFromFasta
  def initialize(filename, chrom_names)
    @filename = filename
    @chrom_names = chrom_names
  end

  def each(&block)
    return enum_for(:each)  unless block_given?
    chr = nil
    seqs = []
    File.open(@filename) do |f|
      f.each_line{|l|
        if l.start_with? '>'
          yield chr, seqs.join  if @chrom_names.include?(chr)
          seqs = []
          chr = l[1..-1].strip
        else
          seqs << l.strip
        end
      }
      yield chr, seqs.join  if @chrom_names.include?(chr)
    end
  end
end

class ChromosomesFromPlainFiles
  def initialize(folder, chrom_names)
    @folder = folder
    @chrom_names = chrom_names
  end

  def each(&block)
    return enum_for(:each)  unless block_given?
    @chrom_names.each{|chr|
      yield chr, File.read(File.join(@folder, "#{chr}.plain"))
    }
  end
end


hsh = Hash.new(0)
human_chromosomes = (1..22).map(&:to_s) + ['X', 'Y']
mouse_chromosomes = (1..19).map(&:to_s) + ['X', 'Y']
human_chromosome_names = human_chromosomes.map{|chr_num| "chr#{chr_num}"}
mouse_chromosome_names = mouse_chromosomes.map{|chr_num| "chr#{chr_num}"}
# chromosome_reader = ChromosomesFromPlainFiles.new('/home/ilya/iogen/genome/mm9/', mouse_chromosome_names)
# chromosome_reader = ChromosomesFromPlainFiles.new('/home/ilya/iogen/genome/mm10/', mouse_chromosome_names)
chromosome_reader = ChromosomesFromPlainFiles.new('/home/ilya/iogen/genome/hg19/', human_chromosome_names)
# chromosome_reader = ChromosomesFromFasta.new('/home/ilya/genome/human/hg38.fa', human_chromosome_names)
chromosome_reader.each do |chr, seq|
# chromosomes.each do |chr|
  $stderr.puts(chr)
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
