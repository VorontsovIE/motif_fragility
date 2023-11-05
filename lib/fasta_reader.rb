class FastaReader
  def initialize(filename)
    @filename = filename
  end

  def each(&block)
    return enum_for(:each)  unless block_given?
    header = nil
    seqs = []
    File.open(@filename) do |f|
      f.each_line{|l|
        l.strip!
        if l.start_with?('>')
          yield header, (seqs.size == 1 ? seqs.first : seqs.join)  if header
          seqs = []
          header = l[1..-1].lstrip
        else
          seqs << l
        end
      }
      yield header, (seqs.size == 1 ? seqs.first : seqs.join)
    end
  end
end
