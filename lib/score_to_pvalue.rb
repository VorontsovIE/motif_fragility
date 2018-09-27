def read_bsearch_table(filename)
  File.readlines(filename).map{|l| l.chomp.split("\t").map(&:to_f) }
end

def pvalue_by_score(requested_score, bsearch_table)
  (bsearch_table.bsearch{|score, pvalue| score >= requested_score } || bsearch_table.last).last
end
