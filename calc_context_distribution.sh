# input file sample: ~/iogen/projects/somatic_mutations/results/AllSNVs/Alexandrov/ALL/cancer.txt 
# run all:
#   find /home/ilya/iogen/projects/somatic_mutations/results/AllSNVs/Alexandrov/ -xtype d | tail -n+2 | xargs -I{} -n1 basename "{}" | xargs -n1 -I{} echo './calc_context_distribution.sh "/home/ilya/iogen/projects/somatic_mutations/results/AllSNVs/Alexandrov/{}/cancer.txt" > "context_distributions/{}.tsv"' | bash
SNV_LIST="$1"
cat "$SNV_LIST" | tail -n+2 | cut -f1 | ruby -e 'ARGF.each_line{|l| puts l.split("@").last }' | sort | uniq -c | ruby -e 'puts readlines.map{|l| l.strip.split.map(&:strip).reverse.join("\t") }'
