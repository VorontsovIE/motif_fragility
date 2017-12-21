MOTIF_PCM=$1
find mutational_signatures/ -xtype f \
  | xargs -n1 -I{} basename -s '.tsv' "{}" | sort \
  | xargs -n1 -I{} echo \
      "echo -n {} $'\t' ;" \
      "java -cp ape.jar ru.autosome.macroape.EvalSimilarity $MOTIF_PCM " \
                        " <( ruby mutate_motif.rb $MOTIF_PCM 'mutational_signatures/{}.tsv' genomic_contexts/hg38_contexts.tsv --flank-length 1 ) " \
                        " --first-pcm --second-pcm --position 0,direct -d 100 " \
      " | grep -Pe '^S\b' | cut -f2" \
  | bash

# batch run:
#   find ~/iogen/motif_collections/motifs/pcm/hocomoco/ -iname '*.pcm' | sort | xargs -n1 -I{} echo "echo '---------'; basename -s '.pcm' {}; ./check_motif_different_cancers.sh {} 11 | ruby -e 'puts readlines.sort_by{|l| l.split(\"\t\").last.to_f}'; " | bash | tee all_motifs.txt

