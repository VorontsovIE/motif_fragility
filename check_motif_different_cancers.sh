MOTIF_PCM=$1
EXPANSION_LENGTH=$2
find context_distributions/ -xtype f | xargs -n1 -I{} basename -s '.tsv' "{}" | sort | xargs -n1 -I{} echo "echo -n {} $'\t' ; java -cp ~/iogen/ape-2.0.3.jar ru.autosome.macroape.EvalSimilarity $MOTIF_PCM <( ruby mutated_motif.rb $MOTIF_PCM 'context_distributions/{}.tsv' $EXPANSION_LENGTH ) --first-pcm --second-pcm --position 0,direct | grep -Pe '^S\b' | cut -f2" | bash
