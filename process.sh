find /home/ilya/iogen/projects/hocomoco11/final_bundle/hocomoco11/full/HUMAN/mono/pcm -xtype f \
  | xargs -n1 -I{} basename -s .pcm '{}' \
  | sort \
  | xargs -n1 -I {MOTIF} \
      echo "find mutational_signatures/ -xtype f | xargs -n1 -I{} basename -s .tsv '{}' | sort | xargs -n1 -I{SIGNATURE} echo \"" \
        "echo -n {MOTIF} $'\t' {SIGNATURE} $'\t'; ruby mutate_motif.rb '/home/ilya/iogen/projects/hocomoco11/final_bundle/hocomoco11/full/HUMAN/mono/pcm/{MOTIF}.pcm' 'mutational_signatures/{SIGNATURE}.tsv' genomic_contexts/hg38_contexts.tsv  --report-site-exposure\"" \
        " | bash " \
  | bash > hocomoco11_human_exposure.txt

find good_motifs -xtype f -iname '*.pcm' | xargs -n1 basename -s .pcm | xargs -I{} echo "ruby get_pwm.rb --pcm --pseudocount log good_motifs/{}.pcm > good_motifs_pwm/{}.pwm " | bash
find good_motifs -xtype f -iname '*.ppm' | xargs -n1 basename -s .ppm | xargs -I{} echo "ruby get_pwm.rb --pfm --word-count 100 --pseudocount log good_motifs/{}.ppm > good_motifs_pwm/{}.pwm " | bash

find good_motifs -xtype f -iname '*.ppm' | xargs -n1 basename -s .ppm | xargs -I{} echo "ruby get_pwm.rb --pcm --word-count 100 --to-pcm good_motifs/{}.ppm > good_motifs_pcm/{}.pcm " | bash
find good_motifs -xtype f -iname '*.pcm' | xargs -n1 basename -s .pcm | xargs -I{} echo "ruby get_pwm.rb --pcm --word-count 100 --to-pcm good_motifs/{}.pcm > good_motifs_pcm/{}.pcm " | bash


for TF in $(find good_motifs -xtype f | xargs -n1 basename | cut -d . -f1 | sort -u ); do
  echo $TF;

  find good_motifs_pcm -xtype f -name "$TF.*" \
  | sort \
  | ruby -e '$stdin.readlines.sort_by{|l| File.basename(l).split("@")[1] }.each{|l| puts l }' \
  | xargs align_motifs --pcm \
  | rvm use 2.6.9 do glue_logos --text-size 10 --logo-shift 200 good_motifs_logo_glued/$TF.png --from-stdin
done
