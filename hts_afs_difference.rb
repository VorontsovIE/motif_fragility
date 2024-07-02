require 'json'
require 'shellwords'
# {
#   "LEUTX": {
#     "LEUTX.FL@AFS.GFPIVT@lanky-jade-raccoon@Halle.Dimont@Motif_1_w20_astrained.ppm": {
#       "CHS": {
#         "THC_0455.Rep-DIANA_0293": {
#           "Peaks": {
#             "LEUTX.FL@CHS@THC_0455.Rep-DIANA_0293@Peaks.thirsty-chestnut-goat.Val.peaks": {
#               "chipseq_pwmeval_ROC": 1,
#               "chipseq_pwmeval_PR": 1,
#               "chipseq_vigg_ROC": 385,
#               "chipseq_centrimo_neglog_evalue": 160,
#               "combined": 7
#             },
#             "combined": 7
#           },
#           "combined": 7
#         },
# ...
data = JSON.parse(File.read('ranks.freeze-approved.json')); nil

# top_data = data.transform_values{|motifs_data|
#   motifs_data.select{|motif, motif_infos|
#     motif_infos['combined'] == 1 || motif_infos.reject{|k,v| k == 'combined' }.any?{|data_type, data_type_info|
#       data_type_info['combined'] == 1
#     }
#   }
# };nil

# File.write('ranks.winners.freeze-approved.json', top_data.to_json)

top_motifs_from_datatype = data.transform_values{|motifs_data|
  motifs_data.group_by{|motif, motif_infos|
    motif.split('@')[1]
  }.transform_values{|motif_pairs|
    motif_pairs.min_by{|motif, motif_info| motif_info['combined'] }.first
  }
}; nil

File.write('top_motifs_from_datatype.json', top_motifs_from_datatype.to_json)

corresponding_motifs_by_tf = top_motifs_from_datatype.transform_values{|tf_info|
  {
    HTS: tf_info.select{|dt| dt.match?(/HTS/) }.values,
    AFS: tf_info.select{|dt| dt.match?(/AFS/) }.values
  }
}.select{|tf, motifs|
  !motifs[:HTS].empty? && !motifs[:AFS].empty?
} 

corresponding_motifs_by_tf.each{|tf, all_datatypes_motifs|
  all_datatypes_motifs.each{|datatype, motifs|
    motifs.each{|motif|
      cmd = [
        "ruby", "attractiveness.rb",
        "--mutational-signature", "TOPMed_10kb_spectra.tsv:*",
        "--genomic_context", "genomic_contexts/hg38_contexts.tsv",
        "--pfm-motif", "motifs_freeze/#{motif}",
        "--pcm-count", "100",
        "--pseudocount", "log",
        "--background-metrics", "background_metrics.json",  # "background_metrics/bg_metrics_len_{length}.json",
      ].shelljoin + " | "  + [
        "ruby", "-r", "json", "-e",
        "data = {tf: '#{tf}', datatype: '#{datatype}',  motif: '#{motif}'}; metrics = JSON.parse($stdin.read); puts data.merge(result: metrics).to_json"
      ].shelljoin
      puts cmd
        # {tf: tf, datatype: datatype,  motif: motif, result: metrics}
      # metrics = JSON.parse(`#{cmd}`)
      # info = {tf: tf, datatype: datatype,  motif: motif, result: metrics}
      # puts info.to_json
    }
  }
}
