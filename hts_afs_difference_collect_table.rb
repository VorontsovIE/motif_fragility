require 'json'

def get_signature_motif_info(motif_info, signature)
  sig_info = motif_info['result'][signature]
  {
    'attractiveness-value' => sig_info['attractiveness']['value'],
    'attractiveness-zscore' => sig_info['attractiveness']['zscore'],
    'pvalue_change-value' => sig_info['pvalue_change']['value'],
    'pvalue_change-zscore' => sig_info['pvalue_change']['zscore'],
  }
end


SIGNATURES = ["comp.1", "comp.2", "comp.3", "comp.4", "comp.5", "comp.6", "comp.7", "comp.8", "comp.9", "comp.10", "comp.11", "comp.12", "comp.13", "comp.14", "offset",]

METRICS = [
  'attractiveness-value', 'attractiveness-zscore',
  'pvalue_change-value', 'pvalue_change-zscore',
]

all_infos = File.readlines('top_hts_afs_motifs_metrics.jsonl').map{|l| JSON.parse(l) }

header = [
  'tf', 'motif_hts', 'motif_afs',
  *SIGNATURES.flat_map{|sig|
    basic_metrics = METRICS.flat_map{|metric|
      ['HTS', 'AFS'].map{|datatype|
        "#{sig}:#{metric}:#{datatype}"
      }
    }
    [
      # *basic_metrics,
      "#{sig}:hts_to_afs_attractiveness",
    ]
  }
]

puts header.join("\t")

all_infos.group_by{|motif_info|
  motif_info.values_at('tf')
}.transform_values{|tf_infos|
  hts_motifs = tf_infos.select{|motif_info| motif_info['datatype'] == 'HTS' }
  afs_motifs = tf_infos.select{|motif_info| motif_info['datatype'] == 'AFS' }
  hts_motifs.product(afs_motifs)
}.flat_map{|tf, hts_afs_pairs|
  # hts_afs_pairs.first(1).map{|hts_motif_info, afs_motif_info|
  hts_afs_pairs.map{|hts_motif_info, afs_motif_info|
    all_sig_row = SIGNATURES.map{|signature|
      hts_sig_info = get_signature_motif_info(hts_motif_info, signature)
      afs_sig_info = get_signature_motif_info(afs_motif_info, signature)
      basic_metrics = METRICS.flat_map{|metric|
        [hts_sig_info[metric], afs_sig_info[metric]]
      }
      [
        # *basic_metrics, 
        hts_sig_info['attractiveness-value'] / afs_sig_info['attractiveness-value']
      ]
    }
    row = [tf, hts_motif_info['motif'], afs_motif_info['motif'], *all_sig_row,]
  }
  # all_infos.map{|motif_info|
  #   all_sig_row = SIGNATURES.flat_map{|signature|
  #     sig_row = get_signature_motif_info(motif_info, signature].values_at(*[
  #       'attractiveness-value', 'attractiveness-zscore',
  #       'pvalue_change-value', 'pvalue_change-zscore',
  #     ])
  #     sig_row
  #   }
  #   row = [motif_info['tf'], motif_info['datatype'], motif_info['motif'], *all_sig_row,]
  #   puts row.join("\t")
  # }
}.each{|row|
  puts row.join("\t")
}
