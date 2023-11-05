require 'set'
require 'fileutils'

FileUtils.mkdir_p 'mutations'

alexandrov_folder = '/home/ilya/iogen/projects/somatic_mutations/source_data/AlexandrovEtAl'
wg_samples_by_cl = File.readlines("#{alexandrov_folder}/samples_summary.txt").map{|l|
  l.chomp.split("\t")
}.select{|cl, sample, seq_type, source|
  seq_type == 'Whole genome'
}.group_by{|cl, sample, seq_type, source|
  cl
}.transform_values{|rows|
  rows.map{|cl, sample, seq_type, source| sample }.to_set
}

wg_samples_by_cl.each{|cl, wg_samples|
  cl_nospace = cl.gsub(' ', '_')
  output_filename = "raw_mutations/#{cl_nospace}.tsv"
  File.open(output_filename, 'w') do |fw|
    cl_muts_filename = "#{alexandrov_folder}/somatic_mutation_data/#{cl}/#{cl}_clean_somatic_mutations_for_signature_analysis.txt"
    File.readlines(cl_muts_filename).map{|l|
      sample, mut_type, chr, pos_from, pos_to, from, to, info = l.chomp.split("\t")
      [sample, mut_type, chr, pos_from.to_i, pos_to.to_i, from, to]
    }.select{|sample, mut_type, chr, pos_from, pos_to, from, to|
      wg_samples.include?(sample) && chr != 'MT' && mut_type == 'subs'
    }.each{|sample, mut_type, chr, pos_from, pos_to, from, to|
      infos = ["chr#{chr}", pos_from - 1, pos_from, "#{cl_nospace};#{sample};chr#{chr}:#{pos_from}@[#{from}/#{to}]"]
      fw.puts infos.join("\t")
    }
  end
}
