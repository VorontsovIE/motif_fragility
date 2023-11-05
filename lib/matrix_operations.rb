require 'tempfile'

def read_matrix(fn, num_columns: 4)
  if fn == 'stdin'
    lines = $stdin.readlines
    name ||= 'motif'
  else
    lines = (fn == 'stdin') ? $stdin.readlines : File.readlines(fn)
    name = File.basename(fn, File.extname(fn))
  end
  lines = lines.map(&:strip)
  rows = lines.map{|l| l.split }
  unless (rows[0].size == num_columns) && rows[0].all?{|x| Float(x, exception: false) }
    hdr = lines.first
    rows.shift
    name = (hdr.start_with?('>') ? hdr[1..-1].strip : hdr.strip).split.first
  end
  matrix = rows.map{|row|
    row.map{|x| Float(x) }
  }
  raise  if matrix.empty?
  raise  unless matrix.all?{|row| row.size == num_columns }
  {name: name, matrix: matrix}
end

#############################

def matrix_as_string(model, transpose_output: false)
  res = [">#{model[:name]}"]
  matrix = transpose_output ? model[:matrix].transpose : model[:matrix]
  res += matrix.map{|row| row.join("\t") }
  res.join("\n")
end

#############################

def calculate_pseudocount(count, pseudocount: :log)
  case pseudocount
  when :log
    Math.log([count, 2].max);
  when :sqrt
    Math.sqrt(count)
  else Numeric
    pseudocount
  end
end

#############################

def pcm2pfm(pcm)
  pfm_matrix = pcm[:matrix].map{|row|
    norm = row.sum
    row.map{|x| x.to_f / norm }
  }
  {name: pcm[:name], matrix: pfm_matrix}
end

def pfm2pcm(pfm, word_count: 100)
  pcm_matrix = pfm[:matrix].map{|row|
    row.map{|el| el * word_count }
  }
  {name: pfm[:name], matrix: pcm_matrix}
end

def pcm2pwm(pcm, pseudocount: :log)
  pwm_matrix = pcm[:matrix].map{|row|
    count = row.sum
    row.map{|el|
      pseudocount_value = calculate_pseudocount(count, pseudocount: pseudocount)
      numerator = el + 0.25 * pseudocount_value
      denominator = 0.25 * (count + pseudocount_value)
      Math.log(numerator / denominator)
    }
  }
  {name: pcm[:name], matrix: pwm_matrix}
end

#############################

def make_thresholds(pwm_fn, thresholds_fn: nil)
  if !thresholds_fn
    thresholds_file = Tempfile.new.tap(&:close)
    thresholds_fn = thresholds_file.path
    $keep_tempfiles ||= []
    $keep_tempfiles << thresholds_file 
  end
  cmd = "java -cp ape.jar ru.autosome.ape.PrecalculateThresholds #{pwm_fn} --single-motif > #{thresholds_fn}"
  system(cmd)
  thresholds_fn
end

# def make_pwm_from_pfm(pfm_fn, pwm_fn: nil)
#   pfm = read_matrix(pfm_fn)
#   pcm = pfm2pcm(pfm, word_count: 100)
#   pwm = pcm2pwm(pcm, pseudocount: :log)

#   if !pwm_fn
#     pwm_file = Tempfile.new{|fw|
#       fw.puts matrix_as_string(pwm)
#     }
#     pwm_fn = pwm_file.path
#     $keep_tempfiles ||= []
#     $keep_tempfiles << pwm_file
#   else
#     File.write(pwm_fn, matrix_as_string(pwm))
#   end

#   pwm_fn
# end

def write_temp_matrix(matrix)
  pwm_file = Tempfile.new.tap(&:close)
  File.write(pwm_file.path, matrix_as_string({name: '>matrix', matrix: matrix.matrix}))
  $keep_tempfiles ||= []
  $keep_tempfiles << pwm_file

  pwm_file.path
end
