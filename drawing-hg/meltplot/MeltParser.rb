#!/usr/bin/env ruby

# Parser of melt data, first hack version, written by JSB and MBP.

prefix = ARGV[0]
tempStart = (ARGV[1]).to_f
tempIncrement = (ARGV[2]).to_f
tempEnd = (ARGV[3]).to_f
outputFile = ARGV[4]
temperatures = Array.new
temp = tempStart
nTemps = 0
while temp <= tempEnd
  temperatures << temp
  temp += tempIncrement
  nTemps += 1
end #while

fUnpaired = Array.new # Fraction of srands unpaired

temperatures.each do |x|
  
  number_of_strands = 0
  n = 0 # Total strand length
  fUnpaired << 0.0
  
  inputFileName = sprintf("./%.1f/%s.fpairs",x.to_f,prefix)

  File.open(inputFileName) do |f|
    seen_id_seq = false
    f.each_line do |l|
      line = l.chomp

      # Fin
      if line =~ /% Number of strands: (\d+)/
        number_of_strands = $1.to_i
        next
      end #if

      # Check to see if we're at the line with % id sequence
      if line =~ /% id sequence/
        seen_id_seq = 1
        next
      end #if

      # Now we're parsing the sequences
      if seen_id_seq
        line =~ /%  #{seen_id_seq} (.+)$/
        sequence = $1
        n += sequence.length
        seen_id_seq += 1
        if seen_id_seq > number_of_strands
          seen_id_seq = false # go back to normal parsing
        end
        next
      end #if

      # Blast past line if it's a comment or newline
      next if line =~ /^%/

      # Split the line if it's got data
      chunks = l.chomp.split(/\s+/)

      if chunks[1].to_i == n+1  # If it's an entry for unpaired
        fUnpaired[fUnpaired.length - 1] += (chunks[2]).to_f
      end #if

    end #do
  end #do

  # Have to divide by n
  fUnpaired[fUnpaired.length - 1] /= n.to_f
  
end #do

# Write results to a file
outfile = File.open(outputFile,"w")
for i in 0..(nTemps-1)
  outfile.printf("%.1f\t%8.7f\n",temperatures[i],fUnpaired[i])
end #for
outfile.close


