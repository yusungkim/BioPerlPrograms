#!/usr/bin/env ruby
# Usage: ruby na2aa.rb DNA_file
# Written by yskim June 14th, 2006

require 'bio'

# environments
# amino acid extract sensitivity
sensitivity = 150 # length of amino acid

# get sequences
lines = ARGF.read
lines.gsub(/^.*ORIGIN\s*/i, '') #remove abstract
lines.gsub(/[\s\W]/, '')  #remove spaces and non-word
lines.downcase!
my_naseq = Bio::Sequence::NA.new(lines)

# translate all possible reading frame
# Caution!! This version do not check the 'complement seq'
# my_aaseq = my_aaseq.complement
# make ORF-concat seq
seq = Bio::Sequence::NA.new('')
1.upto(3) do |frame_shift|
  starting = frame_shift
  my_aaseq = my_naseq.translate(frame_shift, 1, 11)#11=>Bacteria
  orf_count = 0
  #puts "#{"*"*20} Reading Frame Shift = #{frame_shift} #{"*"*20}"
  while (my_aaseq =~ /[mM]\w+\*/) do
    ending = starting + $&.length*3 - 1
    bp = ending - starting + 1
    if ($&.length > sensitivity)
      orf_count = orf_count + 1
      #print "\t\tORF No. #{orf_count}\t#{$&.length-1}aa\t"
      #print "#{bp}bp(#{starting}..#{ending})\n"
      #puts $&
      seq += my_naseq[starting..ending]
    end
    my_aaseq = $'
    starting = ending
  end
  puts
end

#calculate codon usages
codon_usage = Hash.new(0)
seq.window_search(3,3) { |codon|
  codon_usage[codon] += 1
}

SP = 11 #using Bacteria type
codon_table = Bio::CodonTable[SP] #using 11=>Bacteria type

#make aa->na hash table
aa2na = Hash.new
codon_table.each { |codon, aa|
  if (aa2na[aa] == nil)
    aa2na[aa] = Array[codon]
  else
    aa2na[aa].push codon
  end
}

#make codon_usage at each aa by rate
codon_usage_by_aa = Hash.new(0.0)
aa2na.each{ |aa, codons|
  sum = 0
  codons.each{ |codon|
    sum += codon_usage[codon]
  }
  codons.each{ |codon|
    codon_usage_by_aa[codon] = codon_usage[codon].to_f / sum if sum != 0
  }
  puts sum
}

#print result
aa2na.sort.each{ |aa, codons|
  aaThreeLtt = Bio::Sequence::AA.new(aa).codes
  print "#{aaThreeLtt} #{aa}: "
  codons.sort.each{ |codon|
    printf("\t%s => %.3f\n", codon.to_s, codon_usage_by_aa[codon])
  }
}
