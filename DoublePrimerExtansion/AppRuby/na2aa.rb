#!/usr/bin/env ruby
# Usage: ruby na2aa.rb DNA_file
# Written by yskim June 14th, 2006

require 'bio'

# environments
# amino acid extract sensitivity
sensitivity = 28 # length of amino acid

input_seq = ARGF.read
my_naseq = Bio::Sequence::NA.new(input_seq)

# extract NA field
if( my_naseq =~ /ORIGIN\s*/i )
  my_naseq = $'
end
# remove numbers
my_naseq.delete! "0-9"
# remove spaces
my_naseq.delete! "\s"

# translate all possible reading frame
# Caution!! This version do not check the 'complement seq'
# my_aaseq = my_aaseq.complement
1.upto(3) do |frame_shift|
  starting = frame_shift
  my_aaseq = my_naseq.translate(frame_shift)
  orf_count = 0
  puts "#{"*"*20} Reading Frame Shift = #{frame_shift} #{"*"*20}"
  puts my_aaseq
  while (my_aaseq =~ /[mM]\w+\*/) do
    ending = starting + $&.length*3 - 1
    bp = ending - starting + 1
    if ($&.length > sensitivity)
      orf_count = orf_count + 1
      print "\t\tORF No. #{orf_count}\t#{$&.length-1}aa\t"
      print "#{bp}bp(#{starting}..#{ending})\n"
      puts $&
    end
    my_aaseq = $'
    starting = ending
  end
  puts
end
