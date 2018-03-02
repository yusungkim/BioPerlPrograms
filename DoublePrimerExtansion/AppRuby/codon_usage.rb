#!/uer/bin/env ruby
#written by yskim@bio.res.titech.ac.jp
# 2006.06.28

require 'bio'

# get sequences
lines = ARGF.read
lines.gsub(/^.*ORIGIN\s*/, '')	#remove abstract
lines.gsub(/[\s\W]/, '')	#remove spaces and non-word
lines.downcase!
seq = Bio::Sequence::NA.new(lines)

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
}

#print result
aa2na.sort.each{ |aa, codons|
	aaThreeLtt = Bio::Sequence::AA.new(aa).codes
	print "#{aaThreeLtt} #{aa}: "
	codons.sort.each{ |codon|
		printf("\t%s => %.3f\n", codon.to_s, codon_usage_by_aa[codon])
	}
}

# print results
#codon_usage.keys.sort.each { |codon|
#	puts "#{codon} : #{codon_usage[codon]}\n"
#}
