module Bio

#
# bio/util/color_scheme.rb - Popular color codings for nucleic and amino acids
#
# Copyright::  Copyright (C) 2005 Trevor Wennblom <trevor@corevx.com>
# License::    LGPL
#
#  $Id: color_scheme.rb,v 1.2 2005/12/13 14:58:07 trevor Exp $
#
#
#--
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
#
#++
#
#

=begin rdoc
bio/util/color_scheme.rb - Popular color codings for nucleic and amino acids

== Synopsis

The Bio::ColorScheme module contains classes that return popular color codings
for nucleic and amino acids in RGB hex format suitable for HTML code.

The current schemes supported are:
* Buried - Buried index
* Helix - Helix propensity
* Hydropathy - Hydrophobicity
* Nucleotide - Nucelotide color coding
* Strand - Strand propensity
* Taylor - Taylor color coding
* Turn - Turn propensity
* Zappo - Zappo color coding

Planned color schemes include:
* BLOSUM62
* ClustalX
* Percentage Identity (PID)

Color schemes BLOSUM62, ClustalX, and Percentage Identity are all dependent
on the alignment consensus.

This data is currently referenced from the JalView alignment editor.
Clamp, M., Cuff, J., Searle, S. M. and Barton, G. J. (2004), 
"The Jalview Java Alignment Editor," Bioinformatics, 12, 426-7
http://www.jalview.org

Currently the score data for things such as hydropathy, helix, turn, etc. are contained
here but should be moved to bio/data/aa once a good reference is found for these
values.


== Usage

  require 'bio/util/color_scheme'

  seq = 'gattaca'
  scheme = Bio::ColorScheme::Zappo
  postfix = '</span>'
  html = ''
  seq.each_byte do |c|
    color = scheme[c.chr]
    prefix = %Q(<span style="background:\##{color};">)
    html += prefix + c.chr + postfix
  end

  puts html


=== Accessing colors

  puts Bio::ColorScheme::Buried['A']  # 00DC22
  puts Bio::ColorScheme::Buried[:c]   # 00BF3F
  puts Bio::ColorScheme::Buried[nil]  # nil
  puts Bio::ColorScheme::Buried['-']  # FFFFFF
  puts Bio::ColorScheme::Buried[7]    # FFFFFF
  puts Bio::ColorScheme::Buried['junk']  # FFFFFF
  puts Bio::ColorScheme::Buried['t']  # 00CC32


== Author
Trevor Wennblom <trevor@corevx.com>


== Copyright
Copyright (C) 2005 Trevor Wennblom
Licensed under the same terms as BioRuby.

=end
module ColorScheme
  cs_location = 'bio/util/color_scheme'

  # Score sub-classes
  autoload :Buried, "#{cs_location}/buried"
  autoload :Helix, "#{cs_location}/helix"
  autoload :Hydropathy, "#{cs_location}/hydropathy"
  autoload :Strand, "#{cs_location}/strand"
  autoload :Turn, "#{cs_location}/turn"

  # Simple sub-classes
  autoload :Nucleotide, "#{cs_location}/nucleotide"
  autoload :Taylor, "#{cs_location}/taylor"
  autoload :Zappo, "#{cs_location}/zappo"

  # Consensus sub-classes
  # NOTE todo
  # BLOSUM62
  # ClustalX
  # PID

  # A very basic class template for color code referencing.
  class Simple
    def self.[](x)
      return if x.nil?
      # accept symbols and any case
      @colors[x.to_s.upcase]
    end

    def self.colors() @colors end

    #######
    private
    #######

    # Example
    @colors = {
      'A' => '64F73F',
    }
    @colors.default = 'FFFFFF'  # return white by default
  end


  # A class template for color code referencing of color schemes
  # that are score based.  This template is expected to change
  # when the scores are moved into bio/data/aa
  class Score
    def self.[](x)
      return if x.nil?
      # accept symbols and any case
      @colors[x.to_s.upcase]
    end

    def self.min(x) @min end
    def self.max(x) @max end
    def self.scores() @scores end
    def self.colors() @colors end

    #########
    protected
    #########

    def self.percent_to_hex(percent)
      percent = percent.to_f if percent.is_a?(String)
      if (percent > 1.0) or (percent < 0.0) or percent.nil?
        raise 'Percentage must be between 0.0 and 1.0'
      end
      "%02X" % (percent * 255.0)
    end

    def self.rgb_percent_to_hex(red, green, blue)
      percent_to_hex(red) + percent_to_hex(green) + percent_to_hex(blue)
    end

    def self.score_to_percent(score, min, max)
      # .to_f to ensure every operation is float-aware
      percent = (score.to_f - min) / (max.to_f - min)
      percent = 1.0 if percent > 1.0
      percent = 0.0 if percent < 0.0
      percent
    end

    #######
    private
    #######

    # Example
    def self.score_to_rgb_hex(score, min, max)
      percent = score_to_percent(score, min, max)
      rgb_percent_to_hex(percent, 0.0, 1.0-percent)
    end

    @colors = {}
    @scores = {
      'A' => 0.83,
    }
    @min = 0.37
    @max = 1.7
    @scores.each { |k,s| @colors[k] = score_to_rgb_hex(s, @min, @max) }
    @colors.default = 'FFFFFF'  # return white by default

  end


  # NOTE todo
  class Consensus
  end

end  # module ColorScheme
end  # module Bio
