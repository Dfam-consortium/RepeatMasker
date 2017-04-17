#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) dupliconToSVG.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to convert a DupMasker *.duplicons file
##      into a SVG visualization of the data.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
# ChangeLog
#
#     $Log: dupliconToSVG.pl,v $
#     Revision 1.29  2017/02/01 21:01:57  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

dupliconToSVG.pl - Create a SVG file from a DupMasker *.dupicons file.

=head1 SYNOPSIS

  dupliconToSVG.pl [-version] <*.dupliconfile>

=head1 DESCRIPTION

 A utility script to convert a DupMasker *.duplicons file 
 into a SVG visualization of the data.  A separate SVG file 
 is created for each sequence in the duplicons file.  They
 are named as: file.duplicons.1.svg, file.duplicons.2.svg and so on.

 SVG files may be viewed in most current web browsers.  Firefox
 appears to have adequate SVG support for these types of files.
 Once the SVG Print format standard is finalized this utility 
 will be modified to create a file which contains page breaks
 at the appropriate positions.

 The current version only works with the human duplicon lib.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2008 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;

my $Version = "1.0";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

#
# ARGV Processing
#
if ( !@ARGV ) {
  usage();
}

my $dupliconsFile = $ARGV[ 0 ];

# Layout constants
my $docWidth              = 8.5;
my $docHeightIncr         = 11;
my $rulerTicUnitIncr      = 80000;
my $rulerTicLabelDiv      = 1000;
my $rulerTicsPerLine      = 12;
my $rulerTicLabelFontSize = 10;
my $sideMargin            = .5;
my $topBottomMargin       = .5;
my $rulerTicHeight        = .1;
my $dupBlockHeight        = .1;

# TODO: This could be derived
my $rulerLineVerticalIncr = 1;

my %colors = getColors();

##
##  Read duplicons file
##
open IN, "<$dupliconsFile" || die "Could not open $dupliconsFile!\n";
my %sequences = ();
while ( <IN> ) {
  if ( /^\d+/ ) {
    my @fields = split;
    $fields[ 4 ] =~ s/\|/_/g;
    if ( !defined $sequences{ $fields[ 4 ] } ) {
      my $remainingLength = $fields[ 7 ];
      $remainingLength =~ s/[\(\)]//g;
      $sequences{ $fields[ 4 ] }->{'length'} =
          ( $fields[ 6 ] + $remainingLength );
    }
    my $chromosome;
    my $cytoband;
    my $ID = $fields[ 8 ];
    $ID = $fields[ 9 ] if ( $ID eq "C" );
    if ( $ID =~ /\S+\|(\S+)\:\d+-\d+\|(\S+)/ ) {
      $chromosome = $1;
      $cytoband   = $2;
      if ( $cytoband ) {
        if (    $cytoband ne "NA"
             && $cytoband ne "NULL" )
        {
          $cytoband = "$chromosome" . "_$cytoband";
        }
      }
      else {
        $cytoband = "ND";
      }
    }
    else { die "Error: ID Field not in expected format: $ID\n"; }
    push @{ $sequences{ $fields[ 4 ] }->{'duplicons'} },
        {
          'start' => $fields[ 5 ],
          'end'   => $fields[ 6 ],
          'chr'   => $chromosome,
          'cyt'   => $cytoband
        };
  }
}
close IN;

#
# Loop through sequences in duplicons file and
# create a separate .svg file for each.
#
my $seqIdx = 1;
foreach my $sequence ( keys( %sequences ) ) {
  my $outFile = $dupliconsFile . "." . $seqIdx++ . ".svg";
  my $seqSize = $sequences{$sequence}->{'length'};

  # Derived factors
  my $canvasWidth  = $docWidth -      ( 2 * $sideMargin );
  my $canvasHeight = $docHeightIncr - ( 2 * $topBottomMargin );
  my $inchesPerBase = $canvasWidth / ( $rulerTicUnitIncr * $rulerTicsPerLine );
  my $basesPerLine  = $rulerTicUnitIncr * $rulerTicsPerLine;
  my $rulerLines    = $seqSize / $basesPerLine;
  my $docHeight     =
      ( $rulerLineVerticalIncr * ( $rulerLines + 1 ) ) + $topBottomMargin;

  ##
  ## Print SVG Header
  ##
  open OUT, ">$outFile" || die "Could not open $outFile for writing!\n";
  print "Creating file $outFile for sequence $sequence...\n";
  print OUT "<?xml version=\"1.0\" standalone=\"no\"?>\n";
  print OUT
"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n";
  print OUT "<svg width=\"8.5in\" height=\"$docHeight"
      . "in\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";

  ##
  ##  Draw ruler
  ##
  for ( my $i = 0 ; $i <= $rulerLines ; $i++ ) {
    ## Print ruler line
    my $x1 = $sideMargin;
    my $y1 = $topBottomMargin + ( $rulerLineVerticalIncr * $i );
    my $x2 = $sideMargin + $canvasWidth;
    if ( ( $basesPerLine * ( $i + 1 ) ) > $seqSize ) {
      $x2 =
          $sideMargin + ( $seqSize - ( $basesPerLine * $i ) ) * $inchesPerBase;
    }
    my $y2 = $y1;
    print OUT "<line x1=\"$x1"
        . "in\" y1=\"$y1"
        . "in\" x2=\"$x2"
        . "in\" y2=\"$y2"
        . "in\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>\n";

    ## Print tics and labels
    for ( my $j = 0 ; $j < $rulerTicsPerLine ; $j++ ) {
      last
          if ( ( $basesPerLine * $i ) + ( $rulerTicUnitIncr * $j ) > $seqSize );
      $x1 = ( ( $rulerTicUnitIncr * $j ) * $inchesPerBase ) + $sideMargin;
      $y1 = $topBottomMargin + ( $rulerLineVerticalIncr * $i );
      $x2 = $x1;
      $y2 = $y1 + $rulerTicHeight;
      print OUT "<line x1=\"$x1"
          . "in\" y1=\"$y1"
          . "in\" x2=\"$x2"
          . "in\" y2=\"$y2"
          . "in\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>\n";
      my $x = $x1;
      my $y = $y2 + .1;
      print OUT "<text x=\"$x"
          . "in\" y=\"$y"
          . "in\" font-size=\"10\" font-family=\"Verdana\">\n";
      print OUT ""
          . ( ( ( $basesPerLine * $i ) + ( $rulerTicUnitIncr * $j ) ) /
              $rulerTicLabelDiv )
          . "\n";
      print OUT "</text>\n";
    }
  }

  ##
  ## Draw data track
  ##
  foreach my $dup ( @{ $sequences{$sequence}->{'duplicons'} } ) {
    my $color = $colors{ $dup->{'cyt'} };
    my $r;
    my $g;
    my $b;
    if ( $color =~ /([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})/ ) {
      $r = hex( $1 );
      $g = hex( $2 );
      $b = hex( $3 );
    }
    else { die "Color $dup->{'cyt'} is not in recognized format >$color<\n"; }

    my $start     = $dup->{'start'};
    my $end       = $dup->{'end'};
    my $lineStart = int( $start / $basesPerLine );
    my $lineEnd   = int( $end / ( $basesPerLine ) );
    for ( my $i = $lineStart ; $i <= $lineEnd ; $i++ ) {
      my $x =
          ( ( $start - ( $basesPerLine * $i ) ) * $inchesPerBase ) +
          $sideMargin;
      $x = $sideMargin if ( $x < $sideMargin );
      my $y     = ( $i * $rulerLineVerticalIncr ) + $topBottomMargin + .25;
      my $width = ( $end - $start ) * $inchesPerBase;
      $width = ( $sideMargin + $canvasWidth ) - $x
          if ( ( $x + $width ) > ( $sideMargin + $canvasWidth ) );
      my $height = .1;
      print OUT "<rect x=\"$x"
          . "in\" y=\"$y"
          . "in\" width=\"$width"
          . "in\" height=\"$height"
          . "in\" style=\"fill:rgb($r,$g,$b);stroke-width=0\"/>\n";
    }
  }
  print OUT "</svg>\n";
  close OUT;

}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my %colors = getColors();
##
##  Pre-defined colors for human duplicon lib.  NOTE: As new duplicon
##  libraries become available it might make more sense to develop
##  an algorithm for sampling colors from the CIELab space ( a color
##  space where the frequency spectrum is biased towards colors
##  which are easily perceived by the human eye.  Euclidean distance
##  in this space translates to perceived difference between colors. ).
##
##  Returns
##      A hash mapping the human chromosomes and cytogenetic bands
##   to distinctive RGB colors.
##
##-------------------------------------------------------------------------##
sub getColors {
  my %colors = (

    #
    # Chromosome Colors
    #
    "chr1"       => "FFE4E1",
    "chr2"       => "000000",
    "chr3"       => "2F4F4F",
    "chr4"       => "191970",
    "chr5"       => "6495ED",
    "chr6"       => "0000CD",
    "chr7"       => "00BFFF",
    "chr8"       => "4682B4",
    "chr9"       => "00FFFF",
    "chr10"      => "E0FFFF",
    "chr11"      => "006400",
    "chr12"      => "00FF7F",
    "chr13"      => "7CFC00",
    "chr14"      => "FFFF00",
    "chr15"      => "FFD700",
    "chr16"      => "BC8F8F",
    "chr17"      => "CD5C5C",
    "chr18"      => "A52A2A",
    "chr19"      => "FF6347",
    "chr20"      => "CD853F",
    "chr21"      => "FF0000",
    "chr22"      => "FF69B4",
    "chrX"       => "FF1493",
    "chrY"       => "800000",
    "chr_random" => "9400D3",
    "No hit"     => "BEBEBE",
    "ND"         => "BEBEBE",
    "NA"         => "BEBEBE",
    "NULL"       => "BEBEBE",
    "Binary"     => "BEBEBE",
    "small"      => "BEBEBE",
    "No Net"     => "BEBEBE",

    #
    # Cytoband Colors
    #
    "chr1_p36.33"  => "B88444",
    "chr1_p21.1"   => "11287B",
    "chr1_p13.3"   => "4ABB5C",
    "chr1_p13.2"   => "581107",
    "chr1_p13.1"   => "1833C7",
    "chr1_p12"     => "1136C0",
    "chr1_p11.2"   => "598788",
    "chr1_p11.1"   => "EDCAB9",
    "chr1_p36.21"  => "581BD6",
    "chr1_q11"     => "DC58B4",
    "chr1_q12"     => "9A4408",
    "chr1_q21.1"   => "825640",
    "chr1_q21.2"   => "5125DE",
    "chr1_q21.3"   => "5E0474",
    "chr1_q22"     => "C2BD63",
    "chr1_q23.1"   => "B495CD",
    "chr1_q23.2"   => "6B20B2",
    "chr1_p36.13"  => "509915",
    "chr1_q23.3"   => "A995DA",
    "chr1_q24.1"   => "9C0D9B",
    "chr1_q24.2"   => "B236D5",
    "chr1_q24.3"   => "2EE61B",
    "chr1_q25.1"   => "C96585",
    "chr1_q25.2"   => "43BE28",
    "chr1_q25.3"   => "961118",
    "chr1_q31.1"   => "AAEE50",
    "chr1_q31.2"   => "798B1C",
    "chr1_q31.3"   => "D3D9A1",
    "chr1_q32.1"   => "817EA6",
    "chr1_p36.12"  => "408BE6",
    "chr1_q32.2"   => "A87401",
    "chr1_q32.3"   => "A63174",
    "chr1_q41"     => "7D207C",
    "chr1_q42.11"  => "A0BE47",
    "chr1_q42.12"  => "168893",
    "chr1_q42.13"  => "890CD7",
    "chr1_q42.2"   => "90B66D",
    "chr1_p36.32"  => "7047E0",
    "chr1_q42.3"   => "1634E1",
    "chr1_q43"     => "794376",
    "chr1_p36.11"  => "8E9CB0",
    "chr1_q44"     => "EE12B1",
    "chr1_p35.3"   => "6A9B94",
    "chr1_p35.2"   => "8B859B",
    "chr1_p35.1"   => "B85AE1",
    "chr1_p34.3"   => "CB14A2",
    "chr1_p34.2"   => "54886D",
    "chr1_p34.1"   => "620526",
    "chr1_p33"     => "78ED83",
    "chr1_p32.3"   => "D3409E",
    "chr1_p36.31"  => "867A32",
    "chr1_p32.2"   => "2B9D08",
    "chr1_p32.1"   => "D8DE70",
    "chr1_p31.3"   => "8C2CA7",
    "chr1_p31.2"   => "A9E875",
    "chr1_p31.1"   => "5C7AB1",
    "chr1_p36.23"  => "212EE2",
    "chr1_p22.3"   => "566B4C",
    "chr1_p22.2"   => "A5C3D7",
    "chr1_p36.22"  => "8B5A3E",
    "chr1_p22.1"   => "BA45E9",
    "chr1_p21.3"   => "AB4667",
    "chr1_p21.2"   => "84EE86",
    "chr10_p15.3"  => "850363",
    "chr10_q24.31" => "3D8B70",
    "chr10_q24.32" => "73A3B3",
    "chr10_q24.33" => "2D8323",
    "chr10_q25.1"  => "32974B",
    "chr10_q25.2"  => "D45DC2",
    "chr10_q25.3"  => "572E5C",
    "chr10_q26.11" => "736B82",
    "chr10_q26.12" => "46EE4E",
    "chr10_p13"    => "D5B278",
    "chr10_q26.13" => "C47554",
    "chr10_q26.2"  => "357AD6",
    "chr10_q26.3"  => "2AD0A7",
    "chr10_p12.33" => "C64D51",
    "chr10_p12.32" => "7D6BE5",
    "chr10_p12.31" => "601E7E",
    "chr10_p12.2"  => "32AD7A",
    "chr10_p12.1"  => "E8126A",
    "chr10_p11.23" => "506880",
    "chr10_p15.2"  => "510259",
    "chr10_p11.22" => "3B3BDE",
    "chr10_p11.21" => "3B283E",
    "chr10_p15.1"  => "E45220",
    "chr10_p11.1"  => "C72D31",
    "chr10_q11.1"  => "6C641A",
    "chr10_q11.21" => "74302B",
    "chr10_q11.22" => "BAE069",
    "chr10_q11.23" => "22B4A8",
    "chr10_q21.1"  => "65B366",
    "chr10_q21.2"  => "E1B067",
    "chr10_q21.3"  => "562610",
    "chr10_p14"    => "0E3375",
    "chr10_q22.1"  => "745BAA",
    "chr10_q22.2"  => "B0E6CE",
    "chr10_q22.3"  => "9E8749",
    "chr10_q23.1"  => "9A188E",
    "chr10_q23.2"  => "ADC8C2",
    "chr10_q23.31" => "53A12A",
    "chr10_q23.32" => "3656B2",
    "chr10_q23.33" => "BE801D",
    "chr10_q24.1"  => "E74330",
    "chr10_q24.2"  => "329A11",
    "chr11_p15.5"  => "82A736",
    "chr11_q22.2"  => "884606",
    "chr11_q22.3"  => "1DEC79",
    "chr11_p15.3"  => "029727",
    "chr11_q23.1"  => "B4380A",
    "chr11_q23.2"  => "B31503",
    "chr11_q23.3"  => "E28060",
    "chr11_q24.1"  => "35966A",
    "chr11_q24.2"  => "7EEB3E",
    "chr11_p15.2"  => "9B3EA3",
    "chr11_q24.3"  => "19AC0D",
    "chr11_q25"    => "14E1B1",
    "chr11_p15.1"  => "1ACB50",
    "chr11_p14.3"  => "3B5147",
    "chr11_p14.2"  => "E2252E",
    "chr11_p14.1"  => "BD20AD",
    "chr11_p15.4"  => "36EC8B",
    "chr11_p13"    => "B09345",
    "chr11_p12"    => "AD348B",
    "chr11_p11.2"  => "867BD2",
    "chr11_p11.12" => "662BDB",
    "chr11_p11.11" => "05C701",
    "chr11_q11"    => "4DE531",
    "chr11_q12.1"  => "5A05DA",
    "chr11_q12.2"  => "8A6BAA",
    "chr11_q12.3"  => "A7662A",
    "chr11_q13.1"  => "22BB30",
    "chr11_q13.2"  => "D6A78E",
    "chr11_q13.3"  => "0664A2",
    "chr11_q13.4"  => "6495A9",
    "chr11_q13.5"  => "B217E2",
    "chr11_q14.1"  => "EA14A5",
    "chr11_q14.2"  => "83D0D7",
    "chr11_q14.3"  => "450DA1",
    "chr11_q21"    => "1C5BCB",
    "chr11_q22.1"  => "3C8849",
    "chr12_p13.33" => "E08E6A",
    "chr12_q23.2"  => "197E46",
    "chr12_p13.2"  => "B9D51D",
    "chr12_q23.3"  => "EB60CA",
    "chr12_q24.11" => "A78719",
    "chr12_q24.12" => "3AE519",
    "chr12_q24.13" => "34A5B1",
    "chr12_q24.21" => "CB5B4A",
    "chr12_q24.22" => "E9402C",
    "chr12_q24.23" => "561A72",
    "chr12_q24.31" => "705397",
    "chr12_q24.32" => "2D4C2C",
    "chr12_q24.33" => "58E06E",
    "chr12_p13.1"  => "39086A",
    "chr12_p12.3"  => "B4CA9A",
    "chr12_p12.2"  => "574A5E",
    "chr12_p12.1"  => "251D48",
    "chr12_p11.23" => "3CE71C",
    "chr12_p11.22" => "81EE1B",
    "chr12_p11.21" => "E6C78A",
    "chr12_p13.32" => "9136DE",
    "chr12_p11.1"  => "A91BB4",
    "chr12_q11"    => "182E20",
    "chr12_q12"    => "123C1E",
    "chr12_q13.11" => "641383",
    "chr12_q13.12" => "EA46B9",
    "chr12_q13.13" => "7D05E2",
    "chr12_p13.31" => "CDD82C",
    "chr12_q13.2"  => "8E5918",
    "chr12_q13.3"  => "84630C",
    "chr12_q14.1"  => "02460E",
    "chr12_q14.2"  => "B7E9B5",
    "chr12_q14.3"  => "DA9B4C",
    "chr12_q15"    => "E1E02A",
    "chr12_q21.1"  => "7A8500",
    "chr12_q21.2"  => "91C983",
    "chr12_q21.31" => "01ADCA",
    "chr12_q21.32" => "93289D",
    "chr12_q21.33" => "7541BD",
    "chr12_q22"    => "110636",
    "chr12_q23.1"  => "63318A",
    "chr13_p13"    => "2C40B2",
    "chr13_q33.1"  => "03A25A",
    "chr13_q33.2"  => "D0D6C2",
    "chr13_q33.3"  => "0ED23D",
    "chr13_q34"    => "7A4647",
    "chr13_p11.1"  => "DD3D2A",
    "chr13_q11"    => "A690C0",
    "chr13_q12.11" => "98EA7C",
    "chr13_q12.12" => "E30C37",
    "chr13_q12.13" => "30707B",
    "chr13_q12.2"  => "107173",
    "chr13_q12.3"  => "B89432",
    "chr13_q13.1"  => "1ED6E6",
    "chr13_q13.2"  => "ECBE45",
    "chr13_q13.3"  => "B67A51",
    "chr13_p12"    => "8E0468",
    "chr13_q14.11" => "4C8D22",
    "chr13_q14.12" => "1D9ABD",
    "chr13_q14.13" => "712532",
    "chr13_q14.2"  => "430848",
    "chr13_q14.3"  => "D2EAB6",
    "chr13_q21.1"  => "037CE8",
    "chr13_q21.2"  => "C11158",
    "chr13_q21.31" => "D3CBB3",
    "chr13_q21.32" => "DB0359",
    "chr13_q21.33" => "39C0AA",
    "chr13_q22.1"  => "66ACAE",
    "chr13_q22.2"  => "4C7A4A",
    "chr13_q22.3"  => "80245C",
    "chr13_q31.1"  => "AC40E2",
    "chr13_p11.2"  => "71486C",
    "chr13_q31.2"  => "71E041",
    "chr13_q31.3"  => "E3BB8D",
    "chr13_q32.1"  => "65A418",
    "chr13_q32.2"  => "53E211",
    "chr13_q32.3"  => "EBABD8",
    "chr14_p13"    => "A529B4",
    "chr14_q32.31" => "00154C",
    "chr14_q32.32" => "553B93",
    "chr14_q32.33" => "33D8EA",
    "chr14_p11.1"  => "CE1B91",
    "chr14_q11.1"  => "9DE94B",
    "chr14_q11.2"  => "5A9D81",
    "chr14_q12"    => "5084A9",
    "chr14_p12"    => "55A943",
    "chr14_q13.1"  => "0C41C5",
    "chr14_q13.2"  => "CEC8C8",
    "chr14_q13.3"  => "8489EA",
    "chr14_q21.1"  => "D568E5",
    "chr14_q21.2"  => "8A085B",
    "chr14_q21.3"  => "005552",
    "chr14_q22.1"  => "53886E",
    "chr14_q22.2"  => "6D4393",
    "chr14_q22.3"  => "5C524A",
    "chr14_q23.1"  => "A0C210",
    "chr14_q23.2"  => "B7E58B",
    "chr14_q23.3"  => "CE9CAD",
    "chr14_p11.2"  => "118DD2",
    "chr14_q24.1"  => "6158AD",
    "chr14_q24.2"  => "51E78E",
    "chr14_q24.3"  => "332EAA",
    "chr14_q31.1"  => "887D22",
    "chr14_q31.2"  => "8C495A",
    "chr14_q31.3"  => "B45BD7",
    "chr14_q32.11" => "A8C5AE",
    "chr14_q32.12" => "658BCD",
    "chr14_q32.13" => "928854",
    "chr14_q32.2"  => "CDA2AC",
    "chr15_p13"    => "A23814",
    "chr15_p11.1"  => "147825",
    "chr15_q11.1"  => "A5B117",
    "chr15_q11.2"  => "4D8A94",
    "chr15_q12"    => "9C4081",
    "chr15_q13.1"  => "20447B",
    "chr15_q13.2"  => "835A70",
    "chr15_q13.3"  => "83D424",
    "chr15_q14"    => "4DDB6E",
    "chr15_p12"    => "5AA0BD",
    "chr15_q15.1"  => "52673B",
    "chr15_q15.2"  => "38DADD",
    "chr15_q15.3"  => "A529EE",
    "chr15_q21.1"  => "264C4D",
    "chr15_q21.2"  => "AD3039",
    "chr15_q21.3"  => "69469C",
    "chr15_q22.1"  => "8BA9A3",
    "chr15_q22.2"  => "8BD8D1",
    "chr15_q22.31" => "BDAB57",
    "chr15_q22.32" => "ED6353",
    "chr15_q22.33" => "5E124C",
    "chr15_q23"    => "6ECAD4",
    "chr15_q24.1"  => "15AD13",
    "chr15_q24.2"  => "DE2271",
    "chr15_q24.3"  => "63D684",
    "chr15_q25.1"  => "DD343E",
    "chr15_p11.2"  => "9C1165",
    "chr15_q25.2"  => "500AE8",
    "chr15_q25.3"  => "D60D1C",
    "chr15_q26.1"  => "B7CEAE",
    "chr15_q26.2"  => "06E32D",
    "chr15_q26.3"  => "D50A85",
    "chr16_p13.3"  => "7DE883",
    "chr16_p13.13" => "541541",
    "chr16_p13.12" => "4AEDCE",
    "chr16_p13.11" => "F40C0C",
    "chr16_p12.3"  => "3C13C1",
    "chr16_p12.2"  => "C1DDC8",
    "chr16_p12.1"  => "BC9BC9",
    "chr16_p11.2"  => "161010",
    "chr16_p11.1"  => "176312",
    "chr16_q11.1"  => "C4E494",
    "chr16_q11.2"  => "A10A10",
    "chr16_q12.1"  => "A7BCD3",
    "chr16_q12.2"  => "D6AD6A",
    "chr16_q13"    => "D943B4",
    "chr16_q21"    => "C7C614",
    "chr16_p13.2"  => "1372AB",
    "chr16_q22.1"  => "915545",
    "chr16_q22.2"  => "CEBD68",
    "chr16_q22.3"  => "0CBCA6",
    "chr16_q23.1"  => "26C482",
    "chr16_q23.2"  => "B2B368",
    "chr16_q23.3"  => "D5A011",
    "chr16_q24.1"  => "0C916D",
    "chr16_q24.2"  => "38D38D",
    "chr16_q24.3"  => "9D0DEA",
    "chr17_p13.3"  => "8DD2DC",
    "chr17_p12"    => "964C43",
    "chr17_p11.2"  => "5CCE18",
    "chr17_p11.1"  => "692702",
    "chr17_q11.1"  => "94A791",
    "chr17_q11.2"  => "E74131",
    "chr17_q12"    => "18481B",
    "chr17_p13.2"  => "7B6104",
    "chr17_q21.1"  => "081BC4",
    "chr17_q21.2"  => "243353",
    "chr17_q21.31" => "D140B8",
    "chr17_q21.32" => "807627",
    "chr17_q21.33" => "70E4CA",
    "chr17_q22"    => "165363",
    "chr17_q23.1"  => "96D1CA",
    "chr17_q23.2"  => "E55ADB",
    "chr17_q23.3"  => "E71DB0",
    "chr17_q24.1"  => "339A86",
    "chr17_q24.2"  => "0AEE07",
    "chr17_p13.1"  => "DB203A",
    "chr17_q24.3"  => "C688B8",
    "chr17_q25.1"  => "05C44D",
    "chr17_q25.2"  => "A5C615",
    "chr17_q25.3"  => "A62DC7",
    "chr18_p11.32" => "92C22A",
    "chr18_p11.21" => "A52AA3",
    "chr18_p11.1"  => "50A365",
    "chr18_q11.1"  => "076418",
    "chr18_q11.2"  => "64EB0B",
    "chr18_q12.1"  => "3D4DD7",
    "chr18_p11.31" => "A7BB55",
    "chr18_q12.2"  => "42627B",
    "chr18_q12.3"  => "69050D",
    "chr18_q21.1"  => "276048",
    "chr18_q21.2"  => "2E1981",
    "chr18_q21.31" => "0EDE81",
    "chr18_q21.32" => "8493A0",
    "chr18_q21.33" => "126B38",
    "chr18_q22.1"  => "C05240",
    "chr18_q22.2"  => "147B30",
    "chr18_q22.3"  => "A1298A",
    "chr18_q23"    => "44DB00",
    "chr18_p11.23" => "344401",
    "chr18_p11.22" => "884D98",
    "chr19_p13.3"  => "6C6A30",
    "chr19_p13.13" => "341008",
    "chr19_p13.12" => "304BA3",
    "chr19_p13.11" => "B34A17",
    "chr19_p12"    => "E35660",
    "chr19_p11"    => "10A9BB",
    "chr19_q11"    => "9C09A0",
    "chr19_q12"    => "926A96",
    "chr19_q13.11" => "E7AE0E",
    "chr19_q13.12" => "27A105",
    "chr19_q13.13" => "0ACE0D",
    "chr19_q13.2"  => "C6585E",
    "chr19_q13.31" => "A28327",
    "chr19_q13.32" => "CACBE6",
    "chr19_q13.33" => "A58AC2",
    "chr19_q13.41" => "88E2DB",
    "chr19_q13.42" => "41CE6C",
    "chr19_q13.43" => "396801",
    "chr19_p13.2"  => "44571C",
    "chr2_p25.3"   => "78917E",
    "chr2_q12.1"   => "A09737",
    "chr2_q12.2"   => "728D8E",
    "chr2_q12.3"   => "7A8471",
    "chr2_q13"     => "8895D2",
    "chr2_q14.1"   => "113C20",
    "chr2_q14.2"   => "25A8BA",
    "chr2_p24.3"   => "5BC0D3",
    "chr2_q14.3"   => "08263B",
    "chr2_q21.1"   => "4675B3",
    "chr2_q21.2"   => "4D29DC",
    "chr2_q21.3"   => "D15323",
    "chr2_q22.1"   => "E24988",
    "chr2_q22.2"   => "A27283",
    "chr2_q22.3"   => "B7B1D8",
    "chr2_q23.1"   => "C3E17B",
    "chr2_q23.2"   => "40DB62",
    "chr2_q23.3"   => "15024C",
    "chr2_q24.1"   => "0B8475",
    "chr2_q24.2"   => "BC618B",
    "chr2_q24.3"   => "808228",
    "chr2_p24.2"   => "537BC5",
    "chr2_q31.1"   => "AA85EC",
    "chr2_q31.2"   => "3B435D",
    "chr2_q31.3"   => "C6C6D5",
    "chr2_q32.1"   => "D35E24",
    "chr2_q32.2"   => "28469C",
    "chr2_q32.3"   => "C3A7A4",
    "chr2_p24.1"   => "C85EE5",
    "chr2_q33.1"   => "5AC87D",
    "chr2_q33.2"   => "051C31",
    "chr2_q33.3"   => "4CD8B4",
    "chr2_q34"     => "D0066C",
    "chr2_q35"     => "AEBCC6",
    "chr2_q36.1"   => "9EA698",
    "chr2_q36.2"   => "AE0786",
    "chr2_q36.3"   => "3A4460",
    "chr2_q37.1"   => "9850E4",
    "chr2_q37.2"   => "3A4835",
    "chr2_q37.3"   => "AE7C9E",
    "chr2_p23.3"   => "904B32",
    "chr2_p23.2"   => "4D83DB",
    "chr2_p23.1"   => "51440A",
    "chr2_p22.3"   => "D0A201",
    "chr2_p22.2"   => "EB6D3E",
    "chr2_p25.2"   => "B74C75",
    "chr2_p22.1"   => "26EE55",
    "chr2_p21"     => "96A6EE",
    "chr2_p16.3"   => "1B62A6",
    "chr2_p16.2"   => "B84A3C",
    "chr2_p16.1"   => "71143E",
    "chr2_p15"     => "44BCD9",
    "chr2_p14"     => "D96A65",
    "chr2_p13.3"   => "59C400",
    "chr2_p25.1"   => "C83C3C",
    "chr2_p13.2"   => "3BA900",
    "chr2_p13.1"   => "97EC3A",
    "chr2_p12"     => "87982A",
    "chr2_p11.2"   => "0EA7AA",
    "chr2_p11.1"   => "BB3D52",
    "chr2_q11.1"   => "A64664",
    "chr2_q11.2"   => "26805C",
    "chr20_p13"    => "24842C",
    "chr20_p12.1"  => "6118D5",
    "chr20_p11.23" => "3D0396",
    "chr20_p11.22" => "23349A",
    "chr20_p11.21" => "75DC41",
    "chr20_p11.1"  => "C2248E",
    "chr20_q11.1"  => "C3452D",
    "chr20_q11.21" => "7CDB79",
    "chr20_q11.22" => "7EC23C",
    "chr20_q11.23" => "B76C4E",
    "chr20_q12"    => "38C579",
    "chr20_q13.11" => "8E5C7A",
    "chr20_q13.12" => "A8A2DD",
    "chr20_q13.13" => "831496",
    "chr20_q13.2"  => "81E3A8",
    "chr20_p12.3"  => "15CCA0",
    "chr20_q13.31" => "E51AE4",
    "chr20_q13.32" => "9C0CC6",
    "chr20_q13.33" => "D34235",
    "chr20_p12.2"  => "897940",
    "chr21_p13"    => "000E05",
    "chr21_p11.1"  => "ED5C37",
    "chr21_q11.1"  => "3153A1",
    "chr21_q11.2"  => "C540D3",
    "chr21_q21.1"  => "1B29E9",
    "chr21_q21.2"  => "88DE1B",
    "chr21_q21.3"  => "C439D0",
    "chr21_p12"    => "D8D6A4",
    "chr21_q22.11" => "E99DD1",
    "chr21_q22.12" => "2BDC62",
    "chr21_q22.13" => "B3C700",
    "chr21_q22.2"  => "E85422",
    "chr21_q22.3"  => "9A4D52",
    "chr21_p11.2"  => "A8B0B8",
    "chr22_p13"    => "732D6C",
    "chr22_q11.1"  => "13286A",
    "chr22_q11.21" => "61921A",
    "chr22_q11.22" => "1D0108",
    "chr22_q11.23" => "8A1431",
    "chr22_q12.1"  => "3017AA",
    "chr22_q12.2"  => "A624D6",
    "chr22_p12"    => "EC72E2",
    "chr22_q12.3"  => "DE6B5D",
    "chr22_q13.1"  => "0913E4",
    "chr22_q13.2"  => "1297EC",
    "chr22_q13.31" => "8876A1",
    "chr22_q13.32" => "1D8204",
    "chr22_q13.33" => "A05880",
    "chr22_p11.2"  => "407A35",
    "chr22_p11.1"  => "93C176",
    "chr3_p26.3"   => "C4E482",
    "chr3_q12.2"   => "EC8587",
    "chr3_q12.3"   => "625811",
    "chr3_q13.11"  => "433AC5",
    "chr3_q13.12"  => "9B4B73",
    "chr3_q13.13"  => "477750",
    "chr3_q13.2"   => "DDCA3D",
    "chr3_q13.31"  => "A61BB2",
    "chr3_p25.2"   => "5388C7",
    "chr3_q13.32"  => "48419B",
    "chr3_q13.33"  => "DBE11D",
    "chr3_q21.1"   => "98183C",
    "chr3_q21.2"   => "1202BC",
    "chr3_q21.3"   => "2E6A7A",
    "chr3_q22.1"   => "04238B",
    "chr3_p25.1"   => "E8D4B8",
    "chr3_q22.2"   => "349EBA",
    "chr3_q22.3"   => "DBCE84",
    "chr3_q23"     => "3D77B3",
    "chr3_q24"     => "829AEA",
    "chr3_q25.1"   => "14D4DD",
    "chr3_q25.2"   => "D96156",
    "chr3_q25.31"  => "9AE67E",
    "chr3_q25.32"  => "CD4767",
    "chr3_q25.33"  => "63290A",
    "chr3_q26.1"   => "ABDCE9",
    "chr3_p24.3"   => "7BC4D3",
    "chr3_q26.2"   => "031AE6",
    "chr3_q26.31"  => "B2396D",
    "chr3_q26.32"  => "567304",
    "chr3_q26.33"  => "6040D4",
    "chr3_q27.1"   => "D5B70C",
    "chr3_q27.2"   => "36B742",
    "chr3_q27.3"   => "D27DB9",
    "chr3_q28"     => "99BDA2",
    "chr3_q29"     => "AB3439",
    "chr3_p24.2"   => "502760",
    "chr3_p26.2"   => "25E4C2",
    "chr3_p24.1"   => "3361CE",
    "chr3_p23"     => "3C0369",
    "chr3_p22.3"   => "A3CE2B",
    "chr3_p22.2"   => "CB8DA3",
    "chr3_p22.1"   => "9A2D12",
    "chr3_p26.1"   => "C537AD",
    "chr3_p21.33"  => "525B44",
    "chr3_p21.32"  => "314780",
    "chr3_p21.31"  => "67B60E",
    "chr3_p21.2"   => "5BA12D",
    "chr3_p21.1"   => "BA647B",
    "chr3_p14.3"   => "A2E50C",
    "chr3_p14.2"   => "A332CC",
    "chr3_p14.1"   => "AD4692",
    "chr3_p13"     => "A23383",
    "chr3_p12.3"   => "73B0E0",
    "chr3_p12.2"   => "B12911",
    "chr3_p25.3"   => "472802",
    "chr3_p12.1"   => "E16AB6",
    "chr3_p11.2"   => "4E89B5",
    "chr3_p11.1"   => "5895A8",
    "chr3_q11.1"   => "547092",
    "chr3_q11.2"   => "027558",
    "chr3_q12.1"   => "D43B3C",
    "chr4_p16.3"   => "E400E7",
    "chr4_q24"     => "E36A20",
    "chr4_q25"     => "BA7168",
    "chr4_p15.33"  => "9641B7",
    "chr4_q26"     => "33865C",
    "chr4_q27"     => "09446D",
    "chr4_q28.1"   => "600A72",
    "chr4_q28.2"   => "303EAB",
    "chr4_q28.3"   => "A4A3E0",
    "chr4_q31.1"   => "6319E5",
    "chr4_q31.21"  => "86414D",
    "chr4_q31.22"  => "83005B",
    "chr4_q31.23"  => "E453D7",
    "chr4_q31.3"   => "BE0E71",
    "chr4_p15.32"  => "0AA167",
    "chr4_q32.1"   => "47D6C1",
    "chr4_q32.2"   => "5D08EC",
    "chr4_q32.3"   => "309B05",
    "chr4_q33"     => "586597",
    "chr4_q34.1"   => "E7139E",
    "chr4_q34.2"   => "6BC880",
    "chr4_q34.3"   => "BEEAAD",
    "chr4_q35.1"   => "50598D",
    "chr4_p15.31"  => "042832",
    "chr4_q35.2"   => "87A5CB",
    "chr4_p15.2"   => "46103C",
    "chr4_p15.1"   => "559273",
    "chr4_p14"     => "BDEA54",
    "chr4_p13"     => "61E462",
    "chr4_p12"     => "A51949",
    "chr4_p16.2"   => "A4D006",
    "chr4_p11"     => "B155A2",
    "chr4_q11"     => "590601",
    "chr4_q12"     => "CE2A44",
    "chr4_p16.1"   => "910676",
    "chr4_q13.1"   => "EB11F2",
    "chr4_q13.2"   => "57007B",
    "chr4_q13.3"   => "6B12A8",
    "chr4_q21.1"   => "49E1D7",
    "chr4_q21.21"  => "21A381",
    "chr4_q21.22"  => "0D8D55",
    "chr4_q21.23"  => "D6E146",
    "chr4_q21.3"   => "252D8A",
    "chr4_q22.1"   => "E484EC",
    "chr4_q22.2"   => "C1A3C7",
    "chr4_q22.3"   => "54AAA1",
    "chr4_q23"     => "D6E88E",
    "chr5_p15.33"  => "4B8700",
    "chr5_q21.2"   => "2278C5",
    "chr5_q21.3"   => "5E2455",
    "chr5_q22.1"   => "B6B9B3",
    "chr5_q22.2"   => "6B909B",
    "chr5_q22.3"   => "E3B0C7",
    "chr5_q23.1"   => "BBDD53",
    "chr5_q23.2"   => "C37EE4",
    "chr5_q23.3"   => "6ABCB4",
    "chr5_q31.1"   => "0CD4B3",
    "chr5_q31.2"   => "DCE918",
    "chr5_q31.3"   => "D1B7E1",
    "chr5_q32"     => "4E3396",
    "chr5_q33.1"   => "423BE3",
    "chr5_p15.1"   => "223EA2",
    "chr5_q33.2"   => "E9D4E9",
    "chr5_q33.3"   => "228A54",
    "chr5_q34"     => "210DB7",
    "chr5_q35.1"   => "08623B",
    "chr5_q35.2"   => "2000A5",
    "chr5_q35.3"   => "53E9C0",
    "chr5_p14.3"   => "905CCD",
    "chr5_p14.2"   => "52BB21",
    "chr5_p14.1"   => "31CC4D",
    "chr5_p13.3"   => "4083B2",
    "chr5_p13.2"   => "8BAC2D",
    "chr5_p13.1"   => "02B240",
    "chr5_p12"     => "89ABD1",
    "chr5_p15.32"  => "C99219",
    "chr5_p11"     => "8B22CA",
    "chr5_q11.1"   => "D9C381",
    "chr5_q11.2"   => "1429BC",
    "chr5_q12.1"   => "DEEEBD",
    "chr5_q12.2"   => "9C9597",
    "chr5_q12.3"   => "0E1823",
    "chr5_p15.31"  => "34009A",
    "chr5_q13.1"   => "BE59E1",
    "chr5_q13.2"   => "4AAA29",
    "chr5_q13.3"   => "581453",
    "chr5_q14.1"   => "84BE9C",
    "chr5_q14.2"   => "1075E8",
    "chr5_q14.3"   => "D268A9",
    "chr5_q15"     => "588661",
    "chr5_q21.1"   => "88819A",
    "chr5_p15.2"   => "978901",
    "chr6_p25.3"   => "2E9D1E",
    "chr6_q16.3"   => "EB80B5",
    "chr6_q21"     => "EE0509",
    "chr6_p24.2"   => "CE7A04",
    "chr6_q22.1"   => "258E12",
    "chr6_p24.1"   => "428E8D",
    "chr6_q22.2"   => "C12D4B",
    "chr6_q22.31"  => "A47DBC",
    "chr6_q22.32"  => "C796B6",
    "chr6_q22.33"  => "45C431",
    "chr6_q23.1"   => "672C95",
    "chr6_q23.2"   => "C7780E",
    "chr6_p23"     => "5EEA6D",
    "chr6_q23.3"   => "92E448",
    "chr6_q24.1"   => "364A22",
    "chr6_q24.2"   => "4A0865",
    "chr6_q24.3"   => "131591",
    "chr6_q25.1"   => "93E047",
    "chr6_q25.2"   => "B23E26",
    "chr6_p22.3"   => "98825E",
    "chr6_q25.3"   => "89BA1C",
    "chr6_q26"     => "D90D69",
    "chr6_q27"     => "1AD03C",
    "chr6_p25.2"   => "C42E3D",
    "chr6_p22.2"   => "746744",
    "chr6_p22.1"   => "2CEA8A",
    "chr6_p21.33"  => "E3DEB3",
    "chr6_p21.32"  => "A13181",
    "chr6_p21.31"  => "05D060",
    "chr6_p21.2"   => "84E286",
    "chr6_p21.1"   => "6C9E21",
    "chr6_p25.1"   => "4B94D9",
    "chr6_p12.3"   => "962E83",
    "chr6_p12.2"   => "D15563",
    "chr6_p12.1"   => "9AAA09",
    "chr6_p11.2"   => "B080E5",
    "chr6_p11.1"   => "1ED305",
    "chr6_q11.1"   => "1CA3A7",
    "chr6_q11.2"   => "174200",
    "chr6_q12"     => "1905B3",
    "chr6_p24.3"   => "EDECEB",
    "chr6_q13"     => "924DC0",
    "chr6_q14.1"   => "4E5190",
    "chr6_q14.2"   => "146141",
    "chr6_q14.3"   => "60C0B7",
    "chr6_q15"     => "45224E",
    "chr6_q16.1"   => "C1179E",
    "chr6_q16.2"   => "BB4B9C",
    "chr7_p22.3"   => "23325D",
    "chr7_q22.2"   => "D875E5",
    "chr7_q22.3"   => "861473",
    "chr7_q31.1"   => "440AAC",
    "chr7_q31.2"   => "A878AC",
    "chr7_q31.31"  => "52D35A",
    "chr7_q31.32"  => "90DA6B",
    "chr7_q31.33"  => "EBAE3B",
    "chr7_q32.1"   => "933638",
    "chr7_q32.2"   => "219D04",
    "chr7_q32.3"   => "9D35C0",
    "chr7_p21.2"   => "072A7C",
    "chr7_q33"     => "928D27",
    "chr7_q34"     => "84B855",
    "chr7_q35"     => "BC636B",
    "chr7_q36.1"   => "433A55",
    "chr7_p21.1"   => "33465B",
    "chr7_q36.2"   => "14C3B2",
    "chr7_q36.3"   => "1805D4",
    "chr7_p15.3"   => "E17BD2",
    "chr7_p15.2"   => "629298",
    "chr7_p22.2"   => "8A475B",
    "chr7_p15.1"   => "05BED6",
    "chr7_p14.3"   => "632630",
    "chr7_p14.2"   => "705494",
    "chr7_p14.1"   => "16A2DE",
    "chr7_p22.1"   => "7D4884",
    "chr7_p13"     => "86E849",
    "chr7_p12.3"   => "1C116E",
    "chr7_p12.2"   => "762EE9",
    "chr7_p12.1"   => "644BB8",
    "chr7_p11.2"   => "641717",
    "chr7_p11.1"   => "30E932",
    "chr7_q11.1"   => "C0BAA7",
    "chr7_q11.21"  => "BA7DCD",
    "chr7_q11.22"  => "C21A11",
    "chr7_p21.3"   => "8C3B4B",
    "chr7_q11.23"  => "AEB81B",
    "chr7_q21.11"  => "281332",
    "chr7_q21.12"  => "AC6B87",
    "chr7_q21.13"  => "43D7BB",
    "chr7_q21.2"   => "6C7DD5",
    "chr7_q21.3"   => "8D1252",
    "chr7_q22.1"   => "4EAB0A",
    "chr8_p23.3"   => "2188D2",
    "chr8_q22.3"   => "B99DDD",
    "chr8_q23.1"   => "A74882",
    "chr8_q23.2"   => "DE1601",
    "chr8_q23.3"   => "6751C9",
    "chr8_q24.11"  => "CACE64",
    "chr8_q24.12"  => "D2A8BC",
    "chr8_q24.13"  => "5DA2EA",
    "chr8_p22"     => "DE223B",
    "chr8_q24.21"  => "ACABD2",
    "chr8_q24.22"  => "EC0308",
    "chr8_q24.23"  => "A614BD",
    "chr8_q24.3"   => "D383AC",
    "chr8_p21.3"   => "881CC5",
    "chr8_p23.2"   => "03A142",
    "chr8_p21.2"   => "A5C188",
    "chr8_p21.1"   => "0D13CB",
    "chr8_p12"     => "6A4AAC",
    "chr8_p11.23"  => "4629AA",
    "chr8_p11.22"  => "76D8C2",
    "chr8_p11.21"  => "3153D3",
    "chr8_p11.1"   => "BAA236",
    "chr8_q11.1"   => "54906E",
    "chr8_q11.21"  => "440631",
    "chr8_q11.22"  => "18A0D8",
    "chr8_q11.23"  => "AE2393",
    "chr8_q12.1"   => "8D6775",
    "chr8_q12.2"   => "D2D1D5",
    "chr8_p23.1"   => "A9361B",
    "chr8_q12.3"   => "B1C7B4",
    "chr8_q13.1"   => "291548",
    "chr8_q13.2"   => "798410",
    "chr8_q13.3"   => "75A510",
    "chr8_q21.11"  => "DA1575",
    "chr8_q21.12"  => "04A390",
    "chr8_q21.13"  => "910172",
    "chr8_q21.2"   => "3BB031",
    "chr8_q21.3"   => "AD2971",
    "chr8_q22.1"   => "906727",
    "chr8_q22.2"   => "8D02A2",
    "chr9_p24.3"   => "CC2801",
    "chr9_q31.2"   => "E48D26",
    "chr9_q31.3"   => "64C0EB",
    "chr9_q32"     => "BAA9CA",
    "chr9_q33.1"   => "53E6C9",
    "chr9_q33.2"   => "AC5047",
    "chr9_q33.3"   => "2C38E9",
    "chr9_q34.11"  => "607B48",
    "chr9_q34.12"  => "9CE298",
    "chr9_q34.13"  => "A80C97",
    "chr9_q34.2"   => "3465BE",
    "chr9_q34.3"   => "633E29",
    "chr9_p22.3"   => "821AA6",
    "chr9_p22.2"   => "548C8E",
    "chr9_p22.1"   => "44BCB9",
    "chr9_p21.3"   => "234034",
    "chr9_p24.2"   => "8832BA",
    "chr9_p21.2"   => "90241B",
    "chr9_p21.1"   => "264E59",
    "chr9_p13.3"   => "619BBC",
    "chr9_p13.2"   => "4B0E7A",
    "chr9_p13.1"   => "CB1850",
    "chr9_p12"     => "9B388B",
    "chr9_p11.2"   => "E588A0",
    "chr9_p11.1"   => "32EE14",
    "chr9_p24.1"   => "7CA25A",
    "chr9_q11"     => "36C8BE",
    "chr9_q12"     => "24054B",
    "chr9_q13"     => "3CD970",
    "chr9_q21.11"  => "DEDEC3",
    "chr9_q21.12"  => "19D2C0",
    "chr9_q21.13"  => "E20E4E",
    "chr9_q21.2"   => "656A2A",
    "chr9_q21.31"  => "790DE7",
    "chr9_q21.32"  => "941429",
    "chr9_q21.33"  => "1DC414",
    "chr9_q22.1"   => "082248",
    "chr9_q22.2"   => "4952B2",
    "chr9_p23"     => "0D3476",
    "chr9_q22.31"  => "B5332D",
    "chr9_q22.32"  => "6BC6D1",
    "chr9_q22.33"  => "25B538",
    "chr9_q31.1"   => "D591E2",
    "chrX_p22.33"  => "5B5794",
    "chrX_q22.2"   => "9C4220",
    "chrX_q22.3"   => "E9D469",
    "chrX_q23"     => "A3C77B",
    "chrX_q24"     => "EC86D5",
    "chrX_q25"     => "690488",
    "chrX_q26.1"   => "C6880C",
    "chrX_q26.2"   => "EB2916",
    "chrX_q26.3"   => "4B86DD",
    "chrX_q27.1"   => "862A24",
    "chrX_q27.2"   => "D92C0B",
    "chrX_q27.3"   => "B566E7",
    "chrX_q28"     => "4E0734",
    "chrX_p22.13"  => "990B6E",
    "chrX_p22.12"  => "C18BC1",
    "chrX_p22.11"  => "B09A1E",
    "chrX_p21.3"   => "5CB646",
    "chrX_p21.2"   => "25DA8E",
    "chrX_p21.1"   => "C4D165",
    "chrX_p11.4"   => "E68959",
    "chrX_p22.32"  => "476418",
    "chrX_p11.3"   => "2D85A4",
    "chrX_p11.23"  => "3E8129",
    "chrX_p11.22"  => "7803EB",
    "chrX_p11.21"  => "13D5EA",
    "chrX_p11.1"   => "4154C0",
    "chrX_p22.31"  => "648AC9",
    "chrX_q11.1"   => "BA9517",
    "chrX_q11.2"   => "0889E1",
    "chrX_q12"     => "48ED12",
    "chrX_q13.1"   => "10E3E4",
    "chrX_q13.2"   => "69BB08",
    "chrX_q13.3"   => "21AB7A",
    "chrX_q21.1"   => "BA4B8D",
    "chrX_q21.2"   => "974826",
    "chrX_q21.31"  => "21D5D8",
    "chrX_q21.32"  => "49435D",
    "chrX_p22.2"   => "A1826B",
    "chrX_q21.33"  => "55968D",
    "chrX_q22.1"   => "C2A0A6",
    "chrY_p11.32"  => "7A5957",
    "chrY_p11.1"   => "A49189",
    "chrY_q11.1"   => "C1047A",
    "chrY_q11.21"  => "AD3168",
    "chrY_q11.221" => "174B15",
    "chrY_p11.31"  => "80693E",
    "chrY_q11.222" => "550689",
    "chrY_q11.223" => "828872",
    "chrY_q11.23"  => "835393",
    "chrY_q12"     => "4A6014",
    "chrY_p11.2"   => "9EB690"
  );
  return ( %colors );
}

1;
