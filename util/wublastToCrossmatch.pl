#!/usr/local/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) wublastToCrossmatch.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a template for generic perl scripts.  It
##      includes an area for POD documentation (ie. perldoc this-file )
##      and a generic way to store startup parameters in a file
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
#     $Log$
#
###############################################################################
#
# To Do:
#

=head1 NAME

wublastToCrossmatch.pl - Reformat/rescore wublast results in a cm-like manner

=head1 SYNOPSIS

  wublastToCrossmatch.pl [-version]
                         [-matrix="Matrices/crossmatch/14p41g.matrix"]
                         [-gap_init=-40]
                         [-gap_ext=-15]
                         [-minscore=##]
                         <wublast result file>

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2013 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/..";
use Matrix;
use WUBlastSearchEngine;
use SearchResult;
use SearchResultCollection;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name:  $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name:  $';
my $CVSIdTag   =
    '$Id: wublastToCrossmatch.pl,v 1.9 2017/02/01 21:01:58 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',      # print out the version and exit
                    '-matrix=s',
                    '-gap_init=s',
                    '-gap_ext=s',
                    '-minscore=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

my $defaultMatrix = "$FindBin::Bin/../Matrices/crossmatch/14p41g.matrix";
my $matrixName    = $defaultMatrix;
$matrixName = $options{'matrix'} if ( -s $options{'matrix'} );
print STDERR "Using matrix: $matrixName\n";

my $defaultGapInit = -40;
my $defaultGapExt  = -15;
my $gapInit        = $defaultGapInit;
my $gapExt         = $defaultGapExt;
$gapInit = $options{'gap_init'} if ( $options{'gap_init'} );
$gapExt  = $options{'gap_ext'}  if ( $options{'gap_ext'} );
print STDERR "Using gap_open/gap_ext: $gapInit/$gapExt\n";

unless ( -s $ARGV[ 0 ] ) {
  print "\n\nMissing wublast results input file!\n\n";
  usage();
}

#
# Parse matrix file
#
my $matrix = Matrix->new( fileName => $matrixName );

#
# Parse WUBlast results
#
print STDERR "Parsing $ARGV[0] ...\n";
my $searchResultsCollection =
    WUBlastSearchEngine::parseOutput( searchOutput => $ARGV[ 0 ] );
print STDERR "Total results read: " . $searchResultsCollection->size() . "\n";

#
# Create output
#
for ( my $i = 0 ; $i < $searchResultsCollection->size() ; $i++ ) {
  my $result = $searchResultsCollection->get( $i );

  my $rawScore = $result->getScore();
  my (
       $adjScore,        $kimura,                 $cpgsites,
       $percIns,         $percDel,                $positionScores,
       $xdrop_fragments, $wellCharacterizedBases, $trans,
       $transv
      )
      = $result->rescoreAlignment(
                                   scoreMatrix         => $matrix,
                                   gapOpenPenalty      => $gapInit,
                                   gapExtensionPenalty => $gapExt,
                                   complexityAdjust    => 1
      );
  next if ( $options{'minscore'} && $adjScore < $options{'minscore'} );
  $result->setScore( $adjScore );
  print "" . $result->toStringFormatted( SearchResult::AlignWithQuerySeq );
  print "## WUBlast Score: $rawScore\n";
  print "## Well charactized bases: $wellCharacterizedBases\n";
  print "## transI/transV: $trans/$transv\n";
  print "## Kimura Div: " . sprintf( "%0.2f", $kimura ) . "\n";
  print "## CpG Sites: $cpgsites\n";
  print "## Perc Ins/Del: "
      . sprintf( "%0.2f", $percIns ) . "/"
      . sprintf( "%0.2f", $percDel ) . "\n";
  print "## Position scores: " . join( ",", @{$positionScores} ) . "\n\n\n";
}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Methods with the prefix "_" are conventionally considered
##      private.  This bit-o-documentation is not formatted to
##      print out when perldoc is run on this file.
##
##-------------------------------------------------------------------------##
sub _privateMethod {
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " . @{ [ %parameters ] } . "): Called\n" if ( $DEBUG );

}

##-------------------------------------------------------------------------##

=head2 publicMethod()

  Use:  my $retVal = publicMethod( $parameter1 => value, 
                                   $parameter2 => value );

    $parameter1:   A generic scalar parameter
    $parameter2:   A generic scalar parameter

  $retVal contains the scalar result of this subroutine.  This
  is a public function and this documentation will print out
  when perldoc is run on this file.

=cut

##-------------------------------------------------------------------------##·
sub publicMethod {
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " . @{ [ %parameters ] } . "): Called\n" if ( $DEBUG );
}

1;

