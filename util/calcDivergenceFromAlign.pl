#!/usr/local/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) calcDivergenceFromAlign.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to calculate a new divergence measure on the
##      RM alignment files.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2009 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
#
# ChangeLog
#
#     $Log: calcDivergenceFromAlign.pl,v $
#     Revision 1.27  2017/02/01 21:01:56  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

calcDivergenceFromAlign.pl - Recalculate/Summarize the divergences in an align file.

=head1 SYNOPSIS

  calcDivergenceFromAlign.pl [-version] [-s <summary_file>] [-noCpGMod]
                             [-a <new_align_file>] *.align[.gz]

=head1 DESCRIPTION

  A utility script to calculate a new divergence measure on the
  RM alignment files.  Currently we only calculate the Kimura 2-Parameter
  divergence metric.  

  Treat "CG" dinucleotide sites in the consensus sequence as follows:
  Two transition mutations are counted as a single transition, one
  transition is counted as 1/10 of a standard transition, and 
  transversions are counted normally (as the would outside of a CpG
  site).  This modification to the Kimura 2 parameter model accounts
  for the extremely high rate of mutations in at a CpG locus.  


The options are:

=over 4

=item -version

Displays the version of the program

=item -noCpGMod

Do not modify the transition counts at CpG sites. 

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
use FileHandle;

## TODO: Remove this
use lib "/home/rhubley/projects/RepeatMasker";
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use SearchResult;
use CrossmatchSearchEngine;

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
    '$Id: calcDivergenceFromAlign.pl,v 1.27 2017/02/01 21:01:56 rhubley Exp $';
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
                    '-version',    # print out the version and exit
                    '-noCpGMod',
                    '-a=s',
                    '-s=s'
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

if ( !$options{'s'} && !$options{'a'} ) {
  print
"\n\nError: One or more of the options '-a' or '-s' must be supplied!\n\n";
  usage();
}

my $alignFile      = $ARGV[ 0 ];
my $maxDiv         = 70;
my $cntAlign       = 0;
my %repeatMuts     = ();
my %classDivWCLen  = ();
my $prevQueryName  = "";
my $prevQueryBegin = "";
my $prevQueryEnd   = "";
my $prevDiv        = "";
my $prevClass      = "";

my $searchResultsFH = new FileHandle;

if ( $alignFile =~ /.+\.gz/ ) {
  open $searchResultsFH, "gunzip -c $alignFile|"
      or die
      "RepeatLandscape: Could not open gunzip for reading $alignFile: $!\n";
}
else {
  open $searchResultsFH, "<$alignFile"
      or die "RepeatLandscape: Could not open $alignFile for reading: $!\n";
}

if ( $options{'s'} ) {
  open SOUT, ">$options{'s'}"
      or die "Error: Could not open $options{'s'} for writing!\n";
  if ( $options{'noCpGMod'} ) {
    print SOUT "Jukes/Cantor and Kimura subsitution levels\n";
    print SOUT "==========================================\n";
  }
  else {
    print SOUT
        "Jukes/Cantor and Kimura subsitution levels adjusted for CpG sites\n";
    print SOUT
        "=================================================================\n";
  }
  print SOUT "File: " . $alignFile . "\n";
}

#
# Process the alignment file
#
my $outAlign = 0;
if ( $options{'a'} ) {
  open COUT, ">$options{'a'}"
      or die "Could not open $options{'a'} for writing!\n";
  $outAlign = 1;
}

CrossmatchSearchEngine::parseOutput( searchOutput => $searchResultsFH,
                                     callback     => \&processAlignment );

if ( $options{'a'} ) {
  close COUT;
}

if ( $options{'s'} ) {

  print SOUT "Weighted average Kimura divergence for each repeat family\n";
  print SOUT "Class\tRepeat\tabsLen\twellCharLen\tKimura%\n";
  print SOUT "-----\t------\t------\t-----------\t-------\n";
  foreach my $class ( sort keys %repeatMuts ) {
    foreach my $id ( sort keys %{ $repeatMuts{$class} } ) {
      my $kimura = 100;
      if ( $repeatMuts{$class}->{$id}->{'wellCharLen'} > 0 ) {
        $kimura = sprintf( "%4.2f",
                           $repeatMuts{$class}->{$id}->{'sumdiv'} /
                               $repeatMuts{$class}->{$id}->{'wellCharLen'} );

        $kimura = $maxDiv if ( $kimura > $maxDiv );
      }

      if ( $class =~ /Simple|Low_complexity|ARTEFACT/ ) {
        print SOUT "$class\t$id\t"
            . $repeatMuts{$class}->{$id}->{'absLen'} . "\t"
            . $repeatMuts{$class}->{$id}->{'wellCharLen'}
            . "\t----\n";
      }
      else {
        print SOUT "$class\t$id\t"
            . $repeatMuts{$class}->{$id}->{'absLen'} . "\t"
            . $repeatMuts{$class}->{$id}->{'wellCharLen'}
            . "\t$kimura\n";
      }
    }
  }
  print SOUT "\n\n";

  print SOUT "Coverage for each repeat class and divergence (Kimura)\n";
  print SOUT "Div ";
  foreach my $class ( sort keys %repeatMuts ) {
    print SOUT "$class ";
  }
  print SOUT "\n";

  my $j = 0;
  while ( $j <= $maxDiv ) {
    print SOUT "$j ";
    foreach my $class ( sort keys %repeatMuts ) {
      my $label = "$class $j";
      $classDivWCLen{$label} = 0 unless $classDivWCLen{$label};
      print SOUT "$classDivWCLen{$label} ";
    }
    print SOUT "\n";
    ++$j;
  }
  close SOUT;
}

exit;

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my processAlignment( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##
##-------------------------------------------------------------------------##
sub processAlignment {
  my $result = shift;

  return if ( !$result );

  my $hitname;
  my $class;
  my $subjName = $result->getSubjName();
  if ( $subjName =~ /(\S+)\#(\S+)/ ) {
    $hitname = $1;
    $class   = $2;
  }
  else {
    $hitname = $subjName;
    $class   = $result->getSubjType();
  }

  my $seqName    = $result->getQueryName();
  my $queryStart = $result->getQueryStart();
  my $queryEnd   = $result->getQueryEnd();

  # Simple repeats, low complexity and artefacts should not be counted
  #if ( $class =~ /Simple|Low_complexity|ARTEFACT/ )
  #{
  #  if ( $outAlign )
  #  {
  #    print COUT ""
  #      . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
  #  }
  #  return;
  #}

  print STDERR "." if ( $cntAlign++ % 1000 == 0 );

  my ( $div, $transi, $transv, $wellCharBases, $numCpGs );

  my $alen = $queryEnd - $queryStart + 1;
  $wellCharBases = $alen - int( $alen * ( $result->getPctInsert() / 100 ) );

  if ( $class =~ /Simple|Low_complexity|ARTEFACT/ ) {
    $div     = 100;
    $hitname = "combined";
  }
  else {

    # Obtain divergence from modern *.align files directly
    $div = $result->getPctKimuraDiverge();

    if ( $div eq "" ) {

      # Calculate divergence on the fly
      ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
          $result->calcKimuraDivergence( divCpGMod => 1 );
      $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    }
  }

  if ( $prevQueryName eq $seqName ) {
    if ( $prevQueryEnd > $queryStart ) {

      # Overlap
      my $overlapAbsLen = $prevQueryEnd - $queryStart + 1;
      if ( $prevQueryEnd >= $queryEnd ) {
        if ( $outAlign ) {
          print COUT ""
              . $result->toStringFormatted( SearchResult::AlignWithQuerySeq )
              . "\n";
        }
        return;
      }
      if ( $div > $prevDiv ) {

        # Previous gets overlap bases - subtract overlap from this hit
        $wellCharBases -= $overlapAbsLen;
        $alen          -= $overlapAbsLen;
      }
      else {

        # Current gets overlap bases - subtract overlap from previous
        my $key = "$prevClass $prevDiv";
        $classDivWCLen{$key}                        -= $overlapAbsLen;
        $repeatMuts{$class}->{$hitname}->{'sumdiv'} -=
            $prevDiv * $overlapAbsLen;
        $repeatMuts{$class}->{$hitname}->{'wellCharLen'} -= $overlapAbsLen;
        $repeatMuts{$class}->{$hitname}->{'absLen'}      -= $overlapAbsLen;
      }
    }
  }
  $prevQueryName  = $seqName;
  $prevDiv        = $div;
  $prevQueryBegin = $queryStart;
  $prevQueryEnd   = $queryEnd;
  $prevClass      = $class;

  $repeatMuts{$class}->{$hitname}->{'sumdiv'}      += $div * $wellCharBases;
  $repeatMuts{$class}->{$hitname}->{'wellCharLen'} += $wellCharBases;
  $repeatMuts{$class}->{$hitname}->{'absLen'}      += $alen;
  $div = int( $div );
  my $key = "$class $div";
  $classDivWCLen{$key} += $wellCharBases;

  if ( $outAlign ) {
    print COUT ""
        . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
  }

}

1;
