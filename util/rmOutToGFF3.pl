#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) rmOutToGFF3.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to assist to convert old RepeatMasker *.out files
##      to version 3 gff files.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2008 Developed by
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
#     $Log: rmOutToGFF3.pl,v $
#     Revision 1.35  2017/02/01 21:01:58  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

rmOutToGFF3.pl - Convert RepeatMasker OUT files to version 3 GFF files

=head1 SYNOPSIS

  rmOutToGFF3.pl [-version] *.out > new.gff

=head1 DESCRIPTION

  A utility script to convert old RepeatMasker *.out files to 
  the new GFF version 3 standard.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

ReapeatMasker

=head1 COPYRIGHT

Copyright 2008 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Getopt::Long;
use Data::Dumper;
use CrossmatchSearchEngine;
use File::Basename;

#
# Version
#
my $Version = 0.1;
my $DEBUG   = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
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
  exit( 1 );
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

usage() if ( !$ARGV[ 0 ] );

my $annotationFile = $ARGV[ 0 ];

#
# Open up a search results object
#
my $searchResults =
    CrossmatchSearchEngine::parseOutput( searchOutput => $annotationFile );

#
# Read in annotations and throw away the rest
#
print "##gff-version 3\n";
my $currentQueryName;
for ( my $i = 0 ; $i < $searchResults->size() ; $i++ ) {
  my $result = $searchResults->get( $i );

  # First annotation of a region
  if ( $result->getQueryName() ne $currentQueryName ) {
    $currentQueryName = $result->getQueryName();
    print "##sequence-region $currentQueryName 1 "
        . ( $result->getQueryRemaining() + $result->getQueryEnd() ) . "\n";
  }

  # FORMAT:
  #   ##gff-version   3
  #   ##sequence-region   ctg123 1 1497228
  #   SeqID:     QueryName
  #   Source:    Constant - "RepeatMasker"
  #   Type:      similarity => dispersed_repeat
  #   Start:     Query Start
  #   End:       Query End
  #   Score:     New!
  #   Strand:    "+" or "-"
  #   Phase:     0
  #   Attributes: Target=FAM 24 180
  print "" . $currentQueryName . "\t";
  print "RepeatMasker\t";
  print "dispersed_repeat\t";
  print "" . $result->getQueryStart() . "\t";
  print "" . $result->getQueryEnd() . "\t";
  print "" . $result->getScore() . "\t";
  if ( $result->getOrientation() eq "C" ) {
    print "-\t";
  }
  else {
    print "+\t";
  }
  print ".\t";
  print "Target="
      . $result->getSubjName() . " "
      . $result->getSubjStart() . " "
      . $result->getSubjEnd() . "\n";

}

