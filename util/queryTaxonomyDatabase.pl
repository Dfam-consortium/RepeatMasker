#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) queryTaxonomyDatabase.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to assist in querying the RepeatMasker
##      taxonomy database.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
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
#     $Log: queryTaxonomyDatabase.pl,v $
#     Revision 1.77  2017/02/01 21:01:57  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

queryTaxonomyDatabase.pl - Query the RepeatMasker taxonomy database.

=head1 SYNOPSIS

  queryTaxonomyDatabase.pl [-version] -species "species name"
                                      [-isa "species name"]
                                      [-taxDBFile <RepeatMasker Taxonomy Database File>]

=head1 DESCRIPTION

  A utility script to query the RepeatMasker Taxonomy database.
  See Taxonomy.pm for details on the database build process.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=over 4

=item -species "species name"

The full name ( case insensitive ) of the species you would like
to search for in the database.

=back

=over 4

=item -isa "species name"

The full name of a species group.  The utility will determine
if the species is a member of this species group.

=back

=over 4

=item -taxDBFile <RepeatMasker Taxonomy Database File>

The full path name of the RepeatMasker taxonomy database file.  If not
specified the program defaults to looking for the file in relative
to the scripts directory as:

           ../taxonomy.dat

=back

=head1 SEE ALSO

Taxonomy.pm, ReapeatMasker

=head1 COPYRIGHT

Copyright 2004 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/..";
use Getopt::Long;
use Data::Dumper;
use Taxonomy;

#
# Version
#
my $Version = 0.1;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',      # print out the version and exit
                    '-species=s',
                    '-isa=s',
                    '-taxDBFile=s'
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

usage() if ( !exists $options{'species'} );

my $taxFile = "$FindBin::Bin/../Libraries/taxonomy.dat";
if ( defined $options{'taxDBFile'} && -s $options{'taxDBFile'} ) {
  $taxFile = $options{'taxDBFile'};
}
die "queryTaxonmyDatabase: Could not open RepeatMasker Taxonomy "
    . "Database file $taxFile!\n"
    if ( !-s $taxFile );

print "\n\nRepeatMasker Taxonomy Database Utility\n";
print "======================================\n";
my $tax = Taxonomy->new( taxonomyDataFile => $taxFile );

if ( $tax->isSpecies( $options{'species'} ) ne "" ) {
  if ( $options{'species'} && $options{'isa'} ) {
    print "" . $options{'species'} . " is a " . $options{'isa'} . ": ";
    if ( $tax->isA( $options{'species'}, $options{'isa'} ) > 0 ) {
      print "True\n";
    }
    else {
      print "False\n";
    }
  }
  elsif ( $options{'species'} ) {
    print "Species = " . $options{'species'} . "\n";
    print "Lineage = ";
    foreach my $ancestor ( $tax->getLineage( $options{'species'} ) ) {
      print "$ancestor\n          ";
    }
    print "\n";
  }
}
else {
  print "Species "
      . $options{'species'}
      . " is not in the database. "
      . "Here is a list of possible similar\nsounding substitutes:\n";
  foreach
      my $species ( $tax->getSimilarSoundingSpecies( $options{'species'}, 10 ) )
  {
    print "  $species\n";
  }
}

1;
