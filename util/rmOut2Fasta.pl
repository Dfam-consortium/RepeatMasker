#!/usr/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) rmOut2Fasta.pl
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      Using RepeatMasker's batching mechanism generate a particular
##      numbered batch from the input file.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2002-2011 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
#  ChangeLog:
#
#    $Log$
#
###############################################################################
#
# To Do:
#
#

=head1 NAME

rmOut2Fasta.pl - Extract FASTA sequences based on RM *.out annotations

=head1 SYNOPSIS

  rmOut2Fasta.pl [-options] -fasta <*.fa> -out <*.out>

=head1 DESCRIPTION

The options are:

=over 4

=item -fasta <*.fa>

A single FASTA file containing all the sequences referenced in
a RepeatMasker *.out file.

=item -out <*.out>

A single RepeatMasker *.out file.

=back

=head1 SEE ALSO

=over 4

RepeatMasker

=back

=head1 COPYRIGHT

Copyright 2007-2013 Robert Hubley, Arian Smit, Institute for Systems Biology

=head1 AUTHORS

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use Getopt::Long;
use POSIX qw(:sys_wait_h);
use File::Copy;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd;

# RepeatMasker Libraries
use FastaDB;
use CrossmatchSearchEngine;
use SearchResultCollection;

# Debugging flag
my $DEBUG = 0;

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
#  check out the code as "cvs co RepeatMasker" the
#  $Name:  $ macro will be blank and thus we use
#  this as the development version.
#
my $CVSTag = '$Name:  $';
my $version;
if ( $CVSTag =~ /\$\s*Name:\s*open-(\S+)\s*\$/ ) {
  $version = $1;
  $version =~ s/-/./g;
  $version = "open-$version";
}
else {
  $version =
      'development-$Id: rmOut2Fasta.pl,v 1.6 2017/02/01 21:01:57 rhubley Exp $';
}

my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-fasta=s',
                    '-out=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

# Print the internal POD documentation if something is missing
if ( !( $options{'fasta'} && $options{'out'} ) ) {
  print "Missing -fasta or -out options!\n\n";
  exec "pod2text $0";
  die;
}

my $db = FastaDB->new(
                       fileName    => $options{'fasta'},
                       openMode    => SeqDBI::ReadOnly,
                       maxIDLength => 50
);

my $resultCollection =
    CrossmatchSearchEngine::parseOutput( searchOutput => $options{'out'} );

for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  my $result = $resultCollection->get( $i );
  my $qID    = $result->getQueryName();
  my $qBeg   = $result->getQueryStart();
  my $qEnd   = $result->getQueryEnd();
  my $sID    = $result->getSubjName();
  my $sBeg   = $result->getSubjStart();
  my $sEnd   = $result->getSubjEnd();
  my $orient = $result->getOrientation();
  $orient = "+" if ( $orient ne "C" );
  $orient = "-" if ( $orient eq "C" );

  my $seq = $db->getSubstr( $qID, $qBeg - 1, $qEnd - $qBeg + 1 );

  if ( $orient eq "-" ) {
    $seq = reverse( $seq );
    $seq =~ tr/ACGTYRMKHBVDacgtyrmkhbvd/TGCARYKMDVBHtgcarykmdvbh/;
  }
  $seq =~ s/(\S{50})/$1\n/g;
  $seq .= "\n"
      unless ( $seq =~ /.*\n+$/s );

  print ">$qID:$qBeg-$qEnd ( $sID:$sBeg-$sEnd orient=$orient )\n";
  print "$seq";
}

