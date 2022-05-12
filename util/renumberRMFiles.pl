#!/usr/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) renumberRMFiles
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Fix ID field in *.out ( and possibly linked *.align *.cat files )
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

renumberRMFiles - Fix ID field in *.out ( and possibly linked *.align, *.cat files )

=head1 SYNOPSIS

  renumberRMFiles [-version] -out <*.out[.gz]>

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=item -out

The RepeatMasker *.out file to change.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2021 Robert Hubley, Institute for Systems Biology

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
use EMBL;

my $Version  = "1.0";

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
                    '-version',                 # print out the version and exit
                    '-out=s',
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

usage() if ( ! -s $options{'out'} );

#
# Collect data from out file
#
my %newIDs = ();
my $lastID        = -1;
my $lastSeqStart  = 0;
my $lastSeqName   = "";
my $nextID     = 1;

my $file = $options{'out'};
if ( $file =~ /.*\.gz/ ) {
  open IN, "gunzip -c $file|"
      or die "Could not open $file using gunzip!\n";
}
else {
  open IN, "<$file"
      or die "Could not open $file!\n";
}
while ( <IN> ) {
  if ( /^\s*(\d+\s+\d+\.\d+.*)/ ) {
    my $data = $1;
    my @fields  = split(/\s+/,$data);
    my $seqName = $fields[ 4 ];
    my $seqStart= $fields[ 5 ];
    my $id      = $fields[ 14 ];
    my $origID  = $id;
    # Example:
    #   1446 21.7  4.9  1.3 chr20             39784121 39784437 (24659730) + MLT1D           LTR/ERVL-MaLR         1     329   (176) 77945
    #    677 14.4  5.7  0.8 chr20             39784457 39784589 (24659578) + MLT1D           LTR/ERVL-MaLR       366     505     (0) 9655
    if ( $seqName ne $lastSeqName
         || $seqStart < $lastSeqStart
         || ($id < $lastID - 50 && $id < 3)
         || (abs($id-$lastID) > 200) ) {
      %newIDs = ();
    }
    $lastID = $id;
    $lastSeqStart = $seqStart;
    $lastSeqName = $seqName;

    if ( ! defined $id || $id eq "" ) 
    {
      $id = $nextID;
      $newIDs{$id} = $nextID;
      $nextID++;
    }elsif ( exists $newIDs{$id} ) {
      # Already renumbered....keep the same id
      $id = $newIDs{$id};
    }else {
      $newIDs{$id} = $nextID;
      $id = $nextID;
      $nextID++;
    }
    $fields[14] = $id;
    print "" . join(" ",@fields) . "\n";

  }
}
close IN;

