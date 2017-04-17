#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) buildRMLibFromEMBL.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to read the RepeatMasker EMBL format
##      database and automatically create the old RepeatMasker.lib
##      FASTA style databse.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2005 Developed by
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
#     $Log: buildRMLibFromEMBL.pl,v $
#     Revision 1.47  2017/02/01 21:01:56  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME
                                                                                
buildRMLibFromEMBL.pl - Convert new EMBL format to old Fasta database.

=head1 SYNOPSIS
                                                                                
  buildRMLibFromEMBL.pl RepeatMaskerLib.EMBL > RepeatMasker.lib
                                                                                
=head1 DESCRIPTION
                                                                                
  A utility script to convert the new EMBL format repeat database
  back into the original FASTA format.
                                                                                
=head1 SEE ALSO
                                                                                
ReapeatMasker
                                                                                
=head1 COPYRIGHT
                                                                                
Copyright 2005 Robert Hubley, Institute for Systems Biology
                                                                                
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
use RepbaseEMBL;

#
# Version
#
my $Version = 0.1;

#
# Parameters
#
my $inFile;
if ( -f $ARGV[ 0 ] ) {
  $inFile = $ARGV[ 0 ];
}
else {
  die "Missing RepeatMasker EMBL database!";
}

##
## Convert the file
##
my $EMBLFile = "RepeatMaskerLib.embl";
print STDERR "Reading the RepeatMasker EMBL file: $inFile\n";
my $db = RepbaseEMBL->new( fileName => $inFile );

my $seqCount = $db->getRecordCount();

for ( my $i = 0 ; $i < $seqCount ; $i++ ) {
  my $record = $db->getRecord( $i );
  my $id     = $record->getId();
  my $type   = "#" . $record->getRMType();

  if ( $record->getRMSubType() ne "" ) {
    $type .= "/" . $record->getRMSubType();
  }
  my $desc = $record->getDescription();

  my $speciesList = "";
  foreach my $name ( $record->getRMSpeciesArray() ) {
    $speciesList .= "@" . $name . " ";
  }

  my $stageList = "[S:";
  my @stages    = $record->getRMSearchStagesArray();
  foreach my $stage ( @stages ) {
    $stageList .= "$stage,";
  }
  $stageList .= "]";
  $stageList =~ s/,\]/\]/;

  # Write the sequence
  my $seq = $record->getSequence();
  print ">" . $id . "$type $speciesList $stageList $desc\n";
  $seq =~ s/(\S{50})/$1\n/g;
  $seq .= "\n"
      unless ( $seq =~ /.*\n+$/s );
  print $seq;

  # Write the Buffered Sequence
  @stages = $record->getRMBufferStagesArray();
  my %stageHash = ();
  foreach my $stage ( @stages ) {
    if ( $stage =~ /(\d+)\[(\d+)\-(\d+)\]/ ) {
      push @{ $stageHash{"$2-$3"} }, $1;
    }
    elsif ( $stage =~ /(\d+)/ ) {
      push @{ $stageHash{"full"} }, $1;
    }
    else {
      print "RepeatMasker::createLib: Warning buffer stage $stage "
          . "understood!\n";
    }
  }
  foreach my $bufferSeqs ( keys %stageHash ) {
    $seq = $record->getSequence();
    $stageList = "[S:" . join( ",", @{ $stageHash{$bufferSeqs} } ) . "]";
    if ( $bufferSeqs eq "full" ) {
      $type = "#buffer";
    }
    elsif ( $bufferSeqs =~ /(\d+)-(\d+)/ ) {
      $seq = substr( $seq, $1 - 1, $2 - $1 + 1 );
      $type = "_$1" . "_$2#buffer";
    }
    print ">" . $id . "$type $speciesList $stageList $desc\n";
    $seq =~ s/(\S{50})/$1\n/g;
    $seq .= "\n"
        unless ( $seq =~ /.*\n+$/s );
    print $seq;
  }
}

1;
