#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) DFAM.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a simple datastructure for encapsulating the DFAM
##      formatted profile data and RepeatMasker meta-data.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
#
#******************************************************************************
#

=head1 NAME

DFAM


=head1 SYNOPSIS

use DFAM

Usage:

    $DB = DFAM->new( fileName => "dfam-hmmbuild-20111027.hmm" );

=head1 DESCRIPTION


=head1 INSTANCE METHODS

=cut

package DFAM;
use strict;
use Carp;
use Data::Dumper;
use DFAMRecord;
use Text::Wrap;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw();

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = "0.1";
my $CLASS   = "DFAM";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  $this->{'records'} = [];

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  if ( %nameValuePairs ) {
    my %validKeys = ( fileName => 1 );
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      unless ( exists $validKeys{$name} ) {
        croak "$CLASS::new: Error $name is not a valid "
            . "parameter to the constructor!\n\n"
            . "Usage: \$DB = $CLASS->new();\n"
            . "       \$DB = $CLASS->new( fileName=>\"dfam.hmm\" );\n";
      }
      $this->{$name} = $value;
    }
  }

  if ( exists $this->{'fileName'} ) {
    _parseFromFile( $this, $this->{'fileName'} );
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 Use: $DB->writeEMBLFile( $fileName );

  Write out this object as an EMBL file;

=cut

##-------------------------------------------------------------------------##
sub writeEMBLFile {
  my $this     = shift;
  my $fileName = shift;

  open OUT, ">$fileName"
      or die "$CLASS::writeEMBLFile() Unable to open "
      . "file $fileName for writing: $!\n";

  foreach my $record ( @{ $this->{'records'} } ) {
    print OUT "" . &toEMBLString( $record );
  }
  close OUT;

}

sub getRecordCount {
  my $this = shift;

  return ( $#{ $this->{'records'} } + 1 );
}

sub getRecord {
  my $this  = shift;
  my $index = shift;

  return if ( $index < 0 || $index > $#{ $this->{'records'} } );
  return $this->{'records'}->[ $index ];

}

sub getRelease {
  my $this = shift;

  return $this->{'release'};
}

##-------------------------------------------------------------------------##
## Private Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _parseFromFile( $obj, $filename );
##
##   Parse a FASTA output file into this datastructure.
##
##-------------------------------------------------------------------------##
sub _parseFromFile {
  my $this     = shift;
  my $filename = shift;

  open IN, "<$filename"
      || die "$CLASS::_parseFromFile() Unable to open "
      . "Repbase/EMBL file $filename: $!\n";

  my $inRepeatMaskerAnnotBlock = 0;
  my $data                     = "";
  my $recordRef                = DFAMRecord->new();
  my $recordLines              = "";
  my $id                       = "";
  my $acc                      = "";
  while ( <IN> ) {
    ## Get file release if it exists
    if ( /^#\s+Release:\s+(\S+)/ ) {
      $this->{'release'} = $1;
    }

    # IDs are just plain id's - no type/subtype.
    # Need to add back when we create our own libraries.
    if ( /^NAME\s+(\S+)/ ) {
      $id = $1;
    }

    # Accession
    if ( /^ACC\s+(\S+)/ ) {
      $acc = $1;
    }

    #
    if ( /^HMMER(\d+).*/ ) {
      if ( $1 != 3 ) {
        warn "DFAM.pm: Warning HMM format is not recognized: $_\n";
      }
      $recordLines = $_;
      $id          = "";
    }
    else {

      # Create substitution TAG
      if ( /^NAME\s+\S+/ ) {
        $recordLines .= "NAME #NEWID#\n";
      }
      else {
        $recordLines .= $_;
      }
    }

# New TH tag
#TH    TaxId:10090; TaxName:Mus musculus; GA:13.90; TC:33.14; NC:13.85; fdr:0.002;
    if (
/^TH\s+TaxId:\s*(\d+);\s*TaxName:\s*([^;]+);\s*GA:\s*([\d\.]+);\s*TC:\s*([\d\.]+);\s*NC:\s*([\d\.]+);\s*fdr:\s*([\d\.]+);\s*$/
        )
    {
      $recordRef->pushThresh(
                              {
                                'taxid'   => $1,
                                'taxname' => $2,
                                'hit_ga'  => $3,
                                'hit_tc'  => $4,
                                'hit_nc'  => $5,
                                'fdr'     => $6
                              }
      );
    }

    ## Process the end record tag
    if ( /^\/\/\s*$/ || eof( IN ) ) {

      # Create new ID
      my $newId = $id;
      if ( $recordRef->getRMType() ) {
        $newId .= "#" . $recordRef->getRMType();
      }
      if ( $recordRef->getRMSubType() ) {
        $newId .= "/" . $recordRef->getRMSubType();
      }
      $recordLines =~ s/NAME\s+#NEWID#/NAME $newId/g;

      # Save record!
      $recordRef->setRecordLines( $recordLines );
      $recordRef->setId( $id );
      $recordRef->setAcc( $acc );
      push @{ $this->{'records'} }, $recordRef;
      $recordRef = DFAMRecord->new();
      next;
    }

    ## Capture the tag!
    if ( /^CC(?:\s+(\S.*))?/ ) {

      # New tag so overwrite old data
      if ( defined $1 ) {
        $data .= $1 . "\n";
      }
      else {
        $data .= "\n";
      }
      $inRepeatMaskerAnnotBlock = 1;
    }
    else {
      $inRepeatMaskerAnnotBlock = 0;
      if ( $data ne "" ) {

        # Parse last tag
        _parseTagData( $data, $recordRef );
      }
      $data = "";
    }
  }

  close IN;
}

##-------------------------------------------------------------------------##
## Use: my _parseTagData( $tag, $data, $recRef );
##
##  Parse a RepeatMasker annotation block.  This is a helper function
##  for the general parser.  This function expects data for
##  multiline tags to already be compressed into a single
##  blob.  Ie.
##
##     ##  line1
##     ##  line2
##
##  Should be passed to this function as:
##
##     $data = "line1\nline2"
##
##-------------------------------------------------------------------------##
sub _parseTagData {
  my $data   = shift;
  my $recRef = shift;

  if ( $data =~ /RepeatMasker\s+Annotations/i ) {
    if ( $data =~ /Type:[ \t]*(\S+)[ \t]*[\n\r]/i ) {
      $recRef->setRMType( $1 );
    }
    if ( $data =~ /SubType:[ \t]*(\S+)[ \t]*[\n\r]/i ) {
      $recRef->setRMSubType( $1 );
    }

    if ( $data =~ /Species:[ \t]*([\S \t]+)[ \t]*[\n\r]/i ) {
      my $speciesStr = $1;

      $speciesStr =~ s/\s//g;
      $speciesStr =~ s/(.*),$/$1/g;
      my @species = split( /\,/, $speciesStr );
      foreach my $species ( @species ) {
        if ( $species =~ /Search/ ) { die "Hmmm $speciesStr\n"; }
        $species =~ s/\s//g;
        $recRef->pushRMSpecies( $species );
      }
    }
    if ( $data =~ /Refineable/ ) {
      $recRef->setRMRefineable( 1 );
    }
    if ( $data =~ /SearchStages:\s*([\d\s*\,\s*]+)[\n\r]/i ) {
      my $stageStr = $1;
      $stageStr =~ s/\s//g;
      $stageStr =~ s/(.*),$/$1/g;
      my @stages = split( /\,/, $stageStr );
      foreach my $stage ( @stages ) {
        $recRef->pushRMSearchStages( $stage );
      }
    }
    if ( $data =~ /BufferStages:([^\n\r]*)[\n\r]+/i ) {
      my $bufData = $1;
      while ( $bufData =~ s/\s*(\d+\s*(?:\[\s*\d+\s*-\s*\d+\s*\]\s*)*),?//m ) {
        my $bufStr = $1;
        $bufStr =~ s/\s//g;
        if ( $bufStr ne "" ) {
          $recRef->pushRMBufferStages( $bufStr );
        }
      }
    }
    if ( $data =~ /Description:\s+(\S.*)/i ) {
      $recRef->setRMDescription( $1 );
    }
  }

  # Coment line
  # Multi-record/multi-line comment allowed
  $recRef->pushComments( $data );

}

1;
