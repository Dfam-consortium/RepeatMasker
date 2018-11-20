#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) TRFTRFSearchResult.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An object for holding a TRF search result.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2012 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#
#

=head1 NAME

TRFSearchResult

=head1 SYNOPSIS

use TRFSearchResult

Usage: 

    $TRFSearchResultsCollection = TRFSearchResult->new();

  or 

    $TRFSearchResultsCollection = TRFSearchResult->new( 
                                 queryName=>value, subjName=>value,
                                 pctInsert=>value, pctDelete=>value, 
                                 queryStart=>value, queryEnd=>value, 
                                 score=>value, pctDiverge=>value, 
                                 subjRemaining=>value, subjType=>value,
                                 queryRemaining=>value, id=>value,
                                 orientation=>value, queryString=>value,
                                 subjString=>value, matrix=>value,
                                 id=>value, lineageId=>value );

=head1 DESCRIPTION

A class for storing a result from the TRF search engine.

=head1 SEE ALSO

=over 4

TRFSearchResultCollection

=back

=head1 COPYRIGHT

Copyright 2012 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package TRFSearchResult;
use strict;
use SearchResult;
use Data::Dumper;
use Matrix;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(SearchResult Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "TRFSearchResult";
my $DEBUG = 0;

use constant NoAlign           => "trf-1";
use constant AlignWithQuerySeq => "trf-2";

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = $class->SUPER::new( @_ );

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  ## Allow import of values
#if ( %nameValuePairs ) {
#  while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
#    my $method = "set" . _ucFirst( $name );
#    unless ( $this->can( $method ) ) {
#      croak(
#           "TRFSearchResult::add: Instance variable $name doesn't exist." . "" );
#    }
#    $this->$method( $value );
#  }
#}

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a TRFSearchResult *duplicating* all the values of the old
  object in the new one.

=cut

##-------------------------------------------------------------------------##
sub clone {
  my $this = shift;

  my %newHash = %{$this};
  my $newObj  = \%newHash;

  bless $newObj, ref( $this );

  return $newObj;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setConsensus()

  Use: my $value    = getConsensus( );
  Use: my $oldValue = setConsensus( $value );

  Get/Set the TRF consensus sequence.

=cut

##-------------------------------------------------------------------------##
sub getConsensus {
  my $obj = shift;

  my $value = $obj->{'consensus'};

  return $value;
}

sub setConsensus {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'consensus'};
  $obj->{'consensus'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPeriod()

  Use: my $value    = getPeriod( );
  Use: my $oldValue = setPeriod( $value );

  Get/Set the TRF repeat period.

=cut

##-------------------------------------------------------------------------##
sub getPeriod {
  my $obj = shift;

  my $value = $obj->{'period'};

  return $value;
}

sub setPeriod {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'period'};
  $obj->{'period'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setCopyNumber()

  Use: my $value    = getCopyNumber();
  Use: my $oldValue = setCopyNumber( $value );

  Get/Set the TRF copy number for the repeat instance.

=cut

##-------------------------------------------------------------------------##
sub getCopyNumber {
  my $obj = shift;

  my $value = $obj->{'copyNumber'};

  return $value;
}

sub setCopyNumber {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'copyNumber'};
  $obj->{'copyNumber'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setConsensusSize()

  Use: my $value    = getConsensusSize();
  Use: my $oldValue = setConsensusSize( $value );

  Get/Set the length of the TRF consensus.

=cut

##-------------------------------------------------------------------------##
sub getConsensusSize {
  my $obj = shift;

  my $value = $obj->{'consensusSize'};

  return $value;
}

sub setConsensusSize {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'consensusSize'};
  $obj->{'consensusSize'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $cmString = $obj->_toLocalFormat();
##
##   TODO: Document
##
##-------------------------------------------------------------------------##
sub _toLocalFormat {
  my $obj           = shift;
  my $alignmentMode = shift;

  $alignmentMode = TRFSearchResult::NoAlign
      if ( !defined $alignmentMode );
  croak $CLASS
      . "::toStringFormatted: Unknown alignment mode "
      . "( $alignmentMode )\n"
      if (    $alignmentMode != TRFSearchResult::NoAlign
           && $alignmentMode != TRFSearchResult::AlignWithQuerySeq );

  my $retStr = "";
  my $sbjDir;
  my $sbjDirStr;
  my $alignIndex = 0;
  my $alignCol;

  my ( $qryName ) = ( $obj->{qryName} =~ /(\S+).*/ );
  $retStr .=
        "$obj->{score} $obj->{percDiv} $obj->{percDel} "
      . "$obj->{percIns} $qryName $obj->{qryBegin} "
      . "$obj->{qryEnd} ($obj->{qryLeft}) ";
  my ( $sbjName ) = ( $obj->{sbjName} =~ /(\S+).*/ );
  if ( $obj->{sbjOrient} eq "C" ) {
    $retStr .=
        "C $sbjName ($obj->{sbjLeft}) " . "$obj->{sbjEnd} $obj->{sbjBegin}";
  }
  else {
    $retStr .=
        "$sbjName $obj->{sbjBegin} $obj->{sbjEnd} " . "($obj->{sbjLeft})";
  }
  if ( defined $obj->{id} ) {
    $retStr .= " $obj->{id}";
  }
  if ( defined $obj->{lineageId} ) {
    $retStr .= " $obj->{lineageId}";
  }
  if ( defined $obj->{overlap} ) {
    $retStr .= " $obj->{overlap}";
  }
  $retStr .= "\n";

  if (    $alignmentMode != TRFSearchResult::NoAlign
       && defined $obj->{querySeq}
       && $obj->{querySeq} ne "" )
  {
    my $qMasked;
    my $sMasked;
    my $insertions = 0;
    my $deletions  = 0;

    my $query   = $obj->{querySeq};
    my $subject = $obj->{subjSeq};

    $retStr .= "\n";

    #    if (    $obj->{sbjOrient} eq "C"
    #         && $alignmentMode == TRFSearchResult::AlignWithSubjSeq )
    #    {
    #      $query = reverse $query;
    #      $query =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
    #      $subject = reverse $subject;
    #      $subject =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
    #    }

    my $qStart = $obj->{qryBegin};

    #    if (    $obj->{sbjOrient} eq "C"
    #         && $alignmentMode == TRFSearchResult::AlignWithSubjSeq )
    #    {
    #      $qStart = $obj->{qryEnd};
    #    }
    my $qEnd   = 0;
    my $sStart = $obj->{sbjBegin};
    if (    $obj->{sbjOrient} eq "C"
         && $alignmentMode == TRFSearchResult::AlignWithQuerySeq )
    {
      $sStart = $obj->{sbjEnd};
    }
    my $sEnd = 0;
    while ( $query ) {
      $query =~ s/^(.{1,50})//;
      my $qSeq = $1;
      $subject =~ s/^(.{1,50})//;
      my $sSeq = $1;

      $insertions = ( $qSeq =~ tr/-/-/ );
      $deletions  = ( $sSeq =~ tr/-/-/ );

      if ( $sEnd > 0 ) {
        my $qIncr = 0;
        my $sIncr = 0;
        $qIncr = 1 if ( length( $qSeq ) > $insertions );
        $sIncr = 1 if ( length( $sSeq ) > $deletions );
        if ( $obj->{sbjOrient} eq "C" ) {

         #          if ( $alignmentMode == TRFSearchResult::AlignWithSubjSeq ) {
         #            $qStart = $qEnd - $qIncr;
         #            $sStart = $sEnd + $sIncr;
         #          }
         #          else {
          $qStart = $qEnd + $qIncr;
          $sStart = $sEnd - $sIncr;

          #          }
        }
        else {
          $qStart = $qEnd + $qIncr;
          $sStart = $sEnd + $sIncr;
        }
      }

      # Indicate orientation in alignment data by placing a "C" in front
      # of sequences which have been reverse complemented
      #      if (    $obj->{'sbjOrient'} eq "C"
      #           && $alignmentMode == TRFSearchResult::AlignWithSubjSeq )
      #      {
      #        $qEnd = $qStart - length( $qSeq ) + 1 + $insertions;
      #        $retStr .= "C ";
      #      }
      #      else {
      $qEnd = $qStart + length( $qSeq ) - 1 - $insertions;
      $retStr .= "  ";

      #      }
      $qEnd = $qStart if ( length( $qSeq ) == $insertions );

      # Up to 13 characters are allowed from the QueryName/SubjName
      $retStr .= substr( $obj->{'qryName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'qryName'} ) < 13 ? length( $obj->{'qryName'} ) : 13
           )
          );

      # Up to 9 characters for the position, followed by the sequence and
      # end positions
      $retStr .= " " x ( 9 - length( $qStart ) ) . $qStart . " $qSeq $qEnd\n";

      # The modification codes
      $retStr .= " " x 25;
      for ( my $j = 0 ; $j < length( $qSeq ) ; $j++ ) {
        my $qChar = substr( $qSeq, $j, 1 );
        my $sChar = substr( $sSeq, $j, 1 );
        if ( $qChar eq $sChar ) {
          $retStr .= " ";
        }
        elsif ( $qChar eq "-" || $sChar eq "-" ) {
          $retStr .= "-";
        }
        elsif (    $qChar . $sChar eq "CT"
                || $qChar . $sChar eq "TC"
                || $qChar . $sChar eq "AG"
                || $qChar . $sChar eq "GA" )
        {
          $retStr .= "i";
        }
        elsif (    $qChar . $sChar eq "GT"
                || $qChar . $sChar eq "TG"
                || $qChar . $sChar eq "GC"
                || $qChar . $sChar eq "CG"
                || $qChar . $sChar eq "CA"
                || $qChar . $sChar eq "AC"
                || $qChar . $sChar eq "AT"
                || $qChar . $sChar eq "TA" )
        {
          $retStr .= "v";
        }
        elsif (    ( $qChar =~ /[BDHVRYKMSWNX]/ )
                || ( $sChar =~ /[BDHVRYKMSWNX]/ ) )
        {
          $retStr .= "?";
        }
        else {
          $retStr .= " ";
        }
      }
      $retStr .= "\n";

      # Subject orientation/name
      if (    $obj->{'sbjOrient'} eq "C"
           && $alignmentMode == TRFSearchResult::AlignWithQuerySeq )
      {
        $retStr .= "C ";
        $sEnd = $sStart - length( $sSeq ) + 1 + $deletions;
      }
      else {
        $retStr .= "  ";
        $sEnd = $sStart + length( $sSeq ) - 1 - $deletions;
      }
      $retStr .= substr( $obj->{'sbjName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'sbjName'} ) < 13 ? length( $obj->{'sbjName'} ) : 13
           )
          );
      $retStr .= " " x ( 9 - length( $sStart ) ) . $sStart . " $sSeq $sEnd\n";
      $retStr .= "\n";

    }    # while ( $query )

    if ( defined $obj->getMatrixName() ) {
      $retStr .= "Matrix = " . $obj->getMatrixName() . "\n";
    }
    else {
      $retStr .= "Matrix = Unknown\n";
    }
    $retStr .= "Transitions / transversions = ";
    my ( $mismatches, $transitions, $transversions, $numGaps, $totGapLen ) =
        $obj->_getAlignmentStats();
    if ( defined $transitions ) {
      if ( $transversions > 0 ) {
        $retStr .= sprintf( "%0.2f", ( $transitions / $transversions ) );
      }
      else {
        $retStr .= "0.0";
      }
      $retStr .= " ($transitions / $transversions)\n";
      if ( ( $obj->getQueryEnd() - $obj->getQueryStart() ) > 0 ) {
        $retStr .= "Gap_init rate = "
            . sprintf( "%0.2f",
                       $numGaps /
                           ( $obj->getQueryEnd() - $obj->getQueryStart() ) )
            . " ($numGaps / "
            . ( $obj->getQueryEnd() - $obj->getQueryStart() ) . ")";
      }
      else {
        $retStr .= "Gap_init rate = 0.0 ( $numGaps / 0 )";
      }
      if ( $numGaps > 0 ) {
        $retStr .=
              ", avg. gap size = "
            . sprintf( "%0.2f", ( $totGapLen / $numGaps ) )
            . " ($totGapLen / $numGaps)\n\n";
      }
      else {
        $retStr .= ", avg. gap size = 0.0 (0 / 0)\n\n";
      }

    }
    else {
      $retStr .= "Transitions / transversions = Unknown\n";
      $retStr .= "Gap_init rate = Unknown, avg. gap size = Unknown\n\n";
    }
  }

  return $retStr;
}

1;
