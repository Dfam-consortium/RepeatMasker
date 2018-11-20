#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) TRFResult.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A class for holding a single result returned by the
##      trf ( Tandem Repeat Finder ) program.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2005 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless( {
#          'start' => '11376',
#          'end' => '11412',
#          'period' => '16',
#          'copyNumber' => '2.3',
#          'consSize' => '16',
#          'percMatches' => '95',
#          'percIndels' => '0',
#          'score' => '65',
#          'percA' => '40',
#          'percC' => '8',
#          'percG' => '18',
#          'percT' => '32',
#          'entropy' => '1.8',
#          'consensus' => 'GTAAAGATTTGCACAT'
#        }
#     }, 'TRFResult' );
#
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#
#    Consider adding storage for the alignments.
#

=head1 NAME

TRFResult

=head1 SYNOPSIS

use TRFResult

Usage: 

    $TRFResult = TRFResult->new();

  or 

    $TRFResult = TRFResult->new( 
                  start => '11376', end => '11412',
                  period => '16', copyNumber => '2.3',
                  consSize => '16', percMatches => '95',
                  percIndels => '0', score => '65',
                  percA => '40', percC => '8',
                  percG => '18', percT => '32',
                  entropy => '1.8', 
                  consensus => 'GTAAAGATTTGCACAT'
                                           );

=head1 DESCRIPTION

A class for storing a result from the trf search engine.

=head1 SEE ALSO

TRF

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package TRFResult;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "TRFResult";

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( "TRFResult::add: Instance variable $name doesn't exist." . "" );
      }
      $this->$method( $value );
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a TRFResult *duplicating* all the values of the old
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

##---------------------------------------------------------------------##

=head2 get_setSeqID()

  Use: my $value    = getSeqID( );
  Use: my $oldValue = setSeqID( $value );

  Get/Set the SeqID attribute.  

=cut

##---------------------------------------------------------------------##
sub getSeqID {
  my $this = shift;

  return $this->{'seqid'};
}

sub setSeqID {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'seqid'};
  $this->{'seqid'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setStart()

  Use: my $value    = getStart( );
  Use: my $oldValue = setStart( $value );

  Get/Set the Start attribute.  

=cut

##---------------------------------------------------------------------##
sub getStart {
  my $this = shift;

  return $this->{'start'};
}

sub setStart {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'start'};
  $this->{'start'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setEnd()

  Use: my $value    = getEnd( );
  Use: my $oldValue = setEnd( $value );

  Get/Set the End attribute.  

=cut

##---------------------------------------------------------------------##
sub getEnd {
  my $this = shift;

  return $this->{'end'};
}

sub setEnd {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'end'};
  $this->{'end'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPeriod()

  Use: my $value    = getPeriod( );
  Use: my $oldValue = setPeriod( $value );

  Get/Set the Period attribute.  

=cut

##---------------------------------------------------------------------##
sub getPeriod {
  my $this = shift;

  return $this->{'period'};
}

sub setPeriod {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'period'};
  $this->{'period'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setCopyNumber()

  Use: my $value    = getCopyNumber( );
  Use: my $oldValue = setCopyNumber( $value );

  Get/Set the CopyNumber attribute.  

=cut

##---------------------------------------------------------------------##
sub getCopyNumber {
  my $this = shift;

  return $this->{'copyNumber'};
}

sub setCopyNumber {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'copyNumber'};
  $this->{'copyNumber'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setConsSize()

  Use: my $value    = getConsSize( );
  Use: my $oldValue = setConsSize( $value );

  Get/Set the ConsSize attribute.  

=cut

##---------------------------------------------------------------------##
sub getConsSize {
  my $this = shift;

  return $this->{'consSize'};
}

sub setConsSize {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'consSize'};
  $this->{'consSize'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercMatches()

  Use: my $value    = getPercMatches( );
  Use: my $oldValue = setPercMatches( $value );

  Get/Set the PercMatches attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercMatches {
  my $this = shift;

  return $this->{'percMatches'};
}

sub setPercMatches {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percMatches'};
  $this->{'percMatches'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercIndels()

  Use: my $value    = getPercIndels( );
  Use: my $oldValue = setPercIndels( $value );

  Get/Set the PercIndels attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercIndels {
  my $this = shift;

  return $this->{'percIndels'};
}

sub setPercIndels {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percIndels'};
  $this->{'percIndels'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setScore()

  Use: my $value    = getScore( );
  Use: my $oldValue = setScore( $value );

  Get/Set the Score attribute.  

=cut

##---------------------------------------------------------------------##
sub getScore {
  my $this = shift;

  return $this->{'score'};
}

sub setScore {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'score'};
  $this->{'score'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercA()

  Use: my $value    = getPercA( );
  Use: my $oldValue = setPercA( $value );

  Get/Set the PercA attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercA {
  my $this = shift;

  return $this->{'percA'};
}

sub setPercA {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percA'};
  $this->{'percA'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercC()

  Use: my $value    = getPercC( );
  Use: my $oldValue = setPercC( $value );

  Get/Set the PercC attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercC {
  my $this = shift;

  return $this->{'percC'};
}

sub setPercC {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percC'};
  $this->{'percC'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercG()

  Use: my $value    = getPercG( );
  Use: my $oldValue = setPercG( $value );

  Get/Set the PercG attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercG {
  my $this = shift;

  return $this->{'percG'};
}

sub setPercG {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percG'};
  $this->{'percG'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPercT()

  Use: my $value    = getPercT( );
  Use: my $oldValue = setPercT( $value );

  Get/Set the PercT attribute.  

=cut

##---------------------------------------------------------------------##
sub getPercT {
  my $this = shift;

  return $this->{'percT'};
}

sub setPercT {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'percT'};
  $this->{'percT'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setEntropy()

  Use: my $value    = getEntropy( );
  Use: my $oldValue = setEntropy( $value );

  Get/Set the Entropy attribute.  

=cut

##---------------------------------------------------------------------##
sub getEntropy {
  my $this = shift;

  return $this->{'entropy'};
}

sub setEntropy {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'entropy'};
  $this->{'entropy'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setConsensus()

  Use: my $value    = getConsensus( );
  Use: my $oldValue = setConsensus( $value );

  Get/Set the Consensus attribute.  

=cut

##---------------------------------------------------------------------##
sub getConsensus {
  my $this = shift;

  return $this->{'consensus'};
}

sub setConsensus {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'consensus'};
  $this->{'consensus'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setSubjSeq()

  Use: my $value    = getSubjSeq( );
  Use: my $oldValue = setSubjSeq( $value );

  Get/Set the subject sequence attribute.  

=cut

##---------------------------------------------------------------------##
sub getSubjSeq {
  my $this = shift;

  return $this->{'subjseq'};
}

sub setSubjSeq {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'subjseq'};
  $this->{'subjseq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $str  = ""
      . $this->getSeqID() . " - "
      . $this->getStart() . " - "
      . $this->getEnd() . " "
      . $this->getPeriod() . " "
      . $this->getCopyNumber() . " "
      . $this->getConsSize() . " "
      . $this->getPercMatches() . " "
      . $this->getPercIndels() . " "
      . $this->getScore() . " "
      . $this->getPercA() . " "
      . $this->getPercC() . " "
      . $this->getPercG() . " "
      . $this->getPercT() . " "
      . $this->getEntropy() . " "
      . $this->getConsensus() . " "
      . $this->getSubjSeq() . "\n";

  return $str;

}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;
