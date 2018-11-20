#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) FastaDB.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      An implementation of the SeqDBI interface which uses FASTA
##      formatted flat-files as the backing store.
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
#  Object DataStructure
#  ====================
#
#  $this->{'fastaInfo'} =
#        { 'filename' => #The name of the fasta file#,
#          'openMode' = SeqDBI::ReadOnly | SeqDBI::ReadWrite
#          'sequences' => [ { 'id' =>   #The unique & cleaned up seq name#,
#                             'name' => #The untouched seq name#,
#                             'description' => #The original description#,
#                             'length' => #The seq length#,
#                             'subtLength' => #The seq lengh minus XRNYMK#,
#                             'gcLength' => #The number of G&C's in seq#,
#                             'xnLength' => #The seq length minus XN#,
#                             'startLine' => #Fasta line number for seq#,
#                             'startByte' => #Fasta start byte for seq#,
#                             'seqPosIndices' => [ $seqPos, $fByte, $fLine
#                                                  $seqPos2, $fByte2, $fLine2
#                                                  #These are tripples
#                                                  # which can be used
#                                                  # as an index into
#                                                  # the fasta file. The
#                                                  # seqPos are always
#                                                  # positions on the
#                                                  # start of a sequence
#                                                  # line.
#                                                      .. ]
#                           },
#                           { 'id' =>
#                              ..
#                           }
#                         ]
#            'seqIDHash' => { #id# => #index into 'sequences' record#
#                             #id# ...
#                           }
#        };
#
################################################################################
# ChangeLog:
#
#   $Log$
#
################################################################################
#
## To Do:
#
#
#

=head1 NAME

FastaDB - An implementation of SeqDBI.pm interface which uses a FASTA formatted
file as the backing store.

=head1 SEE ALSO

=over 4

RepeatMasker

=back

=head1 COPYRIGHT

Copyright 2004 Robert Hubley Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 METHODS

=cut

#
# Module Dependence and Package Definition
#
package FastaDB;
use strict;
use Data::Dumper;
use Tie::File;
use SeqDBI;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SeqDBI);

@EXPORT = qw();

@EXPORT_OK = qw(min max);

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Constants
#
use constant IndexDistanceDefault => 10000;
use constant FastaSeqWidthDefault => 50;

#
# Version
#
my $VERSION = "";

#
# Globals
#
my $DEBUG        = 0;
my $CLASS        = "FastaDB";
my $seqIndexDist = IndexDistanceDefault;
my $fastaLineLen = FastaSeqWidthDefault;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 new()

  Use: FastaDB->new( [fileName=>"filename.fa"],
                         [openMode=>SeqDBI::ReadOnly|SeqDBI::ReadWrite],
                         [indexDist=>1000] );

  Create a FastaDB object.

=cut

##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->{'fastaInfo'}->{'openMode'} = SeqDBI::ReadWrite;
  $this->{'fastaInfo'}->{'openMode'} = $nameValuePairs{'openMode'}
      if ( defined $nameValuePairs{'openMode'} );

  $seqIndexDist = $nameValuePairs{'indexDist'}
      if ( defined $nameValuePairs{'indexDist'}
           && $nameValuePairs{'indexDist'} =~ /\d+/ );

  if ( defined $nameValuePairs{'maxIDLength'} ) {
    if ( $nameValuePairs{'maxIDLength'} =~ /^\s*(\d+)\s*$/ ) {
      $this->{'fastaInfo'}->{'maxIDLength'} = $1;
    }
    else {
      croak "" . $class
          . "::new(): Attribute maxIDLength is not a number!: "
          . $nameValuePairs{'maxIDLength'} . "\n";
    }
  }

  # Process fasta file if passed to us.
  if ( defined $nameValuePairs{'fileName'} ) {
    $this->{'fastaInfo'}->{'filename'} = $nameValuePairs{'fileName'};
    $this->compact();
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Public Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 getSeqCount()

  Use: my $seqCount =  &getSeqCount( $obj );

  The number of sequences contained in this sequence database.

=cut

##-------------------------------------------------------------------------##
sub getSeqCount {
  my $this = shift;

  return ( 0 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( $#{ $this->{'fastaInfo'}->{'sequences'} } + 1 );

}

##-------------------------------------------------------------------------##

=head2 exists()

  Use: my $bool = &exists( $obj, $seqID );

  Checks to see if the sequence exists in the database.

=cut

##-------------------------------------------------------------------------##
sub exists {
  my $this = shift;
  my $id   = shift;

  return ( 0 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  if ( exists $this->{'fastaInfo'}->{'seqIDHash'}->{$id} ) {
    return ( 1 );
  }

  return ( 0 );

}

##-------------------------------------------------------------------------##

=head2 getSeqLength()

  Use: my $length = &getSeqLength( $obj, [$seqID] );

  The sequence length of the entire database or of a
  particular entry.

=cut

##-------------------------------------------------------------------------##
sub getSeqLength {
  my $this = shift;

  return ( -1 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  my $length = 0;

  if ( @_ ) {
    my $id = shift();
    if ( defined $this->{'fastaInfo'}->{'seqIDHash'}->{$id} ) {
      my $seqIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{$id};
      $length = $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'length'};
    }
  }
  else {
    foreach my $seq ( @{ $this->{'fastaInfo'}->{'sequences'} } ) {
      $length += $seq->{'length'};
    }
  }

  return $length;

}

##-------------------------------------------------------------------------##

=head2 getSeqLengths()

  Use: my @sequenceLength = &getSeqLengths( $obj );

  Get an array of sequence length values for all entries
  in the database.

=cut

##-------------------------------------------------------------------------##
sub getSeqLengths {
  my $this = shift;

  my @seqLengths = ();
  return ( @seqLengths ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( map { $_->{'length'} } @{ $this->{'fastaInfo'}->{'sequences'} } );

}

##-------------------------------------------------------------------------##

=head2 getXNLength()

  Use: my $xnLen = &getXNLength( $obj, [$seqID] );

  Get the length of $seqID sequence minus the masking and
  scaffold characters 'X' and 'N' respectively.  If $seqID is
  omited the value returned is the length of database
  minus all "X" and "N" characters.

=cut

##-------------------------------------------------------------------------##
sub getXNLength {
  my $this = shift;

  return ( -1 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  my $length = 0;

  if ( @_ ) {
    my $seqIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{ shift() };
    $length = $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'xnLength'};
  }
  else {
    foreach my $seq ( @{ $this->{'fastaInfo'}->{'sequences'} } ) {
      $length += $seq->{'xnLength'};
    }
  }

  return $length;

}

##-------------------------------------------------------------------------##

=head2 getXNLengths()

  Use: my @xnLengths = &getXNLengths( $obj );

  Get the length of each sequence minus 'X' and 'N' characters
  ( case insensitive ).  The values are returned in an array.

=cut

##-------------------------------------------------------------------------##
sub getXNLengths {
  my $this = shift;

  my @seqLengths = ();
  return ( @seqLengths ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( map { $_->{'xnLength'} } @{ $this->{'fastaInfo'}->{'sequences'} } );

}

##-------------------------------------------------------------------------##

=head2 getGCLength()

  Use: my $gcLen = &getGCLength( $obj, [$seqID] );

  Get the length of $seqID sequence minus all characters
  except "G" and "C" ( case insensitive ).  If $seqID is
  omited the value returned is the length of database
  minus all characters except "G" and "C".

=cut

##-------------------------------------------------------------------------##
sub getGCLength {
  my $this = shift;

  return ( -1 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  my $length = 0;

  if ( @_ ) {
    my $seqIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{ shift() };
    $length = $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'gcLength'};
  }
  else {
    foreach my $seq ( @{ $this->{'fastaInfo'}->{'sequences'} } ) {
      $length += $seq->{'gcLength'};
    }
  }

  return $length;

}

##-------------------------------------------------------------------------##

=head2 getGCLengths()

  Use: my @gcLengths = &getGCLengths( $obj );

  Get the length of each sequence minus all characters
  except "G" and "C" ( case insensitive ).  The values
  are returned in an array.

=cut

##-------------------------------------------------------------------------##
sub getGCLengths {
  my $this = shift;

  my @seqLengths = ();
  return ( @seqLengths ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( map { $_->{'gcLength'} } @{ $this->{'fastaInfo'}->{'sequences'} } );

}

##-------------------------------------------------------------------------##

=head2 getSubtLength()

  Use my $subtLen = &getSubtLength( $obj, [$seqID] );

  Get length of sequence minus the ambiguous characters.

=cut

##-------------------------------------------------------------------------##
sub getSubtLength {
  my $this = shift;

  return ( -1 ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  my $length = 0;

  if ( @_ ) {
    my $seqIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{ shift() };
    $length =
        $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'subtLength'};
  }
  else {
    foreach my $seq ( @{ $this->{'fastaInfo'}->{'sequences'} } ) {
      $length += $seq->{'subtLength'};
    }
  }

  return $length;

}

##-------------------------------------------------------------------------##

=head2 getIDs()

  Use my @seqIDs = &getIDs( $obj );

  Get an array of all sequence identifiers in the database.

=cut

##-------------------------------------------------------------------------##
sub getIDs {
  my $this = shift;

  my @seqLengths = ();
  return ( @seqLengths ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( map { $_->{'id'} } @{ $this->{'fastaInfo'}->{'sequences'} } );

}

##-------------------------------------------------------------------------##

=head2 getDescriptors()

  Use my @descriptors = &getDescriptors( $obj );

  Get an array of all sequence descriptions contained in the database.

=cut

##-------------------------------------------------------------------------##
sub getDescriptors {
  my $this = shift;

  my @seqLengths = ();
  return ( @seqLengths ) unless ( defined $this->{'fastaInfo'}->{'sequences'} );

  return ( map { $_->{'description'} }
           @{ $this->{'fastaInfo'}->{'sequences'} } );

}

##-------------------------------------------------------------------------##

=head2 getDescription()

  Use: my $desc = &getDescription( $obj, $seqID);

  Get the sequence description text for a given $seqID.

=cut

##-------------------------------------------------------------------------##
sub getDescription {
  my $this  = shift;
  my $seqID = shift;

  return ( -1 )
      unless ( defined $this->{'fastaInfo'}->{'seqIDHash'}->{$seqID} );

  my $index = $this->{'fastaInfo'}->{'seqIDHash'}->{$seqID};
  return $this->{'fastaInfo'}->{'sequences'}->[ $index ]->{'description'};

}

##-------------------------------------------------------------------------##

=head2 compact()

  Use: my &compact();

  Compact database.  Provide a means to remove "holes" or space inefficiencies
  in the database.

=cut

##-------------------------------------------------------------------------##
sub compact {
  my $this = shift;

  print "FastaDB::compact()\n" if ( $DEBUG );

  my $filename = $this->{'fastaInfo'}->{'filename'};
  if ( -s $filename ) {

    # TODO: Put back?
    if ( $this->{'fastaInfo'}->{'openMode'} == SeqDBI::ReadOnly ) {
      &_indexOnly( $this, $filename );
    }
    else {
      &_cleanIndexAndCompact( $this, $filename );
    }
  }
  else {
    croak "FastaDB::compact - Error could not locate file $filename!\n";
  }

}

##-------------------------------------------------------------------------##

=head2 getSequence()

  Use: my $seqString =  &getSequence( $seqID );

  Get the sequence string for a given seqID.

=cut

##-------------------------------------------------------------------------##
sub getSequence {
  my $this  = shift;
  my $seqID = shift;

  my $result = $this->_getFastaRecords( $seqID, $seqID );

  return ( $result->[ 0 ]->{'sequence'} );

}

#sub getSequenceRange{
#  my $this = shift;
#
#  return( &_getFastaRecords( $this, @_ ) );
#
#}

##-------------------------------------------------------------------------##

=head2 getSubstr()

  Use: my $sequence = &getSubstr( $obj, $seqID, $offset, [$length] );

  Get a portion of a sequence from the database.  If the $length parameter
  is left the return value will be everything following the $offset.
  NOTE: $offset is zero-based

=cut

##-------------------------------------------------------------------------##
sub getSubstr {
  my $this   = shift;
  my $seqID  = shift;
  my $offset = shift;
  my $length = shift;

  # Check to see if indices are out of bounds
  my $seqIDLength = &getSeqLength( $this, $seqID );

  # If length is negative then we want all sequence
  # following the offset.
  if ( !defined $length ) {
    $length = $seqIDLength - $offset;
  }
  elsif ( $length <= 0 ) {
    return;
  }

  if ( !defined $seqIDLength || $seqIDLength < 0 ) {
    croak "FastaDB::getSubstr - Error could not locate seqID=$seqID!\n";
  }
  elsif ( $seqIDLength < $offset + $length || $offset < 0 ) {
    croak "FastaDB::getSubstr - Error index out of bounds! "
        . "(SeqID=$seqID, offset=$offset, length=$length actualSeqLen="
        . "$seqIDLength)\n";
  }

  my $fastaDataRec =
      &_getFastaRecords( $this, $seqID, $seqID, $offset,
                         $offset + $length - 1 );

  return ( $fastaDataRec->[ 0 ]->{'sequence'} );

}

##-------------------------------------------------------------------------##

=head2 setSubstr()

  Use: my &setSubstr( $obj, $seqID, $offset, $length, $newSeq);

  Replace a portion of a sequence with a new sequence.  If the
  $newSeq is smaller or larger than the specified length the
  database's sequence will be contracted or expanded accordingly
  ( as the perl substr function ).

=cut

##-------------------------------------------------------------------------##
sub setSubstr {
  my $this   = shift;
  my $seqID  = shift;
  my $offset = shift;
  my $length = shift;
  my $newSeq = shift;

  print $CLASS. "::setSubstr( $seqID, $offset, $length, $newSeq );\n"
      if ( $DEBUG );

  croak $CLASS. "::setSubstr: Attempt to write to a read-only " . "database!\n"
      if ( $this->{'fastaInfo'}->{'openMode'} == SeqDBI::ReadOnly );

  # Check to see if indices are out of bounds
  my $seqIDLength = &getSeqLength( $this, $seqID );
  if ( $seqIDLength < 0 ) {
    croak "FastaDB::substr - Error could not locate seqID=$seqID!\n";
  }
  elsif ( $seqIDLength < ( $offset + $length ) || $offset < 0 ) {
    croak "FastaDB::substr - Error index out of bounds!\n"
        . "(SeqID=$seqID, offset=$offset, length=$length actualSeqLen="
        . "$seqIDLength)\n";
  }
  my $seqIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{$seqID};

  # If length is negative then we want all sequence
  # following the offset.
  if ( $length < 0 ) {
    $length = $seqIDLength - $offset;
  }

  # Figure out rough start line number
  my ( $seqPos, $bytePos, $linePos ) =
      &_getNearestFileIndicesForSeqPos( $this, $seqID, $offset );

  # Save the newSeqLength for later use
  my $newSeqLength = length( $newSeq );

  # parse using Tie_File
  my @lines = ();
  my $file  = $this->{'fastaInfo'}->{'filename'};

  tie @lines, 'Tie::File', $file
      or croak "FastaDB::substr() - Error could not open file $!\n";

  # We now need to determine the real startLine ( and line offset )
  # for our substr offset.
  my $startLine       = -1;
  my $lineStartOffset = -1;
  my $endLine         = -1;
  my $lineEndOffset   = -1;
  my $lineAdj         = 0;

  if ( $seqIDLength > 0 ) {
    while ( $linePos <= $#lines && $endLine < 0 ) {
      $_ = $lines[ $linePos ];
      last if ( /^\>/ );
      s/\s\n\r//g;
      if ( $startLine == -1 && $seqPos + length() - 1 >= $offset ) {

        # We have found the line containing our
        # offset.
        $startLine       = $linePos;
        $lineStartOffset = $offset - $seqPos;
      }
      if (    $endLine == -1
           && $seqPos + length() - 1 >= ( $offset + $length - 1 ) )
      {

        # We have found the line containing our end pos
        $endLine       = $linePos;
        $lineEndOffset = ( $offset + $length ) - $seqPos;
      }
      $linePos++;
      $seqPos += length();
    }
    if ( $offset == $seqIDLength ) {
      ## Special case for appending sequence
      $startLine       = $endLine;
      $lineStartOffset = length( $lines[ $endLine ] );
      $endLine         = $startLine;
      $lineEndOffset   = length( $lines[ $endLine ] );
    }
  }
  else {

    # Add an empty line
    my $recLineNum =
        $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'startLine'};
    splice( @lines, $recLineNum, 0, "" );
    $lineAdj = 1;

    $startLine       = $recLineNum;
    $lineStartOffset = 0;
    $endLine         = $startLine;
    $lineEndOffset   = 0;
  }
  print "FastaDB::setSubstr - startLine = $startLine "
      . "startLineOffset=$lineStartOffset\n            endLine = "
      . "$endLine endLineOffset=$lineEndOffset\n"
      if ( $DEBUG );

  # Calculate stats on the new sequence
  my $newGCLen   = &getGCLength( $this,   $seqID );
  my $newSubtLen = &getSubtLength( $this, $seqID );
  my $newXNLen   = &getXNLength( $this,   $seqID );
  if ( $newSeq ne "" ) {
    $newSubtLen += length( $newSeq ) - ( $newSeq =~ tr/XNRYMK/XNRYMK/ );
    $newXNLen += length( $newSeq );
    while ( $newSeq =~ /([X,N]{20,})/ig ) {
      $newXNLen -= length( $1 );
    }
    $newGCLen += ( $newSeq =~ tr/GC/GC/ );
  }

  $newSeq = substr( $lines[ $startLine ], 0, $lineStartOffset ) . $newSeq;

  if ( $lineEndOffset < length( $lines[ $endLine ] ) ) {
    $newSeq .= substr( $lines[ $endLine ], $lineEndOffset );
  }

  # Should probably just do this with the two lines instead.
  $newSeq =~ s/\n\r\s//g;

  # Create an array of the new lines
  my @newLines = ();
  while ( length( $newSeq ) > $fastaLineLen ) {
    $newSeq =~ s/^(.{$fastaLineLen})//;
    push @newLines, $1;
  }
  if ( length( $newSeq ) > 0 ) {
    push @newLines, $newSeq;
    $newSeq = "";
  }

  # Make the changes
  my @removedLines =
      splice( @lines, $startLine, $endLine - $startLine + 1, @newLines );

  # Prepare return value
  $lineAdj += $#newLines - $#removedLines;
  if ( $startLine == $endLine ) {
    $removedLines[ 0 ] =
        substr( $removedLines[ 0 ], $lineStartOffset, $length );
  }
  else {
    $removedLines[ 0 ] = substr( $removedLines[ 0 ], $lineStartOffset );
    $removedLines[ $#removedLines ] =
        substr( $removedLines[ $#removedLines ], 0, $lineEndOffset );
  }
  my $oldSeq = "";
  if ( @removedLines ) {
    $oldSeq = join( "", @removedLines );
    $oldSeq =~ s/\n\r\s//g;
  }
  undef @removedLines;

  # Adjust stats for removed sequence
  if ( $oldSeq ne "" ) {
    $newSubtLen -= length( $oldSeq ) - ( $oldSeq =~ tr/XNRYMK/XNRYMK/ );
    $newGCLen -= ( $oldSeq =~ tr/GC/GC/ );
    $newXNLen -= length( $oldSeq );
    while ( $oldSeq =~ /([X,N]{20,})/ig ) {
      $newXNLen += length( $1 );
    }
  }

  # Close file and release resources
  untie @lines;

  # Adjust indices and sequence length for change
  my $seqAdj = $newSeqLength - $length;
  $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'length'} += $seqAdj;
  &_adjustFileIndices( $this, $seqID, $offset, $offset + $length - 1,
                       $seqAdj, $lineAdj );

  # Adjust GC and SUBST Lengths
  $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'gcLength'} = $newGCLen;
  $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'subtLength'} =
      $newSubtLen;
  $this->{'fastaInfo'}->{'sequences'}->[ $seqIndex ]->{'XNLength'} = $newXNLen;

  return ( $oldSeq );

}

##-------------------------------------------------------------------------##
## Subroutines and Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
##   my &_adjustFileIndices( $this, $seqID, $startOffset, $endOffset,
##                           $seqAdj, $lineAdj );
##
##  Returns
##
##-------------------------------------------------------------------------##
sub _adjustFileIndices {
  my $this        = shift;
  my $seqID       = shift;
  my $startOffset = shift;
  my $endOffset   = shift;
  my $seqAdj      = shift;
  my $lineAdj     = shift;

  print "FastaDB::_adjustFileIndices( $this, $seqID, $startOffset, "
      . "$endOffset, $seqAdj, $lineAdj)\n"
      if ( $DEBUG );

  my $lineTermByteSize = 1;
  my $byteAdj          = $seqAdj + ( $lineTermByteSize * $lineAdj );

  # Locate the sequence record for $seqID
  my $sequencesArrayRef = $this->{'fastaInfo'}->{'sequences'};
  my $seqIdx            = -1;
  for ( $seqIdx = 0 ; $seqIdx < $#{$sequencesArrayRef} ; $seqIdx++ ) {
    next if ( $sequencesArrayRef->[ $seqIdx ]->{'id'} ne $seqID );
    last;
  }
  croak "FastaDB::_adjustFileIndices - Error..could not locate sequence "
      . "record for seqID=$seqID!\n"
      if ( $sequencesArrayRef->[ $seqIdx ]->{'id'} ne $seqID );
  my $seqLen = $sequencesArrayRef->[ $seqIdx ]->{'length'};

  # Find the first adjustment position and fix
  my $posIndicesArrayRef = $sequencesArrayRef->[ $seqIdx ]->{'seqPosIndices'};
  for ( my $i = $#{$posIndicesArrayRef} ; $i >= 0 ; $i -= 3 ) {
    if ( $posIndicesArrayRef->[ $i - 2 ] >= $startOffset ) {
      if (
          (
               $seqAdj < 0
            && $posIndicesArrayRef->[ $i - 2 ] <=
            ( abs( $seqAdj ) + $startOffset )
          )
          || ( $posIndicesArrayRef->[ $i - 2 ] <= $endOffset )
          || ( ( $posIndicesArrayRef->[ $i - 2 ] + $seqAdj ) > ( $seqLen - 1 ) )
          )
      {
        splice( @{$posIndicesArrayRef}, $i - 2, 3 );
      }
      else {
        $posIndicesArrayRef->[ $i - 2 ] += $seqAdj;
        $posIndicesArrayRef->[ $i - 1 ] += $byteAdj;
        $posIndicesArrayRef->[ $i ]     += $lineAdj;
      }
    }
    else {
      last;
    }
  }
  $seqIdx++;

  # Fix remaining sequences in the file
  for ( my $j = $seqIdx ; $j <= $#{$sequencesArrayRef} ; $j++ ) {
    $sequencesArrayRef->[ $j ]->{'startByte'} += $byteAdj;
    $sequencesArrayRef->[ $j ]->{'startLine'} += $lineAdj;
    $posIndicesArrayRef = $sequencesArrayRef->[ $j ]->{'seqPosIndices'};
    for ( my $i = 0 ; $i <= $#{$posIndicesArrayRef} ; $i += 3 ) {
      $posIndicesArrayRef->[ $i + 1 ] += $byteAdj;
      $posIndicesArrayRef->[ $i + 2 ] += $lineAdj;
    }
  }
}

##-------------------------------------------------------------------------##
##   my &_cleanIndexAndCompact( $obj, $file );
##
##      $file:            The name of the input file containing 1-N
##                        sequences in FASTA format.
##
##  Returns
##
##    Modifies the original file (in-place) in preparation
##    for using it as a database.  The three operations performed
##    on the file are:
##
##     1. Clean - Attempts are made to fix common FASTA format
##        errors including: line termination problems, duplicate sequence
##        names and sequence and header partially appearing on the same
##        line. Multiple sequence records appearing on the same line
##        cause a warning to be printed to STDERR and are skipped
##        (left intact but not indexed for use by the accessor methods).
##        Invalid text in a record cause a similar warning and are
##        also skipped.
##
##     2. Index - Each valid FASTA header is indexed by line number,
##        and byte position.  In addition each FASTA sequence is
##        indexed by sequence position, line number, and byte position.
##        The indexes refer to the first character at the start of
##        a line.  The space between indexes is controlled by setting
##        the indexDistance parameter at object creation time.
##
##     3. Compact - The sequence lines in the FASTA file are
##        reorganized so that each line is maximally N number
##        of bases wide.
##
##    This routine is also responsible for building the
##    fastaInfo datastructure...which includes the following
##    properties of each sequence:
##
##    $length:   The total sequence length of this file.
##    $subtLen:  The total sequence length minus substitution characters.
##    $gcLen:    The total "G" and "C" bases in this file.
##    $xnLen:    The total sequence length minus "x" and "N" characters.
##
##------------------------------------------------------------------------##
sub _cleanIndexAndCompact {
  my $this     = shift;
  my $filename = shift;

  print "FastaDB::_cleanIndexAndCompact()\n" if ( $DEBUG );

  open IN, "+<$filename"
      or croak "FastaDB::_cleanIndexAndCompact() - Error could not open "
      . "file $filename: $!\n";

  my $writeDelaySize = 4096;

  my $seqLine        = "";
  my $writePtr       = 0;
  my $dataOut        = "";
  my $filePos        = 0;
  my $bytesToWrite   = 0;
  my $seekRelative   = 0;
  my $newLineNum     = 0;
  my $newBytePos     = 0;
  my $seqLen         = 0;
  my $badSequence    = 0;
  my $gcLen          = 0;
  my $xnLen          = 0;
  my $subtLen        = 0;
  my $fastaInfoRec   = undef;
  my @fastaInfoRecs  = ();
  my $seqName        = "";
  my $origName       = "";
  my $description    = "";
  my %fastaNameSpace = ();
  my $newHeader      = "";
  my $sameCount      = 1;
  while ( <IN> ) {
    next unless /\S/;    # Eat blank lines

RESTART:

    # Look for sequence headers
    if ( /^\s*\>(.*)/ ) {

      # Fix problems with empty ">" specifiers
      $_ = $1;
      if ( /^\s*$/ ) {
        $_ = "UnnamedSeq";
      }

      if ( $seqLine ) {

        # Place remaining sequence into output buffer
        $dataOut .= "$seqLine\n";

        # This is the first byte position of the new record
        $newBytePos += length( $seqLine ) + 1;
        $seqLine = "";

        # This is the first line of the new record
        $newLineNum++;
      }
      ## Save previous record seqlen here
      if ( $badSequence ) {
        $badSequence = 0;
        pop @fastaInfoRecs;
      }
      else {
        if ( @fastaInfoRecs ) {
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'length'}     = $seqLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'gcLength'}   = $gcLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'xnLength'}   = $xnLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'subtLength'} = $subtLen;
        }
      }
      $seqLen  = 0;
      $gcLen   = 0;
      $xnLen   = 0;
      $subtLen = 0;
      if ( /^(.*?)([ACGTNXacgtnx]{25,})$/ ) {
        print STDERR "It appears as if the fasta header line \n$_\n contains "
            . "the beginning of a sequence.  The program will "
            . "incorporate this into the remaining sequence lines. "
            . "Interrupt with ^c and add a character other "
            . "than [ACGTNXacgtnx] or a space to the end of the "
            . "header line to override this behavior.\n";

        # Sequence on same line as header
        $_ = $1;
        my $ucSeq = uc( $2 );
        $seqLine = $ucSeq;
        my $ucLen = length( $ucSeq );
        $seqLen  += $ucLen;
        $subtLen += $ucLen;
        $xnLen   += $ucLen;
        while ( $ucSeq =~ /([X,N]{20,})/ig ) {
          $xnLen -= length( $1 );
        }
        $subtLen -= ( $ucSeq =~ tr/XNRYMK// );
        $gcLen += ( $ucSeq =~ tr/GC// );
      }

      if ( /\>.*/ && !/\>(\d+)\%/ ) {
        print STDERR "FastaDB:_cleanIndexAndCompact - ERROR: Multiple fasta "
            . "seqs appear on one line ($_) and possibly "
            . "more! Ignoring all entries on this line.\n";
        $dataOut .= ">" . $_ . "\n";
        $newBytePos += length( $_ ) + 2;
        $newLineNum++;
        next;
      }

      ## Save new record info here
      /\s*(\S+)\s*(.*)/;
      $seqName = $1;
      if ( defined $this->{'fastaInfo'}->{'maxIDLength'}
           && length( $seqName ) > $this->{'fastaInfo'}->{'maxIDLength'} )
      {
        croak
"FastaDB::_cleanIndexAndCompact(): Fasta file contains a sequence identifier "
            . "which is too long ( max id length = "
            . $this->{'fastaInfo'}->{'maxIDLength'} . " )\n";
      }
      $description = $2;
      $origName    = $seqName;
      if ( $seqName =~ /\.seq$/ ) {
        $seqName =~ s/\.seq$/_seqqqq_q/;
      }

      # Create a uniq sequence id
      if ( defined( $fastaNameSpace{$seqName} ) ) {
        $sameCount = 1;
        print STDERR "\nThe name $origName is used more than once in the "
            . "fasta file. Second and later occurrences are appended a "
            . "_number\n\n";
        while ( defined( $fastaNameSpace{ "$seqName" . "_$sameCount" } ) ) {
          ++$sameCount;
        }
        $seqName .= "_$sameCount";
      }
      $fastaNameSpace{$seqName} = $#fastaInfoRecs + 1;

      $newHeader = ">$seqName";
      $newHeader .= " " . $description if ( $description ne "" );
      $newHeader .= "\n";
      $dataOut   .= $newHeader;

      # This is the first byte position of the sequence
      $newBytePos += length( $newHeader );

      # This is the first line of the sequence in the record
      $newLineNum++;

      $fastaInfoRec = {
                        'id'            => $seqName,
                        'name'          => $seqName,
                        'description'   => $description,
                        'length'        => 'Unknown',
                        'subtLength'    => 'Unknown',
                        'gcLength'      => 'Unknown',
                        'xnLength'      => 'Unknown',
                        'startLine'     => $newLineNum,
                        'startByte'     => $newBytePos,
                        'seqPosIndices' => []
      };

      ## Save first sequence index here
      push @{ $fastaInfoRec->{'seqPosIndices'} },
          ( 0, $newBytePos, $newLineNum );
      push @fastaInfoRecs, $fastaInfoRec;
      next;
    }

    # Get rid of useless characters
    #   Translates everything but true printable ascii characters
    #   to spaces.  Then we get rid of spaces
    tr/\040-\176/ /c;
    s/[\s\n\r]//g;
    $_ = uc();

    # Check for a data-rich sequence line
    if ( /^[ACGTBDHVRYKMSWNX]+$/ ) {

      # If the first sequence doesn't contain a header
      # make one!
      if ( !@fastaInfoRecs ) {
        $_ = ">UnnamedSeq $_";
        goto RESTART;
      }
      my $tmpSeq = $_;
      $seqLine .= $tmpSeq;
      my $inputLen = length( $tmpSeq );
      $seqLen  += $inputLen;
      $subtLen += $inputLen;
      $xnLen   += $inputLen;
      while ( $tmpSeq =~ /([X,N]{20,})/ig ) {
        $xnLen -= length( $1 );
      }
      $subtLen -= ( $tmpSeq =~ tr/XNRYMK/XNRYMK/ );
      $gcLen += ( $tmpSeq =~ tr/GC/GC/ );
      $tmpSeq = "";
      my $lineLen  = length( $seqLine );
      my $block    = "";
      my $blockLen = 0;

      while ( $lineLen >= $fastaLineLen ) {
        $block = substr( $seqLine, 0, $fastaLineLen );
        substr( $seqLine, 0, $fastaLineLen ) = "";
        $dataOut .= "$block\n";
        $blockLen = length( $block );
        $newBytePos += $blockLen + 1;
        $lineLen -= $blockLen;
        $newLineNum++;
      }
    }
    elsif ( /^([ACGTBDHVRYKMSWNX]+)(\>.*)$/ ) {

      # If the first sequence doesn't contain a header
      # make one!
      if ( !@fastaInfoRecs ) {
        $_ = ">UnnamedSeq $_";
        goto RESTART;
      }

      # A missing line termination - we can correct this one
      # Note the following code is a duplicate of the above
      # for efficiency reasons....do not want to run the large
      # regular expression through the whole file.
      my $tmpSeq = $1;
      my $suffix = $2;
      $seqLine .= $tmpSeq;
      my $inLen = length( $tmpSeq );
      $seqLen  += $inLen;
      $subtLen += $inLen;
      $xnLen   += $inLen;

      while ( $tmpSeq =~ /([X,N]{20,})/ig ) {
        $xnLen -= length( $1 );
      }
      $subtLen -= ( $tmpSeq =~ tr/XNRYMK/XNRYMK/ );
      $gcLen += ( $tmpSeq =~ tr/GC/GC/ );
      $_      = $suffix;
      $tmpSeq = "";
      $suffix = "";
      my $lineLen  = length( $seqLine );
      my $block    = "";
      my $blockLen = 0;

      while ( $lineLen >= $fastaLineLen ) {
        $block .= substr( $seqLine, 0, $fastaLineLen );
        substr( $seqLine, 0, $fastaLineLen ) = "";
        $dataOut .= "$block\n";
        $blockLen = length( $block );
        $newBytePos += $blockLen + 1;
        $lineLen -= $blockLen;
        $newLineNum++;
      }
      goto RESTART;
    }
    else {

      # Unrecognized line!
      print STDERR "FastaDB::_cleanIndexAndCompact - WARNING: RepeatMasker "
          . "encountered a line in an unrecognized format.";

      if ( /([ACGTBDHVRYKMSWNX\s]{0,5}[^ACGTBDHVRYKMSWNX\s>]+.{0,5})/ ) {
        print STDERR "  The offending subsequence is \"$1\". ";
      }

      print STDERR
          "The offending line is \"$_\".  Sequence \"$seqName\" is being "
          . "ignored.\n";
      $badSequence = 1;
    }

    $filePos = tell( IN );

    if ( length( $dataOut ) >= $writeDelaySize
         && ( $filePos - $writePtr ) >= $writeDelaySize )
    {
      if ( ( $filePos - $writePtr ) > length( $dataOut ) ) {
        $bytesToWrite = length( $dataOut );
      }
      else {
        $bytesToWrite = ( $filePos - $writePtr );
      }
      if ( $bytesToWrite ) {
        $seekRelative = $filePos - $writePtr;
        seek IN, -$seekRelative, 1;
        print IN "" . substr( $dataOut, 0, $bytesToWrite );
        substr( $dataOut, 0, $bytesToWrite ) = "";
        seek IN, $seekRelative - $bytesToWrite, 1;
        $writePtr += $bytesToWrite;
      }
    }
  }

  #
  # Trailing Cases
  #
  if ( $badSequence ) {
    pop @fastaInfoRecs;
  }
  else {
    if ( @fastaInfoRecs ) {
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'length'}     = $seqLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'gcLength'}   = $gcLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'xnLength'}   = $xnLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'subtLength'} = $subtLen;
    }
  }

  $dataOut .= $seqLine;
  if ( $dataOut ) {

    # The file can grow due to the addition of new
    # line termination characters.
    seek IN, -( tell( IN ) - $writePtr ), 1;
    print IN $dataOut;
    print IN "\n";
  }

  #   The file may have become shorter due to concatenation
  #   of lines or removal of whitespace.  Adjust the file
  #   size using the truncate command.
  truncate IN, tell( IN );

  # Release file
  close IN;

  # Add indexes for quick retreival
  my $inc = int( $seqIndexDist / $fastaLineLen ) * $fastaLineLen;
  if ( $inc > 0 ) {
    foreach $fastaInfoRec ( @fastaInfoRecs ) {
      my $startBytePos = $fastaInfoRec->{'seqPosIndices'}->[ 1 ];
      my $startLineNum = $fastaInfoRec->{'seqPosIndices'}->[ 2 ];
      for ( my $i = $inc ; $i <= $fastaInfoRec->{'length'} ; $i += $inc ) {
        $newBytePos = $startBytePos + $i + ( $i / $fastaLineLen );
        $newLineNum = $startLineNum + ( $i / $fastaLineLen );
        push @{ $fastaInfoRec->{'seqPosIndices'} },
            ( $i, $newBytePos, $newLineNum );
      }
    }
  }

  # We are outa here!
  $this->{'fastaInfo'}->{'sequences'} = \@fastaInfoRecs;
  $this->{'fastaInfo'}->{'seqIDHash'} = \%fastaNameSpace;

}

##-------------------------------------------------------------------------##
##   my &_indexOnly( $obj, $file );
##
##      $file:            The name of the input file containing 1-N
##                        sequences in FASTA format.
##
##  Returns
##
##   Each valid FASTA header is indexed by line number,
##   and byte position.  In addition each FASTA sequence is
##   indexed by sequence position, line number, and byte position.
##   The indexes refer to the first character at the start of
##   a line.  The space between indexes is controlled by setting
##   the indexDistance parameter at object creation time.
##
##   This routine is also responsible for building the
##   fastaInfo datastructure...which includes the following
##   properties of each sequence:
##
##   $length:   The total sequence length of this file.
##   $subtLen:  The total sequence length minus substitution characters.
##   $gcLen:    The total "G" and "C" bases in this file.
##   $xnLen:    The total sequence length minus 'X' and 'N' characters.
##
##------------------------------------------------------------------------##
sub _indexOnly {
  my $this     = shift;
  my $filename = shift;

  print "FastaDB::_indexOnly()\n" if ( $DEBUG );

  open IN, "<$filename"
      or croak "FastaDB::_indexOnly() - Error could not open "
      . "file $filename: $!\n";

  my $seqLine        = "";
  my $filePos        = 0;
  my $bytesToWrite   = 0;
  my $seekRelative   = 0;
  my $newLineNum     = 0;
  my $newBytePos     = 0;
  my $seqLen         = 0;
  my $badSequence    = 0;
  my $gcLen          = 0;
  my $xnLen          = 0;
  my $subtLen        = 0;
  my $fastaInfoRec   = undef;
  my @fastaInfoRecs  = ();
  my $seqName        = "";
  my $origName       = "";
  my $description    = "";
  my %fastaNameSpace = ();
  my $newHeader      = "";
  my $sameCount      = 1;
  my $indexLen       = 0;
  while ( <IN> ) {

    my $lineLen = length( $_ );

    # This is the current "start_of_line" byte position (0-based)
    $newBytePos = tell( IN ) - length( $_ );

    # This is the current line counter; (1 based)
    $newLineNum++;

    next unless /\S/;    # Eat blank lines

RESTART:

    # Look for sequence headers
    if ( /^\s*\>(.*)/ ) {
      $_ = $1;
      ## Save previous record seqlen here
      if ( $badSequence ) {
        $badSequence = 0;
        pop @fastaInfoRecs;
      }
      else {
        if ( $seqLen ) {
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'length'}     = $seqLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'gcLength'}   = $gcLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'xnLength'}   = $xnLen;
          $fastaInfoRecs[ $#fastaInfoRecs ]->{'subtLength'} = $subtLen;
        }
        else {

          # Empty record
          pop @fastaInfoRecs;
        }
      }
      $indexLen = 0;
      $seqLen   = 0;
      $gcLen    = 0;
      $xnLen    = 0;
      $subtLen  = 0;
      if ( /^(.*?)([ACGTNXacgtnx]{25,})$/ ) {
        print STDERR "It appears as if the fasta header line \n$_\n contains "
            . "the beginning of a sequence.  This sequence will be "
            . "ignored.\n";

        $badSequence = 1;
      }

      if ( /\>.*/ && !/\>(\d+)\%/ ) {
        print STDERR "FastaDB:_cleanIndexAndCompact - ERROR: Multiple fasta "
            . "seqs appear on one line ($_) and possibly "
            . "more! Ignoring all entries on this line.\n";
        $badSequence = 1;
      }

      ## Fix problems with empty ">" specifiers
      if ( /^\s*$/ ) {
        $_ = "UnnamedSeq";
      }

      ## Save new record info here
      /\s*(\S+)\s*(.*)/;
      $seqName = $1;
      if ( defined $this->{'fastaInfo'}->{'maxIDLength'}
           && length( $seqName ) > $this->{'fastaInfo'}->{'maxIDLength'} )
      {
        croak
"FastaDB::_cleanIndexAndCompact(): Fasta file contains a sequence identifier "
            . "which is too long ( max id length = "
            . $this->{'fastaInfo'}->{'maxIDLength'} . " )\n";
      }
      $description = $2;
      $origName    = $seqName;
      if ( $seqName =~ /\.seq$/ ) {
        $seqName =~ s/\.seq$/_seqqqq_q/;
      }

      # Create a uniq sequence id
      if ( defined( $fastaNameSpace{$seqName} ) ) {
        $sameCount = 1;
        print STDERR "\nThe name $origName is used more than once in the "
            . "fasta file. Second and later occurrences are appended a "
            . "_number\n\n";
        while ( defined( $fastaNameSpace{ "$seqName" . "_$sameCount" } ) ) {
          ++$sameCount;
        }
        $seqName .= "_$sameCount";
      }
      $fastaNameSpace{$seqName} = $#fastaInfoRecs + 1;

      $fastaInfoRec = {
                        'id'            => $seqName,
                        'name'          => $seqName,
                        'description'   => $description,
                        'length'        => 'Unknown',
                        'subtLength'    => 'Unknown',
                        'gcLength'      => 'Unknown',
                        'xnLength'      => 'Unknown',
                        'startLine'     => $newLineNum,
                        'startByte'     => $newBytePos + $lineLen,
                        'seqPosIndices' => []
      };

      ## Save first sequence index here
      ##push @{$fastaInfoRec->{'seqPosIndices'}},
      ##                           ( 0, $newBytePos, $newLineNum );
      push @fastaInfoRecs, $fastaInfoRec;
      next;
    }

    # Check for a data-rich sequence line
    if ( /^[ACGTBDHVRYKMSWNXacgtbdhvrykmswnx\s]+$/ ) {

      # If the first sequence doesn't contain a header
      # make one!
      if ( !@fastaInfoRecs ) {
        $fastaInfoRec = {
                          'id'            => "UnnamedSeq",
                          'name'          => "UnnamedSeq",
                          'description'   => "",
                          'length'        => 'Unknown',
                          'subtLength'    => 'Unknown',
                          'gcLength'      => 'Unknown',
                          'xnLength'      => 'Unknown',
                          'startLine'     => -1,
                          'startByte'     => -1,
                          'seqPosIndices' => []
        };

        push @fastaInfoRecs, $fastaInfoRec;

      }
      s/[\s]//g;
      my $tmpSeq  = $_;
      my $lineLen = length( $tmpSeq );
      $seqLen   += $lineLen;
      $subtLen  += $lineLen;
      $indexLen += $lineLen;
      $xnLen    += $lineLen;
      while ( $tmpSeq =~ /([X,N]{20,})/ig ) {
        $xnLen -= length( $1 );
      }
      $subtLen -= ( $tmpSeq =~ tr/XNRYMK/XNRYMK/ );
      $gcLen += ( $tmpSeq =~ tr/GC/GC/ );
      $tmpSeq = "";
      if ( $indexLen >= $seqIndexDist ) {
        $indexLen = 0;
        push @{ $fastaInfoRecs[ $#fastaInfoRecs ]->{'seqPosIndices'} },
            ( ( $seqLen - $lineLen ), $newBytePos - 1, $newLineNum - 1 );

      }
    }
    elsif ( /^([ACGTBDHVRYKMSWNXacgtbdhvrykmswnx]+)(\>.*)$/ ) {
      print STDERR "FastaDB::_indexOnly - Sequence with \> in the "
          . "sequence itself.  Probably missing a line term. "
          . "Ignoring this sequence!\n";
      $badSequence = 1;
    }
    else {

      # Unrecognized line!
      print STDERR "FastaDB::_indexOnly - WARNING: This "
          . "line doesn't appear to be in a format I recognize "
          . "( $_ ). Sequence is being ignored.\n";
      $badSequence = 1;
    }

  }

  #
  # Trailing Cases
  #
  if ( $badSequence ) {
    pop @fastaInfoRecs;
  }
  else {
    if ( $seqLen ) {
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'length'}     = $seqLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'gcLength'}   = $gcLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'xnLength'}   = $xnLen;
      $fastaInfoRecs[ $#fastaInfoRecs ]->{'subtLength'} = $subtLen;
    }
    else {

      # Empty record
      pop @fastaInfoRecs;
    }
  }

  # Release file
  close IN;

  # We are outa here!
  $this->{'fastaInfo'}->{'sequences'} = \@fastaInfoRecs;
  $this->{'fastaInfo'}->{'seqIDHash'} = \%fastaNameSpace;

}

##-------------------------------------------------------------------------##
##   my ($seqPos, $bytePos, $lineNum) =
##                     &_getNearestFileIndicesForSeqPos( $obj, $seqID,
##                                                       $startPos );
##  Returns
##
##-------------------------------------------------------------------------##
sub _getNearestFileIndicesForSeqPos {
  my $this     = shift;
  my $seqID    = shift;
  my $startPos = shift;

  print "FastaDB::_getNearestFileIndicesForSeqPos( $this, $seqID, $startPos "
      . ");\n"
      if ( $DEBUG );

  croak "FastaDB::_getNearestFileIndicesForSeqPos(): SeqID $seqID does not "
      . "exist!\n"
      unless ( defined $this->{'fastaInfo'}->{'seqIDHash'}->{$seqID} );

  my $bytePos        = 0;
  my $seqPos         = -1;
  my $linePos        = -1;
  my $fastaInfo      = $this->{'fastaInfo'}->{'sequences'};
  my $fastaInfoIndex = $this->{'fastaInfo'}->{'seqIDHash'}->{$seqID};
  my $seqRecord      = $fastaInfo->[ $fastaInfoIndex ];

  $seqPos  = 0;
  $bytePos = $seqRecord->{'startByte'};
  $linePos = $seqRecord->{'startLine'};

  my $approxIndexIndex = ( int( $startPos / $seqIndexDist ) - 1 ) * 3;

  $approxIndexIndex = $#{ $seqRecord->{'seqPosIndices'} } - 2
      if ( $approxIndexIndex > $#{ $seqRecord->{'seqPosIndices'} } );

  if ( $approxIndexIndex >= 0 ) {

    if ( $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ] > $startPos ) {

      # search down
      while ($approxIndexIndex >= 0
          && $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ] >= $startPos )
      {
        $approxIndexIndex -= 3;
      }
    }
    elsif ( $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ] < $startPos ) {

      # search up
      while ($approxIndexIndex <= $#{ $seqRecord->{'seqPosIndices'} }
          && $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ] <= $startPos )
      {
        $approxIndexIndex += 3;
      }
      if (    $approxIndexIndex > $#{ $seqRecord->{'seqPosIndices'} }
           || $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ] > $startPos )
      {
        $approxIndexIndex -= 3;
      }
    }

    if (    $approxIndexIndex >= 0
         && $approxIndexIndex <= $#{ $seqRecord->{'seqPosIndices'} } - 2 )
    {
      $seqPos  = $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex ];
      $bytePos = $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex + 1 ];
      $linePos = $seqRecord->{'seqPosIndices'}->[ $approxIndexIndex + 2 ];
    }
  }

  print "FastaDB::_getNearestFileIndicesForSeqPos - Returning ($seqPos, "
      . "$bytePos, $linePos );\n"
      if ( $DEBUG );
  return ( $seqPos, $bytePos, $linePos );
}

##-------------------------------------------------------------------------##
##   my (\@fastaDataRecord) = &_getFastaRecords( $obj, $startSeqID,
##                                               $lastSeqID );
## or
##
##   my (\@fastaDataRecord) = &_getFastaRecords( $obj, $startSeqID,
##                                               $lastSeqID, $startPos,
##                                               $lastPos );
##  Returns
##
##-------------------------------------------------------------------------##
sub _getFastaRecords {
  my $this = shift;

  # Initial defines
  my $startSeqID = undef;    # The first sequence to pull data from.
  my $startPos   = undef;    # The first base position within $startSeqID
                             #  to start extracting from.
  my $lastSeqID  = undef;    # The last sequence (inclusive) to pull from.
  my $lastPos    = undef;    # The last base position (inclusive) to include.
  my $fastaSeqsArrayRef = undef;    # The fasta sequences datastructure
  my $fastaSeqIDHashRef = undef;    # The fasta sequences datastructure
  my $file = $this->{'fastaInfo'}->{'filename'};

  #
  # Determine how this subroutine was called
  #
  if ( @_ == 4 || @_ == 2 ) {
    if ( defined $this->{'fastaInfo'}
         && @{ $this->{'fastaInfo'}->{'sequences'} } > 0 )
    {
      $fastaSeqsArrayRef = $this->{'fastaInfo'}->{'sequences'};
      $fastaSeqIDHashRef = $this->{'fastaInfo'}->{'seqIDHash'};
      if ( @_ == 4 ) {
        $startSeqID = shift;
        $lastSeqID  = shift;
        $startPos   = shift;
        $lastPos    = shift;
        print "FastaDB::_getFastaRecords( $this, $startSeqID, "
            . "$lastSeqID, $startPos, $lastPos );\n"
            if ( $DEBUG );
      }
      else {
        $startSeqID = shift;
        $lastSeqID  = shift;
        print "FastaDB::_getFastaRecords( $this, $startSeqID, "
            . "$lastSeqID );\n"
            if ( $DEBUG );
        $startPos = 0;
        if (  defined $fastaSeqIDHashRef->{$lastSeqID}
           && defined $fastaSeqsArrayRef->[ $fastaSeqIDHashRef->{$lastSeqID} ] )
        {
          $lastPos =
              $fastaSeqsArrayRef->[ $fastaSeqIDHashRef->{$lastSeqID} ]
              ->{'length'};
        }
        else {
          croak "FastaDB::_getFastaRecords: Could not locate "
              . "lastSeqID=$lastSeqID!";
        }
      }
    }
    else {
      croak "FastaDB::_getFastaRecords: Cannot extract sequences "
          . "from empty file!  Called by: _parseFastaFile( obj, "
          . "$startSeqID, $lastSeqID, $startPos, $lastPos );";
    }
  }
  else {
    croak "Usage: FastaDB::_getFastaRecords( \$obj, "
        . "[\$startSeqId, \$lastSeqId, \$startPos, \$lastPos] );";
  }

  # Open the file
  open( IN, $file )
      || croak "FastaDB::_getFastaRecords: Error could not "
      . "open file $file\n";

  ##
  ## Optimization
  ##
  ##  We have a byte-index for sequence lines in the
  ##  fastaInfo record.  This helps us seek to the
  ##  correct location in the file and maybe even
  ##  inside the sequence quickly.
  ##
  my ( $seqPos, $bytePos, $linePos ) =
      &_getNearestFileIndicesForSeqPos( $this, $startSeqID, $startPos );
  if ( $seqPos >= 0 && $bytePos > 0 ) {
    seek IN, $bytePos, 0;
  }

  my $sequence    = "";
  my $seqLen      = 0;
  my $seqName     = $startSeqID;
  my $description = "Unknown";
  if (    defined $fastaSeqIDHashRef->{$startSeqID}
       && defined $fastaSeqsArrayRef->[ $fastaSeqIDHashRef->{$startSeqID} ] )
  {
    $description =
        $fastaSeqsArrayRef->[ $fastaSeqIDHashRef->{$startSeqID} ]
        ->{'description'};
  }
  my @fastaData = ();
  while ( <IN> ) {

    # Better than chomp.
    s/[\n\r]//g;

    # Eat blank lines
    next unless /\S/;

    if ( /^\s*\>\s*(\S+)\s*(.*)/ ) {
      last if ( $seqName eq $lastSeqID );
      ## Could use this to verify we are in the
      ## right place.
      push @fastaData,
          {
            'id'          => $seqName,
            'description' => $description,
            'sequence'    => $sequence
          }
          if ( $seqLen );
      $seqName     = $1;
      $description = $2;
      $seqLen      = 0;
      $sequence    = "";
    }
    elsif ( /^([ACGTBDHVRYKMSWNXacgtbdhvrykmswnx]+)$/ ) {
      $seqLen += length( $1 );
      $sequence .= $1;
      if ( $seqName eq $lastSeqID ) {
        last
            if (
                 $seqLen > $lastPos
                 || ( $startSeqID eq $lastSeqID
                      && ( $seqLen + $seqPos ) > $lastPos )
            );
      }
    }
    else {
      print STDERR "FastaDB::_getFastaRecords: Error could not interpret "
          . "fasta line correctly ( $_ )! Check data before "
          . "proceeding!\n";
      ## TODO: Do some error checking.
    }
  }
  close IN;

  #
  # Trailing Case
  #
  push @fastaData,
      {
        'id'          => $seqName,
        'description' => $description,
        'sequence'    => $sequence
      }
      if ( $seqLen );

  # Truncate data starting with first record
  if ( @fastaData ) {
    if ( $startSeqID eq $lastSeqID ) {
      $fastaData[ 0 ]->{'sequence'} = substr( $fastaData[ 0 ]->{'sequence'},
                                              $startPos - $seqPos,
                                              $lastPos - $startPos + 1 );
    }
    else {
      $fastaData[ 0 ]->{'sequence'} =
          substr( $fastaData[ 0 ]->{'sequence'}, $startPos - $seqPos );
      $fastaData[ $#fastaData ]->{'sequence'} =
          substr( $fastaData[ $#fastaData ]->{'sequence'}, 0, $lastPos + 1 );
    }
  }

  print "FastaDB::_getFastaRecords() Returning..\n" if ( $DEBUG );
  return ( \@fastaData );

}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##        $filename     : A filename to be created
##
##  Returns
##
##      Uses the Data::Dumper module to save out the data
##      structure as a text file.  This text file can be
##      read back into an object of this type.
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

1;
