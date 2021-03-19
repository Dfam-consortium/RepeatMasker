#!/tools/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SimpleBatcher.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      A module for creating simple batches of sequences from a
##      large seqDB.  This implementation will batch complete
##      sequences until it reaches the limit specified by the user.
##      If a single (larger than limit) sequence is encountered it
##      is fragmented with a user specified overlap. Fragment batches
##      are not mixed with complete sequence batches.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
# ChangeLog
#
#   $Log$
#
################################################################################
#
## To Do:
#
#  - Remove transitional datastructures for repeatmasker
#      getSeqNameArray( $batchNum );
#      getDescriptionHash( $batchNum );
#      getSeqWithNameHash( $batchNum );
#
#
#

=head1 NAME

SimpleBatcher - A perl module which assists in the packaging of smaller
                sequence units from a larger seqDB.

=head1 SYNOPSIS

use SimpleBatcher;

  my $seqDBObj = FastaDB->new( fileName => "/usr/lib/seq/test.fa" );
  .. 
  my $batcher = SimpleBatcher->new( $seqDBObj, 
                                    $fragmentLength, 
                                    $overlapLen );

=head1 SEE ALSO

=over 4

SeqDBI, FastaDB, RepeatMasker

=back

=head1 COPYRIGHT

Copyright 2004 Robert Hubley Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence and Package Definition
#
package SimpleBatcher;
use strict;
use Carp;
use Data::Dumper;
use SeqDBI;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw(min max);

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.5;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
##   Create an empty SimpleBatcher suitable for serializing into:
##       SimpleBatcher->new();
##   Creates a FastaBatcher for a fasta file:
##       SimpleBatcher->new( $seqDBObj, $fragmentLength, $overlapLen );
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  croak "Usage: \$batcher = SimpleBatcher->new( [\$seqDBObj, "
      . "\$fragmentLength, \$overlapLength] );"
      if ( @_ != 0 && @_ != 3 );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Process seqDB object if passed to us.
  if ( @_ == 3 ) {
    my $seqDBObj   = shift();
    my $fragLen    = shift();
    my $overlapLen = shift();
    $this->{'fragmentLength'} = $fragLen;
    $this->{'overlapLength'}  = $overlapLen;
    $this->{'seqDBObj'}       = $seqDBObj;
    &_packBatches( $this, $seqDBObj, $fragLen, $overlapLen );
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Public Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 writeBatchFile()

  Use:  $obj->writeBatchFile( $batchNum, $batchFilename );

  NOTE: Batches start at 1
=cut

##-------------------------------------------------------------------------##
sub writeBatchFile {
  my $this          = shift;
  my $batchNum      = shift;
  my $batchFilename = shift;

  return ( -1 )
      unless ( $batchNum > 0
              && defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );

  my $seqDBObj = $this->{'seqDBObj'};

  open OUT, ">$batchFilename"
      || croak "SimpleBatcher::writeBatchFile( "
      . "$batchNum, $batchFilename ): Error "
      . "opening file for writing.";
  my $seqsRef =
      $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ]->{'sequences'};

  ## This is a general case routine.  It could be optimized
  ## by taking advantage of our knowledge of the packBatches
  ## routine.  Since this routine always creates batches which
  ## contain sequences in FASTA file order we could call
  ## _parseFasta once with the first and last seqIDs and positions.
  foreach my $seqRef ( @{$seqsRef} ) {
    my $seqID = $seqRef->{'fastaID'};
    my $seq;
    if ( defined $seqRef->{'startPos'} ) {
      $seq =
          $seqDBObj->getSubstr( $seqID, $seqRef->{'startPos'},
                             $seqRef->{'lastPos'} - $seqRef->{'startPos'} + 1 );
    }
    else {
      $seq = $seqDBObj->getSequence( $seqID );
    }
    print OUT ">" . $seqRef->{'id'};
    my $desc = $seqDBObj->getDescription( $seqID );
    if ( $desc ne "" ) {
      print OUT " " . $desc;
    }
    print OUT "\n";
    $seq =~ s/(\S{50})/$1\n/g;
    $seq .= "\n"
        unless ( $seq =~ /.*\n+$/s );
    print OUT $seq;
  }
  close OUT;

  return ( 1 );
}

##-------------------------------------------------------------------------##

=head2 getSeqDBObj()

  Use: my $seqDBObj = $obj->getSeqDBObj();

  The SeqDB instance which was set when this batcher was
  created.

=cut

##-------------------------------------------------------------------------##
sub getSeqDBObj {
  my $this = shift;

  return $this->{'seqDBObj'};

}

##-------------------------------------------------------------------------##

=head2 translateBatchSeqPositionToFastaSeq()

  Use: $obj->translateBatchSeqPositionToFastaSeq( $seqID, $position );

=cut

##-------------------------------------------------------------------------##
sub translateBatchSeqPositionToFastaSeq {
  my $this     = shift;
  my $seqID    = shift;
  my $position = shift;

  return ( -1 )
      unless ( defined $this->{'batchInfo'}->{'sequences'}->{$seqID} );

  my $batchIDsRef = $this->{'batchInfo'}->{'sequences'}->{$seqID};

  return ( -1 ) if ( @{$batchIDsRef} > 1 );

  my $batchRef = $this->{'batchInfo'}->{'batches'}->[ $batchIDsRef->[ 0 ] - 1 ];
  foreach my $seq ( @{ $batchRef->{'sequences'} } ) {
    next if ( $seq->{'id'} ne $seqID );
    return ( $position + $seq->{'startPos'} );
  }

  return ( -1 );

}

##-------------------------------------------------------------------------##

=head2 getOverlapBoundaries()

  Use: my  = $obj->getOverlapBoundaries();

  Returns: 

  hashRef => [ 
               'fastaID' => [ boundary1, boundary2, ... ],
               ... 
             ] 

  NOTE: Boundaries are not sorted.


=cut

##-------------------------------------------------------------------------##
sub getOverlapBoundaries {
  my $this = shift;

  my %overlapBoundaries = ();

  return ( \%overlapBoundaries )
      unless ( defined $this->{'batchInfo'}->{'batches'} );

  foreach my $batchRef ( @{ $this->{'batchInfo'}->{'batches'} } ) {
    next if ( $batchRef->{'completeSeqs'} );
    foreach my $seqRef ( @{ $batchRef->{'sequences'} } ) {
      if ( $seqRef->{'startPos'} > 0 ) {
        push @{ $overlapBoundaries{ $seqRef->{'fastaID'} } },
            $seqRef->{'startPos'};
      }
      if ( $seqRef->{'lastPos'} > 0 && !$seqRef->{'lastBatch'} ) {
        push @{ $overlapBoundaries{ $seqRef->{'fastaID'} } },
            $seqRef->{'lastPos'};
      }
    }
  }
  return ( \%overlapBoundaries );

}

sub getSeqIDValidRange {
  my $this  = shift;
  my $seqID = shift;

  my $startPos = -1;
  my $endPos   = -1;

  return ( -1, -1 )
      unless ( defined $this->{'batchInfo'}->{'sequences'}->{$seqID} );

  my $batchIDsRef = $this->{'batchInfo'}->{'sequences'}->{$seqID};

  return ( -1, -1 ) if ( @{$batchIDsRef} > 1 );

  my $batchRef = $this->{'batchInfo'}->{'batches'}->[ $batchIDsRef->[ 0 ] - 1 ];
  my $overlapMiddle = sprintf( "%.0f", $this->{'overlapLength'} / 2 );
  foreach my $seq ( @{ $batchRef->{'sequences'} } ) {
    next if ( $seq->{'id'} ne $seqID );
    if ( $seq->{'startPos'} > 0 ) {
      $startPos = $seq->{'startPos'} + $overlapMiddle + 1;
    }
    if ( $seq->{'lastPos'} > 0 && !$seq->{'lastBatch'} ) {
      $endPos = $seq->{'lastPos'} - $overlapMiddle + 1;
    }
    return ( $startPos, $endPos );
  }

  return ( -1, -1 );

}

##-------------------------------------------------------------------------##

=head2 getBatchAverageGC()

  Use: my $avgGC = $obj->getBatchAverageGC( $batchNum );

=cut

##-------------------------------------------------------------------------##
sub getBatchAverageGC {
  my $this     = shift;
  my $batchNum = shift;

  return ( -1 )
      unless ( defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );

  return $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ]->{'averageGC'};

}

##-------------------------------------------------------------------------##

=head2 isBatchFragmented()

  Use: my $bool = $obj->isBatchFragmented( $batchNum );

=cut

##-------------------------------------------------------------------------##
sub isBatchFragmented {
  my $this     = shift;
  my $batchNum = shift;

  return ( -1 )
      unless ( defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );

  return !$this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ]
      ->{'completeSeqs'};

}

##-------------------------------------------------------------------------##

=head2 getSeqIDFromBatchSeqID()

  Use: my $seqID = $obj->getSeqIDFromBatchSeqID( $batchSeqID );

  Return the original sequence identifier which was used
  in the creation of this batch sequence entry.  The
  batch implementation may assign it's own id for each
  sequences because there may not be a one to one mapping
  between the source database and the batched database.

=cut

##-------------------------------------------------------------------------##
sub getSeqIDFromBatchSeqID {
  my $this       = shift;
  my $batchSeqID = shift;

  return ( "" )
      unless ( defined $this->{'batchInfo'}->{'sequences'}->{$batchSeqID} );

  my $batchIDsRef = $this->{'batchInfo'}->{'sequences'}->{$batchSeqID};

  # A batchSeqID is unique for a set of batches.  So a single
  # ID should only be in one batch.  This is to prevent
  # the caller from using the original seqIDs in this call.
  return ( "" ) if ( @{$batchIDsRef} > 1 );

  my $batchRef = $this->{'batchInfo'}->{'batches'}->[ $batchIDsRef->[ 0 ] - 1 ];
  foreach my $seq ( @{ $batchRef->{'sequences'} } ) {
    next if ( $seq->{'id'} ne $batchSeqID );
    return ( $seq->{'fastaID'} );
  }

  return ( "" );

}

##-------------------------------------------------------------------------##

=head2 getBatchLengths()

  Use: my @batchLengths = $obj->getBatchLengths();

=cut

##-------------------------------------------------------------------------##
sub getBatchLengths {
  my $this = shift;

  my @batchLengths = ();
  return ( @batchLengths ) unless ( defined $this->{'batchInfo'}->{'batches'} );

  return ( map { $_->{'sequenceLength'} }
           @{ $this->{'batchInfo'}->{'batches'} } );

}

##-------------------------------------------------------------------------##

=head2 getBatchSeqLength()

  Use: my $seqLen = $obj->getBatchSeqLength( $batchNum );

=cut

##-------------------------------------------------------------------------##
sub getBatchSeqLength {
  my $this     = shift;
  my $batchNum = shift;

  return ( -1 )
      unless ( defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );
  my $batchRef = $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ];
  return ( $batchRef->{'sequenceLength'} );
}

##-------------------------------------------------------------------------##

=head2 getBatchCount()

  Use: my $count = $obj->getBatchCount();

=cut

##-------------------------------------------------------------------------##
sub getBatchCount {
  my $this = shift;

  return ( 0 ) unless ( defined $this->{'batchInfo'}->{'batches'} );

  return ( $#{ $this->{'batchInfo'}->{'batches'} } + 1 );

}

##-------------------------------------------------------------------------##

=head2 getBatchSeqCount()

  Use: my $count = $obj->getBatchSeqCount( $batchNum );

=cut

##-------------------------------------------------------------------------##
sub getBatchSeqCount {
  my $this     = shift;
  my $batchNum = shift;

  return ( 0 )
      unless ( defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );

  my $batchRef = $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ];

  return ( $#{ $batchRef->{'sequences'} } + 1 );

}

##-------------------------------------------------------------------------##

=head2 getBatchSeqIDs()

  Use: my @seqIDs = $obj->getBatchSeqIDs( $batchNum );

=cut

##-------------------------------------------------------------------------##
sub getBatchSeqIDs {
  my $this     = shift;
  my $batchNum = shift;

  return ( () )
      unless ( defined $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ] );

  my $batch = $this->{'batchInfo'}->{'batches'}->[ $batchNum - 1 ];

  return ( map { $_->{'id'} } @{ $batch->{'sequences'} } );

}

##-------------------------------------------------------------------------##
## Subroutines and Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
##  my        &_packBatches( $obj, $seqDBObj, $fragmentLen, $overlapLen );
##
##      $seqDBObj:        A sequence database object
##      $fragmentLen:     The preferred length of a RepeatMasker batch.
##      $overlapLen:      The length of overlap required (in bp) if a
##                        sequence is fragmented.
##
##  Returns
##
##    Reads in a large 1-N record FASTA file and creates K output
##    files named "$file_batch#.masked" ( where # ranges from 1-K )
##    containing approximatly $fragmentLen amount of sequence data
##    each.
##
##    The files may be bigger or smaller depending on how the
##    parser decides to pack the batches.  In general a large
##    sequence is broken down and stored in it's own set of
##    batches.  A batch is created if we are anywhere between
##    .75 and 1.25 the size of the $fragmentLen or < .75 x
##    the fragment length if the next sequence is > $fragmentLen.
##
##    In addition to the creation of all the files this routine also
##    calculates the total sequence len ( $totSeqLen ) and
##    the GC length ( $totGCLen - total "G" and "C" characters in the
##    sequence ).  It also returns a structure containing details on
##    the batch files created:
##
##  \%batchInfo = { 'batches' =>
##                    [
##                      { 'sequenceLength' => #Total Sequence Size#,
##                        'averageGC'      => #Average GC for this batch#,
##                        'completeSeqs'   => #Are all the sequences complete#,
##                        'sequences' => [ {  'id'     => #batch sequence id#,
##                                            'fastaID' => #The id to the fasta record#,
##                                            'startPos'  => #Optional - Used
##                                                            when sequence
##                                                            is fragmented#,
##                                            'lastPos' => #Optional - Used
##                                                           when sequence is
##                                                           fragmented#
##                                         },
##                                         { 'id' =>  ....
##                                         }
##                                       ]
##                      },
##                      { 'sequenceLength' => ...
##                      }
##                    ]
##                  'sequences' =>  { #SequenceID# => [ batch#1, batch#2..],
##                                    #SequenceID# => [ batch#1, batch#2..],
##                                    ...
##                                  }
##                };
##
##------------------------------------------------------------------------##
sub _packBatches {
  my $this        = shift;
  my $seqDBObj    = shift;
  my $fragmentLen = shift;
  my $overlapLen  = shift;

  # No need to do work unless we have sequences.
  my $seqCount = 0;
  return
      unless ( defined $seqDBObj
               && ( $seqCount = $seqDBObj->getSeqCount() ) > 0 );

  ##
  ## Determine batch boundries
  ##
  my %batchInfo      = ();
  my @batches        = ();
  my $completeSeqs   = 1;
  my @batchSequences = ();
  my $totalGC        = 0;
  my %seqToBatch     = ();
  my $batchLenCtr    = 0;
  my $batchNum       = 1;
  my $seqLength      = 0;
  my $gcLength       = 0;

  foreach my $seqID ( $seqDBObj->getIDs() ) {

    $seqLength = $seqDBObj->getSeqLength( $seqID );
    $gcLength  = $seqDBObj->getGCLength( $seqID );

    if (    $batchLenCtr > 0
         && $batchLenCtr + $seqLength > ( $fragmentLen * 1.25 ) )
    {

      # Finalize batch and move on
      push @batches,
          {
            'sequenceLength' => $batchLenCtr,
            'averageGC'      => int( ( $totalGC / $batchLenCtr ) * 100 ),
            'completeSeqs'   => $completeSeqs,
            'sequences'      => [ @batchSequences ]
          };
      @batchSequences = ();
      $batchNum++;
      $totalGC      = 0;
      $completeSeqs = 1;
      $batchLenCtr  = 0;
    }
    if ( $seqLength > $fragmentLen ) {

      # Should break this one up!
      # About how many pieces can we break this into?
      my $divisor = 2;
      while ( ( $seqLength + ( ( $divisor - 1 ) * $overlapLen ) ) / $divisor >
              $fragmentLen )
      {
        $divisor++;
      }

      # Check to see if we need to finish the last batch
      if ( $batchLenCtr > 0 ) {

        # Finalize batch and move on
        push @batches,
            {
              'sequenceLength' => $batchLenCtr,
              'averageGC'      => int( ( $totalGC / $batchLenCtr ) * 100 ),
              'completeSeqs'   => $completeSeqs,
              'sequences'      => [ @batchSequences ]
            };
        @batchSequences = ();
        $batchLenCtr    = 0;
        $totalGC        = 0;
        $completeSeqs   = 1;
        $batchNum++;
      }

      # Alright...now we can break it up into $size pieces
      my $size =
          int( ( $seqLength + ( ( $divisor - 1 ) * $overlapLen ) ) / $divisor );
      $completeSeqs = 0;

      # Break up
      my $i;
      for ( $i = 0 ; $i < $divisor - 1 ; $i++ ) {

        # Create a batch
        push @batchSequences,
            {
              'id' => $seqID . "frag-" . ( $i + 1 ),
              'fastaID'  => $seqID,
              'startPos' => ( $i * ( $size - $overlapLen ) ),
              'lastPos'  => ( $i * ( $size - $overlapLen ) ) + $size - 1
            };
        push @{ $seqToBatch{$seqID} }, $batchNum;
        push @{ $seqToBatch{ $seqID . "frag-" . ( $i + 1 ) . "" } }, $batchNum;
        my $gcSeq =
            $seqDBObj->getSubstr( $seqID, ( $i * ( $size - $overlapLen ) ),
                                  $size );
        $gcLength = ( $gcSeq =~ tr/GCS// );
        my $gcSize = $size - ( $gcSeq =~ tr/XNRYMK// );
        my $gcLevel = 0;
        $gcLevel = sprintf( "%.0f", ( ( $gcLength / $gcSize ) * 100 ) )
            if ( $gcSize > 0 );
        push @batches,
            {
              'sequenceLength' => $size,
              'averageGC'      => $gcLevel,
              'completeSeqs'   => $completeSeqs,
              'sequences'      => [ @batchSequences ]
            };
        @batchSequences = ();
        $batchNum++;
      }

      # Trailing batch will have remainder sequence attached
      $size =
          ( $size +
            ( ( $seqLength + ( ( $divisor - 1 ) * $overlapLen ) ) % $divisor )
          );
      push @{ $seqToBatch{$seqID} }, $batchNum;
      push @batchSequences,
          {
            'id' => $seqID . "frag-" . ( $i + 1 ) . "",
            'fastaID'   => $seqID,
            'lastBatch' => 1,
            'startPos'  => $seqLength - $size,
            'lastPos'   => $seqLength - 1
          };
      push @{ $seqToBatch{ $seqID . "frag-" . ( $i + 1 ) . "" } }, $batchNum;
      my $gcSeq = $seqDBObj->getSubstr( $seqID, $seqLength - $size, $size );
      $gcLength = ( $gcSeq =~ tr/GCS// );
      my $gcSize = $size - ( $gcSeq =~ tr/XNRYMK// );
      my $gcLevel = 0;
      $gcLevel = sprintf( "%.0f", ( ( $gcLength / $gcSize ) * 100 ) )
          if ( $gcSize > 0 );
      push @batches,
          {
            'sequenceLength' => $size,
            'averageGC'      => $gcLevel,
            'completeSeqs'   => $completeSeqs,
            'sequences'      => [ @batchSequences ]
          };
      @batchSequences = ();
      $completeSeqs   = 1;
      $batchLenCtr    = 0;
      $totalGC        = 0;
      $batchNum++;
    }
    else {
      if (    $seqLength + $batchLenCtr > ( $fragmentLen * 1.25 )
           && $batchLenCtr > 0 )
      {

        # Finish old batch first
        push @batches,
            {
              'sequenceLength' => $batchLenCtr,
              'averageGC'      => int( ( $totalGC / $batchLenCtr ) * 100 ),
              'completeSeqs'   => $completeSeqs,
              'sequences'      => [ @batchSequences ]
            };
        @batchSequences = ();
        $batchLenCtr    = 0;
        $totalGC        = 0;
        $completeSeqs   = 1;
        $batchNum++;
      }

      # Small enough to add on
      push @batchSequences,
          {
            'id'      => $seqID,
            'fastaID' => $seqID
          };
      push @{ $seqToBatch{$seqID} }, $batchNum;
      $totalGC     += $gcLength;
      $batchLenCtr += $seqLength;

      # Mark the data structure
    }
  }
  if ( $batchLenCtr > 0 ) {

    # Finish last batch
    push @batches,
        {
          'sequenceLength' => $batchLenCtr,
          'averageGC'      => int( ( $totalGC / $batchLenCtr ) * 100 ),
          'completeSeqs'   => $completeSeqs,
          'sequences'      => [ @batchSequences ]
        };
  }

  $this->{'batchInfo'}->{'batches'}   = \@batches;
  $this->{'batchInfo'}->{'sequences'} = \%seqToBatch;

}

##-------------------------------------------------------------------------##

=head2 reintegrateBatchFile()

  Use:  my $obj->reintegrateBatchFile(
                                   $batchFileName,
                                   $outputFileName 
                                     );

    $batchFileNameArrayRef:     The fasta files to reintegrate

    **DEPRECATED**

    This routine was originally used to generate the *.mask
    output file.  It reintegrated the individual batch *.mask
    files using the SimpleBatch datastructure. 

=cut

##-------------------------------------------------------------------------##
sub reintegrateBatchFile {
  my $this           = shift;
  my $batchFile      = shift;
  my $outputFileName = shift;

  my $combBatchesSeqDB = FastaDB->new( fileName => $batchFile,
                                       openMode => SeqDBI::ReadOnly );
  my $seqDBObj = $this->{'seqDBObj'};

  open OUT, ">$outputFileName"
      or die "SimpleBatcher::reintegrateBatchFile "
      . "Error could not open file "
      . "$outputFileName!\n";
  my $fastaDataRef = undef;
  my @seqIDs       = $seqDBObj->getIDs();
  foreach my $seqID ( @seqIDs ) {
    my $description = $seqDBObj->getDescription( $seqID );
    my $batchNumRef = $this->{'batchInfo'}->{'sequences'}->{$seqID};
    if ( @{$batchNumRef} == 1 ) {
      my $seq = $combBatchesSeqDB->getSequence( $seqID );

      print OUT ">" . $seqID;
      if ( $description ne "" ) {
        print OUT " " . $description;
      }
      print OUT "\n";
      $seq =~ s/(\S{50})/$1\n/g;
      $seq .= "\n"
          unless ( $seq =~ /.*\n+$/s );
      print OUT $seq;
    }
    else {
      print OUT ">" . $seqID;
      if ( $description ne "" ) {
        print OUT " " . $description;
      }
      print OUT "\n";

      my $seq;
      for ( my $i = 0 ; $i < @{$batchNumRef} ; $i++ ) {
        $seq =
            $combBatchesSeqDB->getSequence(
                                           $seqID . "frag-" . ( $i + 1 ) . "" );

        # TODO: Make overlap length a parameter of object.

        if ( $i == 0 ) {
          $seq =~ s/\w{500}$//;
        }
        elsif ( $i == ( @{$batchNumRef} - 1 ) ) {
          $seq =~ s/^\w{500}//;
        }
        else {
          $seq =~ s/^\w{500}(\w+)\w{500}$/$1/;
        }
        $seq =~ s/(\S{50})/$1\n/g;
        print OUT $seq;

      }
      print OUT "\n" unless ( $seq =~ /.*\n+$/s );
      $seq = "";
    }
  }
  close OUT;

}

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
