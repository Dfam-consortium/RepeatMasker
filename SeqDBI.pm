#!/tools/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SeqDBI.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      A simple sequence database interface. Support for
##      adding/removing and editing of sequences and simple
##      description type annotations.
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
#     Provide a mechanism for providing alphabet distributions
#     This will help with determining if all characters are
#     masked.  Also helps with GC content etc.
#

=head1 NAME

SeqDBI

=head1 SYNOPSIS

use SeqDBI;

Usage:

   $mySeqDBI = $SeqDBI->new();

=head1 DESCRIPTION

A simple sequence database interface. Support for
adding/removing and editing of sequences and simple description 
type annotations.

=head1 SEE ALSO

=over 4

FastaDB, RepeatMasker

=back

=head1 COPYRIGHT

Copyright 2004 Robert Hubley Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence and Package Definition
#
package SeqDBI;
use strict;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

#
# Constants
#
use constant ReadOnly  => 1;
use constant ReadWrite => 2;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw(min max);

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.1;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 new()

  Use:  my $seqDB = SeqDBI->new();   

  Must be implemented by concrete classes.

=cut

##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  croak "Usage: \$seqDB = SeqDBI->new(); " if ( @_ != 0 );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  return $this;
}

##-------------------------------------------------------------------------##
## Public Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 addSequence()

  Use: my $seqID = addSequence( $obj, $name, $description, $sequence );

  Add a sequence to the database.

=cut

##-------------------------------------------------------------------------##
sub addSequence {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 removeSequence()

  Use: my removeSequence( $obj, $seqID );

  Remove a sequence from the database.

=cut

##-------------------------------------------------------------------------##
sub removeSequence {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSeqCount()

  Use: my $seqCount =  &getSeqCount( $obj );

  The number of sequences contained in this sequence database.

=cut

##-------------------------------------------------------------------------##
sub getSeqCount {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSeqLength()

  Use: my $length = &getSeqLength( $obj, [$seqID] );

  The sequence length of the entire database or of a 
  particular entry.

=cut

##-------------------------------------------------------------------------##
sub getSeqLength {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSeqLengths()

  Use: my @sequenceLength = &getSeqLengths( $obj );

  Get an array of sequence length values for all entries
  in the database.

=cut

##-------------------------------------------------------------------------##
sub getSeqLengths {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
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
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getXNLengths()

  Use: my @xnLengths = &getXNLengths( $obj );

  Get the length of each sequence minus 'X' and 'N' characters
  ( case insensitive ).  The values are returned in an array.

=cut

##-------------------------------------------------------------------------##
sub getXNLengths {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
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
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
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
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSubtLength()

  Use my $subtLen = &getSubtLength( $obj, [$seqID] );

  Get length of sequence minus the ambiguous characters.

=cut

##-------------------------------------------------------------------------##
sub getSubtLength {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getIDS()

  Use my @seqIDs = &getIDS( $obj );

  Get an array of all sequence identifiers in the database.

=cut

##-------------------------------------------------------------------------##
sub getIDs {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getDescriptors()

  Use my @descriptors = &getDescriptors( $obj );

  Get an array of all sequence descriptions contained in the database.

=cut

##-------------------------------------------------------------------------##
sub getDescriptors {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getDescription()

  Use: my $desc = &getDescription( $obj, $seqID);

  Get the sequence description text for a given $seqID.

=cut

##-------------------------------------------------------------------------##
sub getDescription {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 compact()

  Use: my &compact();

  Compact database.  Provide a means to remove "holes" or space inefficiencies
  in the database.   

=cut

##-------------------------------------------------------------------------##
sub compact {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSequence()

  Use: my $seqString =  &getSequence( $seqID );

  Get the sequence string for a given seqID.

=cut

##-------------------------------------------------------------------------##
sub getSequence {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

##-------------------------------------------------------------------------##

=head2 getSubstr()

  Use: my $sequence = &getSubstr( $obj, $seqID, $offset, [$length] );

  Get a portion of a sequence from the database.  If the $length parameter
  is left the return value will be everything following the $offset.

=cut

##-------------------------------------------------------------------------##
sub getSubstr {
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
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
  die "Interface SeqDBI does not implement this function\n"
      . "Try instantiating an implementation of this interface\n"
      . "instead!\n";
}

1;
