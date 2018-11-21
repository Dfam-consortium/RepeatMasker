#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) ArrayListIterator.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of a list iterator.  This iterator acts
##      as an observer to an underlying ArrayList object.  It will
##      adjust it's indexes to accomodate remove and addition of
##      elements to the underlying list.  NOTE: This class is not
##      thread safe.
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
# Implementation Details:
#
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

ArrayListIterator

=head1 SYNOPSIS

use ArrayList;

Usage: 

    my $elementsCollection = ArrayList->new();
    my $elementIterator = $elementsCollection->getIterator();

=head1 DESCRIPTION

  An implementation of a list iterator.  This iterator acts
  as an observer to an underlying ArrayList object.  It will
  adjust it's indexes to accomodate remove and addition of 
  elements to the underlying list.  NOTE: This class is not
  thread safe.

=head1 INSTANCE METHODS

=cut 

package ArrayListIterator;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "ArrayListIterator";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class         = shift;
  my $simpleListRef = shift;
  my $startPosition = shift;

  # Create ourself as a hash
  my $this = {};

  $this->{ $CLASS . '_simpleListObj' }     = $simpleListRef;
  $this->{ $CLASS . '_nextIndex' }         = $startPosition;
  $this->{ $CLASS . '_lastReturnedIndex' } = -1;

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  return $this;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my signalElementRangeRemoved( $indexStart, $count );
##
##	  $indexStart	: The start of the removal
##    $count      : The number elements removed
##
##  Returns
##
##    A private method for use as a callback by the observed ArrayList
##    object.  Used when the ArrayList has shrunk and the indexes
##    need to be checked.
##
##-------------------------------------------------------------------------##
sub signalElementRangeRemoved {
  my $this       = shift;
  my $indexStart = shift;
  my $count      = shift;

  croak $CLASS. "::signalElementRangeRemoved() indexStart missing.\n"
      unless ( defined $indexStart );
  croak $CLASS. "::signalElementRangeRemoved() count missing.\n"
      unless ( defined $indexStart );

  # Adjust current index;
  if ( $this->{ $CLASS . '_nextIndex' } > ( $indexStart + $count ) ) {

    # We are after the removal
    $this->{ $CLASS . '_nextIndex' }         -= $count;
    $this->{ $CLASS . '_lastReturnedIndex' } -= $count;
  }
  elsif ( $this->{ $CLASS . '_nextIndex' } < $indexStart ) {

    # We are before the removal...nothing to do
  }
  else {

    # We are in the middle...set to the start of the removal
    $this->{ $CLASS . '_nextIndex' } = $indexStart;

    # Invalidate the remove operation
    $this->{ $CLASS . '_lastReturnedIndex' } = -1;
  }
}

##-------------------------------------------------------------------------##
## Use: my signalElementRangeAdded( $indexStart, $count );
##
##	  $indexStart	: The start of the insertion
##    $count      : The number elements inserted
##
##  Returns
##
##    A private method for use as a callback by the observed ArrayList
##    object.  Used when the ArrayList has grown and the indexes
##    need to be checked.
##
##-------------------------------------------------------------------------##
sub signalElementRangeAdded {
  my $this       = shift;
  my $indexStart = shift;
  my $count      = shift;

  croak $CLASS. "::signalElementRangeAdded() indexStart missing.\n"
      unless ( defined $indexStart );
  croak $CLASS. "::signalElementRangeAdded() count missing.\n"
      unless ( defined $indexStart );

  # Adjust current index;
  if ( $this->{ $CLASS . '_nextIndex' } >= $indexStart ) {

    # We are after the addition
    $this->{ $CLASS . '_nextIndex' }         += $count;
    $this->{ $CLASS . '_lastReturnedIndex' } += $count;
  }
  else {

    # We are before the addition...nothing to do
  }

}

##-------------------------------------------------------------------------##

=head2 next()

    Use: my $value = $obj->next();

      Obtain the next element in the enumeration.

=cut

##-------------------------------------------------------------------------##
sub next {
  my $this = shift;
  if ( $this->hasNext() ) {
    my $retVal =
        $this->{ $CLASS . '_simpleListObj' }
        ->get( $this->{ $CLASS . '_nextIndex' } );
    $this->{ $CLASS . '_lastReturnedIndex' } = $this->{ $CLASS . '_nextIndex' };
    $this->{ $CLASS . '_nextIndex' }++;
    return $retVal;
  }
  else {
    warn $CLASS . "::next() At the end of the datastructure!\n";
    return;
  }
}

##-------------------------------------------------------------------------##

=head2 hasNext()

    Use: my $boolean = $obj->hasNext();

      Determine if there are any more entities in the enumeration.

=cut

##-------------------------------------------------------------------------##
sub hasNext {
  my $this         = shift;
  my $highestIndex = $this->{ $CLASS . '_simpleListObj' }->size();
  if ( $this->{ $CLASS . '_nextIndex' } >= $highestIndex ) {
    return 0;
  }
  else {
    return 1;
  }
}

##-------------------------------------------------------------------------##

=head2 previous()

    Use: my $value = $obj->previous();

      Obtain the previous element in the enumeration.

=cut

##-------------------------------------------------------------------------##
sub previous {
  my $this = shift;
  if ( $this->hasPrevious() ) {
    $this->{ $CLASS . '_nextIndex' }--;
    my $retVal =
        $this->{ $CLASS . '_simpleListObj' }
        ->get( $this->{ $CLASS . '_nextIndex' } );
    $this->{ $CLASS . '_lastReturnedIndex' } = $this->{ $CLASS . '_nextIndex' };
    return $retVal;
  }
  else {
    warn "$CLASS::previous() At the begining of the datastructure!\n";
    return;
  }
}

##-------------------------------------------------------------------------##

=head2 hasPrevious()

    Use: my $boolean = $obj->hasPrevious();

      Determine if there are any more entities previous to 
      the current position in the enumeration.

=cut

##-------------------------------------------------------------------------##
sub hasPrevious {
  my $this = shift;
  if ( $this->{ $CLASS . '_nextIndex' } == 0 ) {
    return 0;
  }
  else {
    return 1;
  }
}

##-------------------------------------------------------------------------##

=head2 remove()

    Use: $obj->remove();

      Remove the last element retrieved from the enumeration.

=cut

##-------------------------------------------------------------------------##
sub remove {
  my $this = shift;
  if ( $this->{ $CLASS . '_lastReturnedIndex' } >= 0 ) {
    my $tmpIndex = $this->{ $CLASS . '_lastReturnedIndex' };
    $this->{ $CLASS . '_lastReturnedIndex' }--;
    return $this->{ $CLASS . '_simpleListObj' }->remove( $tmpIndex );
  }
  else {
    return ( undef );
  }
}

##-------------------------------------------------------------------------##

=head2 insert()

    Use: $obj->insert( $entity );

     Insert the entity after the last element retrieved from 
     the enumeration.

=cut

##-------------------------------------------------------------------------##
sub insert {
  my $this  = shift;
  my $value = shift;
  return $this->{ $CLASS . '_simpleListObj' }
      ->insert( $value, $this->{ $CLASS . '_nextIndex' } );
}

sub getIterator {
  my $this = shift;
  return $this->{ $CLASS . '_simpleListObj' }
      ->getIterator( $this->{ $CLASS . '_nextIndex' } );
}

sub getIndex {
  my $this = shift;
  return ( $this->{ $CLASS . '_nextIndex' } );
}

sub release {
  my $this = shift;
  $this->{ $CLASS . '_simpleListObj' }->releaseIterator( $this );
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

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
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
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
