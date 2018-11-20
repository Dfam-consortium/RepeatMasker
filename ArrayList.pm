#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) ArrayList.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of a generic list collection utilizing
##      a perl array to hold data.
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
# bless( { 'ArrayList_elements' => [ $element1,
#                          $element2,
#                         ..
#                       ]
#        }
#     }, 'ArrayList' );
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

ArrayList

=head1 SYNOPSIS

use ArrayList

Usage: 

    $elementsCollection = ArrayList->new();

=head1 DESCRIPTION

A class for storing list type data.

=head1 INSTANCE METHODS

=cut 

package ArrayList;
use strict;
use Data::Dumper;
use ArrayListIterator;
use Scalar::Util qw(weaken isweak blessed);
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "ArrayList";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class    = shift;
  my $arrayRef = shift;

  # Create ourself as a hash
  my $this = {};

  $this->{ $CLASS . '_elements' } = [];

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  if (    defined $arrayRef
       && ref( $arrayRef ) eq "ARRAY"
       && $#{$arrayRef} >= 0 )
  {
    $this->addAll( $arrayRef );
  }

  return $this;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 add()

  Use: $obj->add( $value );

  Add a element to the collection.

=cut

##-------------------------------------------------------------------------##
sub add {
  my $this  = shift;
  my $value = shift;

  croak $CLASS. "::add() value missing.\n" unless ( defined $value );

  push @{ $this->{ $CLASS . '_elements' } }, $value;

}

##-------------------------------------------------------------------------##

=head2 addAll()

  Use: $obj->addAll( $ArrayListRef );
   or
  Use: $obj->addAll( \@arrayRef );

  Append a set of values from a primitive array or an 
  entire ArrayList collection to this one. 

=cut

##-------------------------------------------------------------------------##
sub addAll {
  my $this = shift;

  croak $CLASS. "::addAll() array/arrayList missing.\n"
      unless ( @_ );

  croak $CLASS
      . "::addAll() value is not a primitive "
      . "perl array or an ArrayList.  Cannot append "
      . "to the collection "
      if (
           $#_ > 0
           || (    ( !ref( $_[ 0 ] ) eq "ARRAY" )
                && ( !ref( $_[ 0 ] ) eq "ArrayList " ) )
      );

  if ( blessed( $_[ 0 ] ) && $_[ 0 ]->isa( "ArrayList" ) ) {
    for ( my $i = 0 ; $i < $_[ 0 ]->size() ; $i++ ) {
      $this->add( $_[ 0 ]->get( $i ) );
    }
  }
  else {
    for ( my $i = 0 ; $i <= $#{ $_[ 0 ] } ; $i++ ) {
      $this->add( $_[ 0 ]->[ $i ] );
    }
  }

}

##-------------------------------------------------------------------------##

=head2 insert()

  Use: $obj->insert( $value, $index );

  Insert value into list at index.  This shift all elements
  starting at index down by one.

=cut

##-------------------------------------------------------------------------##
sub insert {
  my $obj   = shift;
  my $value = shift;
  my $index = shift;

  croak $CLASS. "::insert() value missing.\n"
      unless ( defined $value );
  croak $CLASS. "::insert() index missing.\n"
      unless ( defined $index );

  splice( @{ $obj->{ $CLASS . '_elements' } }, $index, 0, ( $value ) );

  # Notify our iterators
  if ( defined $obj->{ $CLASS . '_iterators' } ) {
    for ( my $i = $#{ $obj->{ $CLASS . '_iterators' } } ; $i >= 0 ; $i-- ) {
      my $it = $obj->{ $CLASS . '_iterators' }->[ $i ];
      if ( ref( $it ) ne "REF" ) {
        splice( @{ $obj->{ $CLASS . '_iterators' } }, $i, 1 );
      }
      else {
        $$it->signalElementRangeAdded( $index, 1 );
      }
    }
  }

  return $value;
}

##-------------------------------------------------------------------------##

=head2 remove()

  Use: my $oldValueHash = $obj->remove( $elementIndex );

  Remove a element from the collection.

=cut

##-------------------------------------------------------------------------##
sub remove {
  my $obj          = shift;
  my $elementIndex = shift;

  croak $CLASS. "::remove: Index out of bounds ( $elementIndex )."
      if (    $elementIndex < 0
           || $elementIndex > $#{ $obj->{ $CLASS . '_elements' } } );

  # Notify our iterators
  if ( defined $obj->{ $CLASS . '_iterators' } ) {
    for ( my $i = $#{ $obj->{ $CLASS . '_iterators' } } ; $i >= 0 ; $i-- ) {
      my $it = $obj->{ $CLASS . '_iterators' }->[ $i ];
      if ( ref( $it ) ne "REF" ) {
        splice( @{ $obj->{ $CLASS . '_iterators' } }, $i, 1 );
      }
      else {
        $$it->signalElementRangeRemoved( $elementIndex, 1 );
      }
    }
  }

  return splice( @{ $obj->{ $CLASS . '_elements' } }, $elementIndex, 1 );

}

##-------------------------------------------------------------------------##

=head2 set()

  Use:  my  $oldResultRef = $obj->set( $index, $value );

  Sets the element at the specified position in this collection
  and returns the element previously stored at this position.

=cut

##-------------------------------------------------------------------------##
sub set {
  my $this  = shift;
  my $index = shift;
  my $value = shift;

  croak $CLASS. "::set() index missing.\n" unless ( defined $index );

  croak $CLASS. "::set( $index, ) value missing.\n"
      unless ( defined $value );

  croak $CLASS. "set( $index ) Index out of bounds!\n"
      if ( $index < 0 || $index >= $this->size() );

  my $retVal = $this->{ $CLASS . '_elements' }->[ $index ];
  $this->{ $CLASS . '_elements' }->[ $index ] = $value;

  return ( $retVal );

}

##-------------------------------------------------------------------------##

=head2 get()

  Use:  my $value = $obj->get( $index );

  Returns the element at the specified position in this collection.

=cut

##-------------------------------------------------------------------------##
sub get {
  my $this  = shift;
  my $index = shift;

  croak $CLASS. "::get( $index ) Index out of bounds!\n"
      if ( $index < 0 || $index >= $this->size() );

  return ( $this->{ $CLASS . '_elements' }->[ $index ] );

}

##-------------------------------------------------------------------------##

=head2 clear()

  Use: $obj->clear();

  Removes all of the elements from this collection.

=cut

##-------------------------------------------------------------------------##
sub clear {
  my $this = shift;

  undef $this->{ $CLASS . '_elements' };
  $this->{ $CLASS . '_elements' } = [];

}

##-------------------------------------------------------------------------##

=head2 size()

  Use:  my $numRecords = $obj->size();

  The number of records stored in this object.

=cut

##-------------------------------------------------------------------------##
sub size {
  my $this = shift;

  if ( defined $this->{ $CLASS . '_elements' } ) {
    return ( $#{ $this->{ $CLASS . '_elements' } } + 1 );
  }
  else {
    return 0;
  }

}

##-------------------------------------------------------------------------##

=head2 sort()

  Use: $obj->sort( $comparatorRef );

  Apply a one-time sort to this collection using the supplied
  compartor.  NOTE: The comparator code reference is assumed
  to be a typical PERL sort code block.  This means that
  the two elements being compared will already be assigned
  to the $a and $b parameters inside the comparator.  The
  "<=>" and "cmp" operators can be used to rank numerically
  and alphabetically.

  ie.

  $obj->sort( sub ($$) { 
                          $_[0] <=> $_[1] 
                       } );

=cut

##-------------------------------------------------------------------------##
sub sort {
  my $this          = shift;
  my $comparatorRef = shift;

  croak $CLASS. "::sort() comparatorRef missing.\n"
      unless ( defined $comparatorRef );

  croak $CLASS
      . "::sort( "
      . ref( $comparatorRef )
      . " ) comparatorRef is "
      . "the wrong type. It should be a procedure CODE.\n"
      unless ( ref( $comparatorRef ) eq "CODE" );

  if ( $this->size() ) {
    my @sortedResults =
        sort $comparatorRef @{ $this->{ $CLASS . '_elements' } };
    if ( $#sortedResults == $this->size() - 1 ) {
      $this->{ $CLASS . '_elements' } = \@sortedResults;
    }
    else {
      croak $CLASS. "::sort() error occurred in sort.\n";
    }
  }
}

##-------------------------------------------------------------------------##

=head2 toArrayRef()

  Use: $obj->toArrayRef();

  A helper function which returns an ArrayList object as
  a native perl array.

=cut

##-------------------------------------------------------------------------##
sub toArrayRef {
  my $this = shift;
  return $this->{ $CLASS . '_elements' };
}

##-------------------------------------------------------------------------##

=head2 getIterator()

  Use: $obj->getIterator( $index );

  Returns an interator over the ArrayList.  If $index is supplied
  the iterator will begin at element $index+1 otherwise it will 
  begin at the first element in the list.

=cut

##-------------------------------------------------------------------------##
sub getIterator {
  my $this  = shift;
  my $index = shift;

  my $it;
  if ( defined $index ) {
    if ( $index >= 0 && $index <= $this->size() ) {
      $it = ArrayListIterator->new( $this, $index );
    }
    else {
      croak $CLASS. "::getIterator() index out of bounds!\n";
    }
  }
  else {
    $it = ArrayListIterator->new( $this, 0 );
  }

  # Free up some memory if we have an excessive amount of
  # iterator slots.
  if ( $#{ $this->{ $CLASS . '_iterators' } } > 50 ) {
    for ( my $i = $#{ $this->{ $CLASS . '_iterators' } } ; $i >= 0 ; $i-- ) {
      my $it = $this->{ $CLASS . '_iterators' }->[ $i ];
      if ( ref( $it ) ne "REF" ) {
        splice( @{ $this->{ $CLASS . '_iterators' } }, $i, 1 );
      }
    }
  }
  my $weakRef = $it;

  # Weaken the reference so we don't prohibit the iterator from
  # being eaten by the garbage collector.
  &weaken( $weakRef );
  push @{ $this->{ $CLASS . '_iterators' } }, \$weakRef;
  return $it;
}

## DEBUG ROUTINE...remove
sub numIters {
  my $this = shift;
  print "Num iterators = " . $#{ $this->{ $CLASS . '_iterators' } } . "\n";
  foreach my $it ( @{ $this->{ $CLASS . '_iterators' } } ) {
    print "ref=" . ref( $it ) . "\n";
  }
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
