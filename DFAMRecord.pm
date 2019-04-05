#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) DFAMRecord.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2011 Developed by
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
#
# To Do:
#
#

=head1 NAME

DFAMRecord

=head1 SYNOPSIS

use DFAMRecord;

Usage:

    myDFAMRecord = DFAMRecord->new();

=head1 DESCRIPTION

=head1 SEE ALSO

=over 4

=back

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut

package DFAMRecord;
use strict;
##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my ( $self ) = @_;

  # instance-variables:
  my $id;
  my @comments;
  my $RMType;
  my $RMSubType;
  my $RMRefineable;
  my @RMSpecies;
  my @RMSearchStages;
  my @RMBufferStages;
  my $RMDescription;
  my @THLines;
  my $recordLines;

  $self = bless {
                  _id       => $id,
                  _comments => \@comments,
                  @comments,
                  _RMType       => $RMType,
                  _RMSubType    => $RMSubType,
                  _RMRefineable => $RMRefineable,
                  _RMSpecies    => \@RMSpecies,
                  @RMSpecies,
                  _RMSearchStages => \@RMSearchStages,
                  @RMSearchStages,
                  _RMBufferStages => \@RMBufferStages,
                  @RMBufferStages,
                  _RMDescription => $RMDescription,
                  _THLines       => \@THLines,
                  _recordLines   => $recordLines
      },
      ref( $self ) || $self;

  return $self;
}

##-------------------------------------------------------------------------##
## Class Methods
##-------------------------------------------------------------------------##

sub specific {
  my ( $self ) = @_;
}

sub inherit_from {
  my ( $self, $base_blessed ) = @_;
  my @l = keys %$base_blessed;
  foreach ( @l ) {
    $self->{$_} = $base_blessed->{$_};
  }
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##

=head2 getId()

  Use: my $value = getId();

  Get the value of Id.

=cut

##-------------------------------------------------------------------------##
sub getId {
  my ( $self ) = @_;
  $self->{_id};
}

##-------------------------------------------------------------------------##

=head2 getCommentsElement()

  Use: my $value = getCommentsElement( $index );

  Get an element from the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub getCommentsElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getCommentsRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getRMType()

  Use: my $value = getRMType();

  Get the value of RMType.

=cut

##-------------------------------------------------------------------------##
sub getRMType {
  my ( $self ) = @_;
  $self->{_RMType};
}

##-------------------------------------------------------------------------##

=head2 getRMSubType()

  Use: my $value = getRMSubType();

  Get the value of RMSubType.

=cut

##-------------------------------------------------------------------------##
sub getRMSubType {
  my ( $self ) = @_;
  $self->{_RMSubType};
}

##-------------------------------------------------------------------------##

=head2 getRMRefineable()

  Use: my $value = getRMRefineable();

  Get the value of RMRefineable

=cut

##-------------------------------------------------------------------------##
sub getRMRefineable {
  my ( $self ) = @_;
  $self->{_RMRefineable};
}

##-------------------------------------------------------------------------##

=head2 getRecordLines()

  Use: my $value = getRecordLines();

  Get the value of RecordLines

=cut

##-------------------------------------------------------------------------##
sub getRecordLines {
  my ( $self ) = @_;
  $self->{_recordLines};
}

##-------------------------------------------------------------------------##

=head2 getRMDescription()

  Use: my $value = getRMDescription();

  Get the value of RMDescription.

=cut

##-------------------------------------------------------------------------##
sub getRMDescription {
  my ( $self ) = @_;
  $self->{_RMDescription};
}

##-------------------------------------------------------------------------##

=head2 getRMSpeciesElement()

  Use: my $value = getRMSpeciesElement( $index );

  Get an element from the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSpeciesElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getRMSearchStagesElement()

  Use: my $value = getRMSearchStagesElement( $index );

  Get an element from the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSearchStagesElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getRMBufferStagesElement()

  Use: my $value = getRMBufferStagesElement( $index );

  Get an element from the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMBufferStagesElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getThreshElement()

  Use: my $value = getThreshElement( $index );

  Get an element from the HMM TH Array.

=cut

##-------------------------------------------------------------------------##
sub getThreshElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getThreshRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getThreshArray()

  Use: my @value = getThreshArray();

  Get a copy of the TH Array.

=cut

##-------------------------------------------------------------------------##
sub getThreshArray {
  my ( $self ) = @_;
  my $listRef = $self->getThreshRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getCommentsArray()

  Use: my @value = getCommentsArray();

  Get a copy of the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub getCommentsArray {
  my ( $self ) = @_;
  my $listRef = $self->getCommentsRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getRMSpeciesArray()

  Use: my @value = getRMSpeciesArray();

  Get a copy of the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSpeciesArray {
  my ( $self ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getRMSearchStagesArray()

  Use: my @value = getRMSearchStagesArray();

  Get a copy of the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSearchStagesArray {
  my ( $self ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getRMBufferStagesArray()

  Use: my @value = getRMBufferStagesArray();

  Get a copy of the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMBufferStagesArray {
  my ( $self ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getRMSpeciesRef()

  Use: my $value = getRMSpeciesRef();

  Get a reference to the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSpeciesRef {
  my ( $self ) = @_;
  my $listRef = $self->{_RMSpecies};
}

##-------------------------------------------------------------------------##

=head2 getThreshRef()

  Use: my $value = getThreshRef();

  Get a reference to the Thresholds Array.

=cut

##-------------------------------------------------------------------------##
sub getThreshRef {
  my ( $self ) = @_;
  my $listRef = $self->{_THLines};
}

##-------------------------------------------------------------------------##

=head2 getCommentsRef()

  Use: my $value = getCommentsRef();

  Get a reference to the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub getCommentsRef {
  my ( $self ) = @_;
  my $listRef = $self->{_comments};
}

##-------------------------------------------------------------------------##

=head2 getRMSearchStagesRef()

  Use: my $value = getRMSearchStagesRef();

  Get a reference to the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSearchStagesRef {
  my ( $self ) = @_;
  my $listRef = $self->{_RMSearchStages};
}

##-------------------------------------------------------------------------##

=head2 getRMBufferStagesRef()

  Use: my $value = getRMBufferStagesRef();

  Get a reference to the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMBufferStagesRef {
  my ( $self ) = @_;
  my $listRef = $self->{_RMBufferStages};
}

##-------------------------------------------------------------------------##

=head2 clearComments()

  Use: clearComments();

  Clear the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub clearComments {
  my ( $self ) = @_;
  my $listRef = $self->getCommentsRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearId()

  Use: clearId();

  Clear the contents of Id.

=cut

##-------------------------------------------------------------------------##
sub clearId {
  my ( $self ) = @_;
  my $v = $self->setId( undef );
}

##-------------------------------------------------------------------------##

=head2 clearRMType()

  Use: clearRMType();

  Clear the contents of RMType.

=cut

##-------------------------------------------------------------------------##
sub clearRMType {
  my ( $self ) = @_;
  my $v = $self->setRMType( undef );
}

##-------------------------------------------------------------------------##

=head2 clearRMSubType()

  Use: clearRMSubType();

  Clear the contents of RMSubType.

=cut

##-------------------------------------------------------------------------##
sub clearRMSubType {
  my ( $self ) = @_;
  my $v = $self->setRMSubType( undef );
}

##-------------------------------------------------------------------------##

=head2 clearRMRefineable()

  Use: clearRMRefineable();

  Clear the contents of RMRefineable.

=cut

##-------------------------------------------------------------------------##
sub clearRMRefineable {
  my ( $self ) = @_;
  my $v = $self->setRMRefineable( undef );
}

##-------------------------------------------------------------------------##

=head2 clearRMDescription()

  Use: clearRMDescription();

  Clear the contents of RMDescription.

=cut

##-------------------------------------------------------------------------##
sub clearRMDescription {
  my ( $self ) = @_;
  my $v = $self->setRMDescription( undef );
}

##-------------------------------------------------------------------------##

=head2 clearRMSpecies()

  Use: clearRMSpecies();

  Clear the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub clearRMSpecies {
  my ( $self ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearRMSearchStages()

  Use: clearRMSearchStages();

  Clear the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub clearRMSearchStages {
  my ( $self ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearThresh()

  Use: clearThresh();

  Clear the TH Array.

=cut

##-------------------------------------------------------------------------##
sub clearThresh {
  my ( $self ) = @_;
  my $listRef = $self->getThreshRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearRMBufferStages()

  Use: clearRMBufferStages();

  Clear the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub clearRMBufferStages {
  my ( $self ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getRMSpeciesCount()

  Use: getRMSpeciesCount();

  Get a count of the elements in the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSpeciesCount {
  my ( $self ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getRMSearchStagesCount()

  Use: getRMSearchStagesCount();

  Get a count of the elements in the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMSearchStagesCount {
  my ( $self ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getThreshCount()

  Use: getThreshCount();

  Get a count of the elements in the TH Array.

=cut

##-------------------------------------------------------------------------##
sub getThreshCount {
  my ( $self ) = @_;
  my $listRef = $self->getThreshRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getRMBufferStagesCount()

  Use: getRMBufferStagesCount();

  Get a count of the elements in the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub getRMBufferStagesCount {
  my ( $self ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 popRMSpecies()

  Use: my $value = popRMSpecies();

  Pop a value off of the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub popRMSpecies {
  my ( $self ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popThresh()

  Use: my $value = popThresh();

  Pop a value off of the TH Array.

=cut

##-------------------------------------------------------------------------##
sub popThresh {
  my ( $self ) = @_;
  my $listRef = $self->getThreshRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popRMSearchStages()

  Use: my $value = popRMSearchStages();

  Pop a value off of the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub popRMSearchStages {
  my ( $self ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popRMBufferStages()

  Use: my $value = popRMBufferStages();

  Pop a value off of the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub popRMBufferStages {
  my ( $self ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 pushComments()

  Use: my $value = pushComments( $value );

  Push a $value on to the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub pushComments {
  my ( $self, $value ) = @_;
  my $listRef = $self->getCommentsRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushRMSpecies()

  Use: my $value = pushRMSpecies( $value );

  Push a $value on to the RMSpecies Array.

=cut

##-------------------------------------------------------------------------##
sub pushRMSpecies {
  my ( $self, $value ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushThresh()

  Use: my $value = pushThresh( $value );

  Push a $value on to the TH Array.

=cut

##-------------------------------------------------------------------------##
sub pushThresh {
  my ( $self, $value ) = @_;
  my $listRef = $self->getThreshRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushRMSearchStages()

  Use: my $value = pushRMSearchStages( $value );

  Push a $value on to the RMSearchStages Array.

=cut

##-------------------------------------------------------------------------##
sub pushRMSearchStages {
  my ( $self, $value ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushRMBufferStages()

  Use: my $value = pushRMBufferStages( $value );

  Push a $value on to the RMBufferStages Array.

=cut

##-------------------------------------------------------------------------##
sub pushRMBufferStages {
  my ( $self, $value ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 setId()

  Use: my $value = setId( $value );

  Set the value of Id.

=cut

##-------------------------------------------------------------------------##
sub setId {
  my ( $self, $value ) = @_;
  $self->{_id} = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMType()

  Use: my $value = setRMType( $value );

  Set the value of RMType.

=cut

##-------------------------------------------------------------------------##
sub setRMType {
  my ( $self, $value ) = @_;
  $self->{_RMType} = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMSubType()

  Use: my $value = setRMSubType( $value );

  Set the value of RMSubType.

=cut

##-------------------------------------------------------------------------##
sub setRMSubType {
  my ( $self, $value ) = @_;
  $self->{_RMSubType} = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMRefineable()

  Use: my $value = setRMRefineable( $value );

  Set the value of RMRefineable.

=cut

##-------------------------------------------------------------------------##
sub setRMRefineable {
  my ( $self, $value ) = @_;
  $self->{_RMRefineable} = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMDescription()

  Use: my $value = setRMDescription( $value );

  Set the value of RMDescription.

=cut

##-------------------------------------------------------------------------##
sub setRMDescription {
  my ( $self, $value ) = @_;
  $self->{_RMDescription} = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMSpecies()

  Use: my $oldValue = setRMSpecies( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setRMSpecies {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getRMSpeciesRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setThresh()

  Use: my $oldValue = setThresh( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setThresh {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getThreshRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMSearchStages()

  Use: my $oldValue = setRMSearchStages( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setRMSearchStages {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getRMSearchStagesRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setRMBufferStages()

  Use: my $oldValue = setRMBufferStages( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setRMBufferStages {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getRMBufferStagesRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setRecordLines()

  Use: my $oldValue = setRecordLines( $value );

  Set the value of recordLines.

=cut

##-------------------------------------------------------------------------##
sub setRecordLines {
  my ( $self, $value ) = @_;
  $self->{_recordLines} = $value;
}

##-------------------------------------------------------------------------##

=head2 get_setAcc()

  Use: my $value    = getAcc();
  Use: my $oldValue = setAcc( $value );

  Get/Set the accession.

=cut

##-------------------------------------------------------------------------##
sub getAcc {
  my $obj = shift;

  my $value = $obj->{'Acc'};

  return $value;
}

sub setAcc {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'Acc'};
  $obj->{'Acc'} = $value;

  return $oldValue;
}

1;

