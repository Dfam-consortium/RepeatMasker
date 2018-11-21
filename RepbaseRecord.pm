#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RepbaseRecord.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
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

RepbaseRecord

=head1 SYNOPSIS

use RepbaseRecord;

Usage:

    myRepbaseRecord = RepbaseRecord->new();

=head1 DESCRIPTION

=head1 SEE ALSO

=over 4

=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut

package RepbaseRecord;
use strict;
##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my ( $self ) = @_;

  # instance-variables:
  my $id;
  my $dataclass;
  my $molecule;
  my @divisions;
  my @accessions;
  my $dateCreated;
  my $dateUpdated;
  my @comments;
  my $RMType;
  my $RMSubType;
  my $RMRefineable;
  my @RMSpecies;
  my @RMSearchStages;
  my @RMBufferStages;
  my $RMDescription;
  my @keywords;
  my $description;
  my $speciesName;
  my $commonSpeciesName;
  my @classification;
  my @references;
  my $xref;
  my $length;
  my %composition;
  my $sequence;
  my $ftLines;

  $self = bless {
                  _id        => $id,
                  _dataclass => $dataclass,
                  _molecule  => $molecule,
                  _divisions => \@divisions,
                  @divisions,
                  _accessions => \@accessions,
                  @accessions,
                  _dateCreated => $dateCreated,
                  _dateUpdated => $dateUpdated,
                  _comments    => \@comments,
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
                  _keywords      => \@keywords,
                  @keywords,
                  _description       => $description,
                  _speciesName       => $speciesName,
                  _commonSpeciesName => $commonSpeciesName,
                  _classification    => \@classification,
                  @classification,
                  _references => \@references,
                  @references,
                  _xref        => $xref,
                  _length      => $length,
                  _composition => \%composition,
                  %composition,
                  _sequence => $sequence,
                  _ftLines  => $ftLines
      },
      ref( $self ) || $self;

  #$self->inherit_from($self->your_base::new());	# adapt when inheriting
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

=head2 getDataclass()

  Use: my $value = getDataclass();

  Get the value of Dataclass.

=cut

##-------------------------------------------------------------------------##
sub getDataclass {
  my ( $self ) = @_;
  $self->{_dataclass};
}

##-------------------------------------------------------------------------##

=head2 getMolecule()

  Use: my $value = getMolecule();

  Get the value of Molecule.

=cut

##-------------------------------------------------------------------------##
sub getMolecule {
  my ( $self ) = @_;
  $self->{_molecule};
}

##-------------------------------------------------------------------------##

=head2 getDateCreated()

  Use: my $value = getDateCreated();

  Get the value of DateCreated.

=cut

##-------------------------------------------------------------------------##
sub getDateCreated {
  my ( $self ) = @_;
  $self->{_dateCreated};
}

##-------------------------------------------------------------------------##

=head2 getDateUpdated()

  Use: my $value = getDateUpdated();

  Get the value of DateUpdated.

=cut

##-------------------------------------------------------------------------##
sub getDateUpdated {
  my ( $self ) = @_;
  $self->{_dateUpdated};
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

=head2 getDescription()

  Use: my $value = getDescription();

  Get the value of Description.

=cut

##-------------------------------------------------------------------------##
sub getDescription {
  my ( $self ) = @_;
  $self->{_description};
}

##-------------------------------------------------------------------------##

=head2 getSpeciesName()

  Use: my $value = getSpeciesName();

  Get the value of SpeciesName.

=cut

##-------------------------------------------------------------------------##
sub getSpeciesName {
  my ( $self ) = @_;
  $self->{_speciesName};
}

##-------------------------------------------------------------------------##

=head2 getCommonSpeciesName()

  Use: my $value = getCommonSpeciesName();

  Get the value of CommonSpeciesName.

=cut

##-------------------------------------------------------------------------##
sub getCommonSpeciesName {
  my ( $self ) = @_;
  $self->{_commonSpeciesName};
}

##-------------------------------------------------------------------------##

=head2 getXref()

  Use: my $value = getXref();

  Get the value of Xref.

=cut

##-------------------------------------------------------------------------##
sub getXref {
  my ( $self ) = @_;
  $self->{_xref};
}

##-------------------------------------------------------------------------##

=head2 getLength()

  Use: my $value = getLength();

  Get the value of Length.

=cut

##-------------------------------------------------------------------------##
sub getLength {
  my ( $self ) = @_;
  $self->{_length};
}

##-------------------------------------------------------------------------##

=head2 getSequence()

  Use: my $value = getSequence();

  Get the value of Sequence.

=cut

##-------------------------------------------------------------------------##
sub getSequence {
  my ( $self ) = @_;
  $self->{_sequence};
}

##-------------------------------------------------------------------------##

=head2 getFTLines()

  Use: my $value = getFTLines();

  Get the value of FT lines.

=cut

##-------------------------------------------------------------------------##
sub getFTLines {
  my ( $self ) = @_;
  $self->{_ftLines};
}

##-------------------------------------------------------------------------##

=head2 getDivisionsElement()

  Use: my $value = getDivisionsElement( $index );

  Get an element from the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub getDivisionsElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getDivisionsRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getAccessionsElement()

  Use: my $value = getAccessionsElement( $index );

  Get an element from the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub getAccessionsElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getAccessionsRef();
  return $$listRef[ $index ];
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

=head2 getKeywordsElement()

  Use: my $value = getKeywordsElement( $index );

  Get an element from the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub getKeywordsElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getKeywordsRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getClassificationElement()

  Use: my $value = getClassificationElement( $index );

  Get an element from the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub getClassificationElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getClassificationRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getReferencesElement()

  Use: my $value = getReferencesElement( $index );

  Get an element from the References Array.

=cut

##-------------------------------------------------------------------------##
sub getReferencesElement {
  my ( $self, $index ) = @_;
  my $listRef = $self->getReferencesRef();
  return $$listRef[ $index ];
}

##-------------------------------------------------------------------------##

=head2 getCompositionHash()

  Use: my %value = getCompositionHash();

  Get a copy of the Composition Hash.

=cut

##-------------------------------------------------------------------------##
sub getCompositionHash {
  my ( $self ) = @_;
  my $hashRef = $self->getCompositionRef();
  return %$hashRef;
}

##-------------------------------------------------------------------------##

=head2 getCompositionElement()

  Use: my $value = getCompositionElement( $key );

  Get element from  the Composition Hash indexed by $key.

=cut

##-------------------------------------------------------------------------##
sub getCompositionElement {
  my ( $self, $key ) = @_;
  my $hashRef = $self->getCompositionRef();
  return $$hashRef{$key};
}

##-------------------------------------------------------------------------##

=head2 getCompositionKeys()

  Use: my $value = getCompositionKeys();

  Get keys for the Composition Hash.

=cut

##-------------------------------------------------------------------------##
sub getCompositionKeys {
  my ( $self ) = @_;
  my %h = $self->getComposition();
  return keys %h;
}

##-------------------------------------------------------------------------##

=head2 getDivisionsArray()

  Use: my @value = getDivisionsArray();

  Get a copy of the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub getDivisionsArray {
  my ( $self ) = @_;
  my $listRef = $self->getDivisionsRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getAccessionsArray()

  Use: my @value = getAccessionsArray();

  Get a copy of the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub getAccessionsArray {
  my ( $self ) = @_;
  my $listRef = $self->getAccessionsRef();
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

=head2 getKeywordsArray()

  Use: my @value = getKeywordsArray();

  Get a copy of the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub getKeywordsArray {
  my ( $self ) = @_;
  my $listRef = $self->getKeywordsRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getClassificationArray()

  Use: my @value = getClassificationArray();

  Get a copy of the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub getClassificationArray {
  my ( $self ) = @_;
  my $listRef = $self->getClassificationRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getReferencesArray()

  Use: my @value = getReferencesArray();

  Get a copy of the References Array.

=cut

##-------------------------------------------------------------------------##
sub getReferencesArray {
  my ( $self ) = @_;
  my $listRef = $self->getReferencesRef();
  return @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getCompositionRef()

  Use: my $value = getCompositionRef();

  Get a reference to the Composition Hash.

=cut

##-------------------------------------------------------------------------##
sub getCompositionRef {
  my ( $self ) = @_;
  my $hashRef = $self->{_composition};
}

##-------------------------------------------------------------------------##

=head2 getDivisionsRef()

  Use: my $value = getDivisionsRef();

  Get a reference to the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub getDivisionsRef {
  my ( $self ) = @_;
  my $listRef = $self->{_divisions};
}

##-------------------------------------------------------------------------##

=head2 getAccessionsRef()

  Use: my $value = getAccessionsRef();

  Get a reference to the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub getAccessionsRef {
  my ( $self ) = @_;
  my $listRef = $self->{_accessions};
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

=head2 getKeywordsRef()

  Use: my $value = getKeywordsRef();

  Get a reference to the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub getKeywordsRef {
  my ( $self ) = @_;
  my $listRef = $self->{_keywords};
}

##-------------------------------------------------------------------------##

=head2 getClassificationRef()

  Use: my $value = getClassificationRef();

  Get a reference to the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub getClassificationRef {
  my ( $self ) = @_;
  my $listRef = $self->{_classification};
}

##-------------------------------------------------------------------------##

=head2 getReferencesRef()

  Use: my $value = getReferencesRef();

  Get a reference to the References Array.

=cut

##-------------------------------------------------------------------------##
sub getReferencesRef {
  my ( $self ) = @_;
  my $listRef = $self->{_references};
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

=head2 clearDataclass()

  Use: clearDataclass();

  Clear the contents of Dataclass.

=cut

##-------------------------------------------------------------------------##
sub clearDataclass {
  my ( $self ) = @_;
  my $v = $self->setDataclass( undef );
}

##-------------------------------------------------------------------------##

=head2 clearMolecule()

  Use: clearMolecule();

  Clear the contents of Molecule.

=cut

##-------------------------------------------------------------------------##
sub clearMolecule {
  my ( $self ) = @_;
  my $v = $self->setMolecule( undef );
}

##-------------------------------------------------------------------------##

=head2 clearDateCreated()

  Use: clearDateCreated();

  Clear the contents of DateCreated.

=cut

##-------------------------------------------------------------------------##
sub clearDateCreated {
  my ( $self ) = @_;
  my $v = $self->setDateCreated( undef );
}

##-------------------------------------------------------------------------##

=head2 clearDateUpdated()

  Use: clearDateUpdated();

  Clear the contents of DateUpdated.

=cut

##-------------------------------------------------------------------------##
sub clearDateUpdated {
  my ( $self ) = @_;
  my $v = $self->setDateUpdated( undef );
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

=head2 clearDescription()

  Use: clearDescription();

  Clear the contents of Description.

=cut

##-------------------------------------------------------------------------##
sub clearDescription {
  my ( $self ) = @_;
  my $v = $self->setDescription( undef );
}

##-------------------------------------------------------------------------##

=head2 clearSpeciesName()

  Use: clearSpeciesName();

  Clear the contents of SpeciesName.

=cut

##-------------------------------------------------------------------------##
sub clearSpeciesName {
  my ( $self ) = @_;
  my $v = $self->setSpeciesName( undef );
}

##-------------------------------------------------------------------------##

=head2 clearCommonSpeciesName()

  Use: clearCommonSpeciesName();

  Clear the contents of CommonSpeciesName.

=cut

##-------------------------------------------------------------------------##
sub clearCommonSpeciesName {
  my ( $self ) = @_;
  my $v = $self->setCommonSpeciesName( undef );
}

##-------------------------------------------------------------------------##

=head2 clearXref()

  Use: clearXref();

  Clear the contents of Xref.

=cut

##-------------------------------------------------------------------------##
sub clearXref {
  my ( $self ) = @_;
  my $v = $self->setXref( undef );
}

##-------------------------------------------------------------------------##

=head2 clearLength()

  Use: clearLength();

  Clear the contents of Length.

=cut

##-------------------------------------------------------------------------##
sub clearLength {
  my ( $self ) = @_;
  my $v = $self->setLength( undef );
}

##-------------------------------------------------------------------------##

=head2 clearSequence()

  Use: clearSequence();

  Clear the contents of Sequence.

=cut

##-------------------------------------------------------------------------##
sub clearSequence {
  my ( $self ) = @_;
  my $v = $self->setSequence( undef );
}

##-------------------------------------------------------------------------##

=head2 clearFTLines()

  Use: clearFTLines();

  Clear the contents of the FT Lines string.

=cut

##-------------------------------------------------------------------------##
sub clearFTLines {
  my ( $self ) = @_;
  my $v = $self->setFTLines( undef );
}

##-------------------------------------------------------------------------##

=head2 clearComposition()

  Use: clearComposition();

  Clear the Composition Hash.

=cut

##-------------------------------------------------------------------------##
sub clearComposition {
  my ( $self ) = @_;
  my $hashRef = $self->getCompositionRef();
  undef %$hashRef;
}

##-------------------------------------------------------------------------##

=head2 deleteCompositionElement()

  Use: deleteCompositionElement( $key );

  Delete the element mapped to the key ( $key ) in the Composition Hash.

=cut

##-------------------------------------------------------------------------##
sub deleteCompositionElement {
  my ( $self, $key ) = @_;
  my $hashRef = $self->getCompositionRef();
  delete $$hashRef{$key};
}

##-------------------------------------------------------------------------##

=head2 clearDivisions()

  Use: clearDivisions();

  Clear the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub clearDivisions {
  my ( $self ) = @_;
  my $listRef = $self->getDivisionsRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearAccessions()

  Use: clearAccessions();

  Clear the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub clearAccessions {
  my ( $self ) = @_;
  my $listRef = $self->getAccessionsRef();
  undef @$listRef;
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

=head2 clearKeywords()

  Use: clearKeywords();

  Clear the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub clearKeywords {
  my ( $self ) = @_;
  my $listRef = $self->getKeywordsRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearClassification()

  Use: clearClassification();

  Clear the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub clearClassification {
  my ( $self ) = @_;
  my $listRef = $self->getClassificationRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 clearReferences()

  Use: clearReferences();

  Clear the References Array.

=cut

##-------------------------------------------------------------------------##
sub clearReferences {
  my ( $self ) = @_;
  my $listRef = $self->getReferencesRef();
  undef @$listRef;
}

##-------------------------------------------------------------------------##

=head2 getDivisionsCount()

  Use: getDivisionsCount();

  Get a count of the elements in the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub getDivisionsCount {
  my ( $self ) = @_;
  my $listRef = $self->getDivisionsRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getAccessionsCount()

  Use: getAccessionsCount();

  Get a count of the elements in the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub getAccessionsCount {
  my ( $self ) = @_;
  my $listRef = $self->getAccessionsRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getCommentsCount()

  Use: getCommentsCount();

  Get a count of the elements in the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub getCommentsCount {
  my ( $self ) = @_;
  my $listRef = $self->getCommentsRef();
  return ( $#$listRef + 1 );
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

=head2 getKeywordsCount()

  Use: getKeywordsCount();

  Get a count of the elements in the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub getKeywordsCount {
  my ( $self ) = @_;
  my $listRef = $self->getKeywordsRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getClassificationCount()

  Use: getClassificationCount();

  Get a count of the elements in the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub getClassificationCount {
  my ( $self ) = @_;
  my $listRef = $self->getClassificationRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 getReferencesCount()

  Use: getReferencesCount();

  Get a count of the elements in the References Array.

=cut

##-------------------------------------------------------------------------##
sub getReferencesCount {
  my ( $self ) = @_;
  my $listRef = $self->getReferencesRef();
  return ( $#$listRef + 1 );
}

##-------------------------------------------------------------------------##

=head2 popDivisions()

  Use: my $value = popDivisions();

  Pop a value off of the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub popDivisions {
  my ( $self ) = @_;
  my $listRef = $self->getDivisionsRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popAccessions()

  Use: my $value = popAccessions();

  Pop a value off of the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub popAccessions {
  my ( $self ) = @_;
  my $listRef = $self->getAccessionsRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popComments()

  Use: my $value = popComments();

  Pop a value off of the Comments Array.

=cut

##-------------------------------------------------------------------------##
sub popComments {
  my ( $self ) = @_;
  my $listRef = $self->getCommentsRef();
  return pop @$listRef;
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

=head2 popKeywords()

  Use: my $value = popKeywords();

  Pop a value off of the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub popKeywords {
  my ( $self ) = @_;
  my $listRef = $self->getKeywordsRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popClassification()

  Use: my $value = popClassification();

  Pop a value off of the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub popClassification {
  my ( $self ) = @_;
  my $listRef = $self->getClassificationRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 popReferences()

  Use: my $value = popReferences();

  Pop a value off of the References Array.

=cut

##-------------------------------------------------------------------------##
sub popReferences {
  my ( $self ) = @_;
  my $listRef = $self->getReferencesRef();
  return pop @$listRef;
}

##-------------------------------------------------------------------------##

=head2 pushDivisions()

  Use: my $value = pushDivisions( $value );

  Push a $value on to the Divisions Array.

=cut

##-------------------------------------------------------------------------##
sub pushDivisions {
  my ( $self, $value ) = @_;
  my $listRef = $self->getDivisionsRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushAccessions()

  Use: my $value = pushAccessions( $value );

  Push a $value on to the Accessions Array.

=cut

##-------------------------------------------------------------------------##
sub pushAccessions {
  my ( $self, $value ) = @_;
  my $listRef = $self->getAccessionsRef();
  push @$listRef, $value;
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

=head2 pushKeywords()

  Use: my $value = pushKeywords( $value );

  Push a $value on to the Keywords Array.

=cut

##-------------------------------------------------------------------------##
sub pushKeywords {
  my ( $self, $value ) = @_;
  my $listRef = $self->getKeywordsRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushClassification()

  Use: my $value = pushClassification( $value );

  Push a $value on to the Classification Array.

=cut

##-------------------------------------------------------------------------##
sub pushClassification {
  my ( $self, $value ) = @_;
  my $listRef = $self->getClassificationRef();
  push @$listRef, $value;
}

##-------------------------------------------------------------------------##

=head2 pushReferences()

  Use: my $value = pushReferences( $value );

  Push a $value on to the References Array.

=cut

##-------------------------------------------------------------------------##
sub pushReferences {
  my ( $self, $value ) = @_;
  my $listRef = $self->getReferencesRef();
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

=head2 setDataclass()

  Use: my $value = setDataclass( $value );

  Set the value of Dataclass.

=cut

##-------------------------------------------------------------------------##
sub setDataclass {
  my ( $self, $value ) = @_;
  $self->{_dataclass} = $value;
}

##-------------------------------------------------------------------------##

=head2 setMolecule()

  Use: my $value = setMolecule( $value );

  Set the value of Molecule.

=cut

##-------------------------------------------------------------------------##
sub setMolecule {
  my ( $self, $value ) = @_;
  $self->{_molecule} = $value;
}

##-------------------------------------------------------------------------##

=head2 setDateCreated()

  Use: my $value = setDateCreated( $value );

  Set the value of DateCreated.

=cut

##-------------------------------------------------------------------------##
sub setDateCreated {
  my ( $self, $value ) = @_;
  $self->{_dateCreated} = $value;
}

##-------------------------------------------------------------------------##

=head2 setDateUpdated()

  Use: my $value = setDateUpdated( $value );

  Set the value of DateUpdated.

=cut

##-------------------------------------------------------------------------##
sub setDateUpdated {
  my ( $self, $value ) = @_;
  $self->{_dateUpdated} = $value;
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

=head2 setDescription()

  Use: my $value = setDescription( $value );

  Set the value of Description.

=cut

##-------------------------------------------------------------------------##
sub setDescription {
  my ( $self, $value ) = @_;
  $self->{_description} = $value;
}

##-------------------------------------------------------------------------##

=head2 setSpeciesName()

  Use: my $value = setSpeciesName( $value );

  Set the value of SpeciesName.

=cut

##-------------------------------------------------------------------------##
sub setSpeciesName {
  my ( $self, $value ) = @_;
  $self->{_speciesName} = $value;
}

##-------------------------------------------------------------------------##

=head2 setCommonSpeciesName()

  Use: my $value = setCommonSpeciesName( $value );

  Set the value of CommonSpeciesName.

=cut

##-------------------------------------------------------------------------##
sub setCommonSpeciesName {
  my ( $self, $value ) = @_;
  $self->{_commonSpeciesName} = $value;
}

##-------------------------------------------------------------------------##

=head2 setXref()

  Use: my $value = setXref( $value );

  Set the value of Xref.

=cut

##-------------------------------------------------------------------------##
sub setXref {
  my ( $self, $value ) = @_;
  $self->{_xref} = $value;
}

##-------------------------------------------------------------------------##

=head2 setLength()

  Use: my $value = setLength( $value );

  Set the value of Length.

=cut

##-------------------------------------------------------------------------##
sub setLength {
  my ( $self, $value ) = @_;
  $self->{_length} = $value;
}

##-------------------------------------------------------------------------##

=head2 setSequence()

  Use: my $value = setSequence( $value );

  Set the value of Sequence.

=cut

##-------------------------------------------------------------------------##
sub setSequence {
  my ( $self, $value ) = @_;
  $self->{_sequence} = $value;
}

##-------------------------------------------------------------------------##

=head2 setFTLines()

  Use: my $value = setFTLines( $value );

  Set the value of FT lines.

=cut

##-------------------------------------------------------------------------##
sub setFTLines {
  my ( $self, $value ) = @_;
  $self->{_ftLines} = $value;
}

##-------------------------------------------------------------------------##

=head2 setComposition()

  Use: my $value = setComposition( $key, $value );

  Set the value of hash element $key

=cut

##-------------------------------------------------------------------------##
sub setComposition {
  my ( $self, $key, $value ) = @_;
  my $hashRef = $self->getCompositionRef();
  $$hashRef{$key} = $value;
}

##-------------------------------------------------------------------------##

=head2 setDivisions()

  Use: my $oldValue = setDivisions( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setDivisions {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getDivisionsRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setAccessions()

  Use: my $oldValue = setAccessions( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setAccessions {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getAccessionsRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setComments()

  Use: my $oldValue = setComments( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setComments {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getCommentsRef();
  $$listRef[ $index ] = $value;
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

=head2 setKeywords()

  Use: my $oldValue = setKeywords( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setKeywords {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getKeywordsRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setClassification()

  Use: my $oldValue = setClassification( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setClassification {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getClassificationRef();
  $$listRef[ $index ] = $value;
}

##-------------------------------------------------------------------------##

=head2 setReferences()

  Use: my $oldValue = setReferences( $index, $value );

  Set the value of Array element $index.

=cut

##-------------------------------------------------------------------------##
sub setReferences {
  my ( $self, $index, $value ) = @_;
  my $listRef = $self->getReferencesRef();
  $$listRef[ $index ] = $value;
}

1;

