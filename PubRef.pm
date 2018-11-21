#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) PubRef.pm
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

PubRef

=head1 SYNOPSIS

use PubRef;

Usage:

    myPubRef = PubRef->new();

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

package PubRef;
use strict;
##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my ( $self ) = @_;

  # instance-variables:
  my $number;
  my $comment;
  my $xref;
  my $authors;
  my $title;
  my $location;
  my $position;

  $self = bless {
                  _number   => $number,
                  _comment  => $comment,
                  _xref     => $xref,
                  _authors  => $authors,
                  _title    => $title,
                  _location => $location,
                  _position => $position,
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

=head2 getNumber()

  Use: my $value = getNumber();

  Get the value of Number.

=cut

##-------------------------------------------------------------------------##
sub getNumber {
  my ( $self ) = @_;
  $self->{_number};
}

##-------------------------------------------------------------------------##

=head2 getComment()

  Use: my $value = getComment();

  Get the value of Comment.

=cut

##-------------------------------------------------------------------------##
sub getComment {
  my ( $self ) = @_;
  $self->{_comment};
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

=head2 getAuthors()

  Use: my $value = getAuthors();

  Get the value of Authors.

=cut

##-------------------------------------------------------------------------##
sub getAuthors {
  my ( $self ) = @_;
  $self->{_authors};
}

##-------------------------------------------------------------------------##

=head2 getTitle()

  Use: my $value = getTitle();

  Get the value of Title.

=cut

##-------------------------------------------------------------------------##
sub getTitle {
  my ( $self ) = @_;
  $self->{_title};
}

##-------------------------------------------------------------------------##

=head2 getLocation()

  Use: my $value = getLocation();

  Get the value of Location.

=cut

##-------------------------------------------------------------------------##
sub getLocation {
  my ( $self ) = @_;
  $self->{_location};
}

##-------------------------------------------------------------------------##

=head2 getPosition()

  Use: my $value = getPosition();

  Get the value of Position.

=cut

##-------------------------------------------------------------------------##
sub getPosition {
  my ( $self ) = @_;
  $self->{_position};
}

##-------------------------------------------------------------------------##

=head2 clearNumber()

  Use: clearNumber();

  Clear the contents of Number.

=cut

##-------------------------------------------------------------------------##
sub clearNumber {
  my ( $self ) = @_;
  my $v = $self->setNumber( undef );
}

##-------------------------------------------------------------------------##

=head2 clearComment()

  Use: clearComment();

  Clear the contents of Comment.

=cut

##-------------------------------------------------------------------------##
sub clearComment {
  my ( $self ) = @_;
  my $v = $self->setComment( undef );
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

=head2 clearAuthors()

  Use: clearAuthors();

  Clear the contents of Authors.

=cut

##-------------------------------------------------------------------------##
sub clearAuthors {
  my ( $self ) = @_;
  my $v = $self->setAuthors( undef );
}

##-------------------------------------------------------------------------##

=head2 clearTitle()

  Use: clearTitle();

  Clear the contents of Title.

=cut

##-------------------------------------------------------------------------##
sub clearTitle {
  my ( $self ) = @_;
  my $v = $self->setTitle( undef );
}

##-------------------------------------------------------------------------##

=head2 clearLocation()

  Use: clearLocation();

  Clear the contents of Location.

=cut

##-------------------------------------------------------------------------##
sub clearLocation {
  my ( $self ) = @_;
  my $v = $self->setLocation( undef );
}

##-------------------------------------------------------------------------##

=head2 clearPosition()

  Use: clearPosition();

  Clear the contents of Position.

=cut

##-------------------------------------------------------------------------##
sub clearPosition {
  my ( $self ) = @_;
  my $v = $self->setPosition( undef );
}

##-------------------------------------------------------------------------##

=head2 setNumber()

  Use: my $value = setNumber( $value );

  Set the value of Number.

=cut

##-------------------------------------------------------------------------##
sub setNumber {
  my ( $self, $value ) = @_;
  $self->{_number} = $value;
}

##-------------------------------------------------------------------------##

=head2 setComment()

  Use: my $value = setComment( $value );

  Set the value of Comment.

=cut

##-------------------------------------------------------------------------##
sub setComment {
  my ( $self, $value ) = @_;
  $self->{_comment} = $value;
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

=head2 setAuthors()

  Use: my $value = setAuthors( $value );

  Set the value of Authors.

=cut

##-------------------------------------------------------------------------##
sub setAuthors {
  my ( $self, $value ) = @_;
  $self->{_authors} = $value;
}

##-------------------------------------------------------------------------##

=head2 setTitle()

  Use: my $value = setTitle( $value );

  Set the value of Title.

=cut

##-------------------------------------------------------------------------##
sub setTitle {
  my ( $self, $value ) = @_;
  $self->{_title} = $value;
}

##-------------------------------------------------------------------------##

=head2 setLocation()

  Use: my $value = setLocation( $value );

  Set the value of Location.

=cut

##-------------------------------------------------------------------------##
sub setLocation {
  my ( $self, $value ) = @_;
  $self->{_location} = $value;
}

##-------------------------------------------------------------------------##

=head2 setPosition()

  Use: my $value = setPosition( $value );

  Set the value of Position.

=cut

##-------------------------------------------------------------------------##
sub setPosition {
  my ( $self, $value ) = @_;
  $self->{_position} = $value;
}

1;

