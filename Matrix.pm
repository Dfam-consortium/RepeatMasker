#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) Matrix.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An for holding a generic biological sequence similarity
##      search result.
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
# bless( {
#        }
#     }, 'Matrix' );
#
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#

=head1 NAME

Matrix

=head1 SYNOPSIS

use Matrix

Usage: 

    my $matrix = Matrix->new();

  or 

    my $matrix = Matrix->new( fileName => foo.matrix );

=head1 DESCRIPTION

A class for storing a RepeatMasker matrix

=head1 SEE ALSO

=over 4

=back

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package Matrix;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

# Matrix Type
use constant crossmatchMatrix => 1;
use constant blastMatrix      => 2;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "Matrix";

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

    # Special parameters first
    if ( defined $nameValuePairs{'fileName'} ) {
      $this->parseFromFile( $nameValuePairs{'fileName'},
                            $nameValuePairs{'matrixType'} );
    }
    else {

      # Default calls to set routines
      while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
        my $method = "set" . _ucFirst( $name );
        unless ( $this->can( $method ) ) {
          croak( "Matrix::add: Instance variable $name doesn't exist." . "" );
        }
        $this->$method( $value );
      }
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a Matrix *duplicating* all the values of the old
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

=head2 get_setMatrixName()

  Use: my $value    = getMatrixName( );
  Use: my $oldValue = setMatrixName( $value );

  Get/Set the name of the matrix.

=cut

##-------------------------------------------------------------------------##
sub getMatrixName {
  my $obj = shift;

  my $value = $obj->{'matrixName'};

  return $value;
}

sub setMatrixName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'matrixName'};
  $obj->{'matrixName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAlphabetRef()

  Use: my $alphabetArrayRef = getAlphabetRef();
  Use: my $oldAlphabetArrayRef = setAlphabetRef( $alphabetArrayRef );

  Get/Set the alphabet of the matrix.

=cut

##-------------------------------------------------------------------------##
sub getAlphabetRef {
  my $obj = shift;

  my $value = $obj->{'alphabetArrayRef'};

  return $value;
}

sub setAlphabetRef {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'alphabetArrayRef'};
  $obj->{'alphabetArrayRef'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 getMatrixValuesRef()

  Use: my $alphabetArrayRef = getMatrixValuesRef();

  Get the matrix.

=cut

##-------------------------------------------------------------------------##
sub getMatrixValuesRef {
  my $obj = shift;

  my $value = $obj->{'matrixValues'};

  return $value;
}

##-------------------------------------------------------------------------##

=head2 getMatrixFreqsRef()

  Use: my $alphabetArrayRef = getMatrixFreqsRef();

  Get the matrix.

=cut

##-------------------------------------------------------------------------##
sub getMatrixFreqsRef {
  my $obj = shift;

  my $value = $obj->{'alphabetFreqArrayRef'};

  return $value;
}

##-------------------------------------------------------------------------##

=head2 getLambda()

  Use: my $lambda = getLambda();

  Get the 

=cut

##-------------------------------------------------------------------------##
sub getLambda {
  my $obj = shift;

  my $value = $obj->{'lambda'};

  return $value;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2

  Use: toString();

  Create an string representation of the matrix.

=cut

##-------------------------------------------------------------------------##
sub toString {
  my $obj = shift;
}

##-------------------------------------------------------------------------##

=head2

  Use: parseFromFile( $csvString );

  Populate object with values stored in CSV format.

=cut

##-------------------------------------------------------------------------##
sub parseFromFile {
  my $this       = shift;
  my $fileName   = shift;
  my $matrixType = shift;

  $matrixType = Matrix::blastMatrix if ( !defined $matrixType );

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  open MATRIX, "<$fileName"
      or croak( "$subroutine: Could not open $fileName for reading!" );

  my $row = 0;
  while ( <MATRIX> ) {
    $_ = uc( $_ );
    chomp;

    # FREQS A 0.325 C 0.175 G 0.175 T 0.325
    if ( /FREQS\s+(.*)/ ) {
      $this->{'freqHashRef'} = { split( " ", $1 ) };
    }

    #
    #   A  C  G  T
    # A 1 -1 -1 -1
    if ( /^\s*[A-Z]\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+/ ) {
      s/ //g;
      $this->{'alphabetArrayRef'} = [ split( //, $_ ) ];
    }
    elsif ( $this->{'alphabetArrayRef'}
            && /^\s*([^\d]\s+)?[\d-]+\s+[\d-]+\s+[\d-]+\s+[\d-]+/ )
    {
      my @rowValues = split;
      if ( $1 ) {
        shift @rowValues;
      }
      $this->{'matrixValues'}->[ $row++ ] = [ @rowValues ];
    }

  }
  close MATRIX;

  if ( !defined $this->{'freqHashRef'} ) {
    $this->{'freqHashRef'}->{"A"} = 0.25;
    $this->{'freqHashRef'}->{"C"} = 0.25;
    $this->{'freqHashRef'}->{"G"} = 0.25;
    $this->{'freqHashRef'}->{"T"} = 0.25;
  }

  for ( my $i = 0 ; $i <= $#{ $this->{'alphabetArrayRef'} } ; $i++ ) {
    $this->{'alphabetFreqArrayRef'}->[ $i ] =
        $this->{'freqHashRef'}->{ $this->{'alphabetArrayRef'}->[ $i ] } || 0;
  }

  $this->_calculateLambda();

}

##-------------------------------------------------------------------------##
## Use: my transposeMatrix();
##
## Returns
##      Transposes the matrix values ( ie. rows -> cols ) inside
##      the object and recalculates the lambda value.
##-------------------------------------------------------------------------##
sub transposeMatrix {
  my $this    = shift;
  my @tMatrix = ();
  for ( my $i = 0 ; $i <= $#{ $this->{'alphabetArrayRef'} } ; $i++ ) {
    for ( my $j = 0 ; $j <= $#{ $this->{'alphabetArrayRef'} } ; $j++ ) {
      push @{ $tMatrix[ $i ] }, $this->{'matrixValues'}->[ $j ]->[ $i ];
    }
  }
  $this->{'matrixValues'} = \@tMatrix;
  $this->_calculateLambda();
}

##-------------------------------------------------------------------------##
## Use: my _calculateLambda( $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $lambda : The lambda parameter derived from the matrix
##                and the matrix alphabet frequencies.  This
##                is derived from Phil Green's swat/cross_match
##                programs.
##-------------------------------------------------------------------------##
sub _calculateLambda {
  my $this = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $lambda_upper = 0;
  my $lambda_lower = 0;
  my $lambda       = 0.5;

  my $matFreqsRef  = $this->{'alphabetFreqArrayRef'};
  my $matScoresRef = $this->{'matrixValues'};
  my $sum          = 0;

  do {
    $sum = 0; # 9/8/2021 Fixed a major bug in lambda calculation
    # NOTE: This fix only had an impact on a single RM matrix
    #       ( identity.matrix ) as this is the only RM matrix
    #       that exceeded a lambda of 1.0.
    my $check = 0;
    for ( my $i = 0 ; $i <= $#$matFreqsRef ; $i++ ) {
      for ( my $j = 0 ; $j <= $#$matFreqsRef ; $j++ ) {
        if ( $$matFreqsRef[ $i ] && $$matFreqsRef[ $j ] ) {
          $sum += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ] *
              exp( $lambda * $$matScoresRef[ $i ][ $j ] );
          $check += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ];
        }
      }
    }
    croak "$subroutine: Sanity check failed lowerLambda freqSum = $check"
        . " but should have been between close to 1 ( .999 - 1.001 )"
        if (    $check > 1.001
             || $check < .999 );

    if ( $sum < 1.0 ) {
      $lambda_lower = $lambda;
      $lambda *= 2.0;
    }
  } while ( $sum < 1.0 );
  $lambda_upper = $lambda;

  while ( $lambda_upper - $lambda_lower > .00001 ) {
    $lambda = ( $lambda_lower + $lambda_upper ) / 2.0;
    my $sum   = 0;
    my $check = 0;
    for ( my $i = 0 ; $i <= $#$matFreqsRef ; $i++ ) {
      for ( my $j = 0 ; $j <= $#$matFreqsRef ; $j++ ) {
        if ( $$matFreqsRef[ $i ] && $$matFreqsRef[ $j ] ) {
          $sum += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ] *
              exp( $lambda * $$matScoresRef[ $i ][ $j ] );
          $check += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ];
        }
      }
    }
    croak "$subroutine: Sanity check failed upperLambda freqSum = $check"
        . " but should have been between close to 1 ( .999 - 1.001 )"
        if (    $check > 1.001
             || $check < .999 );

    if ( $sum >= 1.0 ) {
      $lambda_upper = $lambda;
    }
    else {
      $lambda_lower = $lambda;
    }
  }

  $this->{'lambda'} = $lambda;
}

1;
