#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) Taxonomy.pm
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
# bless(
#      'Taxonomy' );
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

Taxononmy

=head1 SYNOPSIS

use Taxononmy;

Usage: 

# Build initial taxonomy cached database
my $taxDB = Taxonomy->new( taxonomyDataFile=>"taxonomy.dat",
                           namesDmpFile=>"NCBITaxonomyDB/names.dmp",
                           nodesDmpFile=>"NCBITaxonomyDB/nodes.dmp",
                           rmDatabaseFile=>"Libraries/ReapeatMasker.lib" );

if ( $taxDB->isA( "Mouse", "Mammalia" ) ) {
  print "A Mouse is a Mammal!\n";
}else {
  print "Error...a Mouse should be Mammal!\n";
}

# All future executions may use this cached database to quickly
# load
my $taxDB = Taxonomy->new( taxonomyDataFile=>"taxonomy.dat" );

=head1 DESCRIPTION

A general object for querying the NCBI Taxonomy database.  The
database may be in the original dump format or in an application
specific "frozen" format.  

The NCBI Taxonomy database is filtered for use solely with the
RepeatMasker library.  Names in the taxonomy database are removed
if they are not related ( child or ancestor ) to a species contained 
in the RepeatMasker.lib database.

=head1 INSTANCE METHODS

=cut 

package Taxonomy;
use strict;
use FindBin;
use Data::Dumper;
use FastaDB;
use EMBL;
use Text::Soundex;
use Carp;

use Storable qw(nstore retrieve);
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.1;
my $CLASS   = "Taxononmy";

##-------------------------------------------------------------------------##

=begin

=over 4

=item my $instance = Taxonomy->new( taxonomyDataFile=>"filename",
                                    [namesDmpFile=>"filename",
                                     nodesDmpFile=>"filename",
                                     rmDatabaseFile=>"filename"] );

Construct a new Taxonomy object.  Use the taxononmy data file 
specified to populate the object or build a new taxonomy data
file using the NCBI/RepeatMasker data files:

  namesDmpFile:  NCBI Taxonomy DB names.dmp file
  nodesDmpFile:  NCBI Taxonomy DB nodes.dmp file
  rmDatabaseFile: RepeatMasker monolithic repeat database file

=back

=cut

##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  my $this = {};

  my %supplementalSynonyms = (
                               'anopheles'   => "anopheles genus",
                               'carnivore'   => "carnivora",
                               'artiodactyl' => "cetartiodactyla",
                               'puffer'      => "fugu",
                               'fruitfly'    => "drosophila fruit fly genus",
                               'drosophila'  => "drosophila fruit fly genus",
                               'worm'        => "caenorhabditis",
                               'elegans'     => "caenorhabditis",
                               'mustard'     => "arabidopsis",
                               'corn'        => "zea",
                               'maize'       => "zea",
                               'theria'      => "theria mammalia",
                               "mus"         => "mouse",
                               'cionaint'    => "ciona intestinalis",
                               'cionasav'    => "ciona savignyi",
                               'diatom'      => "diatomea",
                               'rodent'      => "rodentia",
                               'rat'         => "rattus",
                               'mammal'      => "mammalia"
  );

  if (    defined $nameValuePairs{'namesDmpFile'}
       && defined $nameValuePairs{'nodesDmpFile'} )
  {
    if (    -s $nameValuePairs{'namesDmpFile'}
         && -s $nameValuePairs{'nodesDmpFile'} )
    {

      my $dataFile = "taxonomy.dat";
      if ( defined $nameValuePairs{'taxonomyDataFile'} ) {
        $dataFile = $nameValuePairs{'taxonomyDataFile'};
      }

      # Bless this hash in the name of the father, the son...
      bless $this, $class;

      if ( defined $nameValuePairs{'rmDatabaseFile'}
           && -s $nameValuePairs{'rmDatabaseFile'} )
      {
        $this->_buildFromNCBIDB(
                                 $nameValuePairs{'namesDmpFile'},
                                 $nameValuePairs{'nodesDmpFile'},
                                 $dataFile,
                                 \%supplementalSynonyms,
                                 $nameValuePairs{'rmDatabaseFile'}
        );
      }
      else {
        $this->_buildFromNCBIDB(
                                 $nameValuePairs{'namesDmpFile'},
                                 $nameValuePairs{'nodesDmpFile'},
                                 $dataFile,
                                 \%supplementalSynonyms
        );
      }
    }
    else {
      croak $CLASS. "::new() NCBI Taxonomy DB dump files do not exist!\n";
    }

  }
  elsif (    defined $nameValuePairs{'taxonomyArray'}
          && defined $nameValuePairs{'synonymHash'} )
  {

    #print join(", ", @INC ) . "\n";
    print "Storable version: $Storable::VERSION\n";

    # Bless this hash in the name of the father, the son...
    bless $this, $class;

    $this->_buildFromStructures(
                                 $nameValuePairs{'taxonomyArray'},
                                 $nameValuePairs{'synonymHash'},
                                 "taxonomy.dat",
                                 \%supplementalSynonyms
    );

  }
  elsif ( defined $nameValuePairs{'taxonomyDataFile'}
          && -s $nameValuePairs{'taxonomyDataFile'} )
  {

    # Read in a serialized version of ourselves
    $this = retrieve( $nameValuePairs{'taxonomyDataFile'} );

    # Bless this hash in the name of the father, the son...
    bless $this, $class;

  }
  elsif ( -s "taxonomy.dat" ) {

    # Read in a serialized version of ourselves
    $this = retrieve( "taxonomy.dat" );

    # Bless this hash in the name of the father, the son...
    bless $this, $class;

  }
  else {
    croak $CLASS. "::new() Could not locate the taxonomy data file!\n";
  }

  $this->{'isACache'} = {};

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=begin

=over 4

=item my @speciesLineage = getLineage( $obj, $species );

Return an array of species names which defines the lineage 
from $species to the ancestral "root" of the Taxonomy database.
The array is filed with $species at element 0 and "root" at 
element n.  A null list is returned if $species is not found
in the database.

=back

=cut

##-------------------------------------------------------------------------##
sub getLineage {
  my $this    = shift;
  my $species = shift;

  $species =~ s/[\<\>,'\)\(]//g;
  $species = lc( $species );
  my @lineage = ();

  if ( defined $this->{'synonyms'}->{$species} ) {
    my $id = $this->{'synonyms'}->{$species};

    push @lineage, $this->{'taxonomy'}->[ $id ]->{'name'};

    while ( $this->{'taxonomy'}->[ $id ]->{'parent'} != $id ) {
      $id = $this->{'taxonomy'}->[ $id ]->{'parent'};
      push @lineage, $this->{'taxonomy'}->[ $id ]->{'name'};
    }
  }

  return @lineage;

}

##-------------------------------------------------------------------------##

=begin

=over 4

=item my @children = getChildren( $obj, $clade );

Return an array of species names which are direct children
of the clade specified.  A null list is returned if 
clade is not found in the database.

=back

=cut

##-------------------------------------------------------------------------##
sub getChildren {
  my $this    = shift;
  my $species = shift;

  $species =~ s/[\<\>,'\)\(]//g;
  $species = lc( $species );

  my @children = ();
  if ( defined $this->{'synonyms'}->{$species} ) {
    my $id = $this->{'synonyms'}->{$species};

    # Gather your children around you
    for ( my $i = 0 ; $i <= $#{ $this->{'taxonomy'} } ; $i++ ) {
      next if ( !defined $this->{'taxonomy'}->[ $i ]->{'name'} );
      if ( $this->{'taxonomy'}->[ $i ]->{'parent'} == $id ) {
        push @children, $this->{'taxonomy'}->[ $i ]->{'name'};
      }
    }
  }

  return @children;
}

##-------------------------------------------------------------------------##

=begin

=over 4

=item my $species = isSpecies( $obj, $species );

Return species scientific name if the $species (synonym or otherwise) exists
in the database null otherwise.

NOTE: The NCBI Taxonomy database allows many special characters in 
      authoratative name of the species.  We cannot expect these
      to be used as-is in every part of the RM software.  We there
      fore sanitive the names in a predictable way:

      ">", "<", ",", ")", and "(", and "'" get stripped from the name
      Leading spaces get stripped
      Trailing spaces get stripped

=back

=cut

##-------------------------------------------------------------------------##
sub isSpecies {
  my $this    = shift;
  my $species = shift;

  $species =~ s/[\<\>,'\)\(]//g;
  $species =~ s/^\s+(\S.*)$/$1/;
  $species =~ s/(.*\S)\s+$/$1/;
  if ( exists $this->{'synonyms'}->{ lc( $species ) } ) {
    return (
             lc(
                $this->{'taxonomy'}->[ $this->{'synonyms'}->{ lc( $species ) } ]
                    ->{'name'}
             )
    );
  }

}

##-------------------------------------------------------------------------##

=begin

=over 4

=item my @species = getSimilarSoundingSpecies( $obj, $species );

Use the soundex algorithm to find similar sounding species in the
Taxonomy database.  

=back

=cut

##-------------------------------------------------------------------------##
sub getSimilarSoundingSpecies {
  my $this     = shift;
  my $species  = shift;
  my $maxCount = shift;

  my $count = 1;
  $species =~ s/[\<\>,'\)\(]//g;
  $species =~ s/^\s+(\S.*)$/$1/;
  $species =~ s/(.*\S)\s+$/$1/;
  my $speciesCode    = soundex( lc( $species ) );
  my @codes          = soundex( keys( %{ $this->{'synonyms'} } ) );
  my @similarSpecies = ();
  for ( my $i = 0 ; $i <= $#codes ; $i++ ) {

    if ( $codes[ $i ] eq $speciesCode ) {
      push @similarSpecies, ( keys( %{ $this->{'synonyms'} } ) )[ $i ];
      $count++;
      last if ( defined $maxCount && $count > $maxCount );
    }
  }
  return @similarSpecies;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=begin

=over 4 

=item my $bool = isA( $obj, $species, $group );

Returns true (1) if $species isA member of (in the clade of) $group,
false (0) if not, and (-1) if either $species or $group is not in
the Taxonomy database.

=back

=cut

##-------------------------------------------------------------------------##
sub isA {
  my $this    = shift;
  my $species = shift;
  my $group   = shift;

  $species =~ s/[\<\>,'\)\(]//g;
  $group   =~ s/[\<\>,'\)\(]//g;
  my $taxonomy = $this->{'taxonomy'};

  # Make sure group is not a synonym
  if ( exists $this->{'synonyms'}->{ lc( $group ) } ) {
    my $groupID = $this->{'synonyms'}->{ lc( $group ) };
    $group = lc( $taxonomy->[ $groupID ]->{'name'} );
  }
  else {
    return ( -1 );
  }

  unless ( exists $this->{'synonyms'}->{ lc( $species ) } ) {
    return ( -1 );
  }

  my $ID = $this->{'synonyms'}->{ lc( $species ) };

  if ( not defined $this->{'isACache'}->{ lc( $species ) } ) {

    # Must traverse the tree
    $this->{'isACache'}->{ lc( $species ) } =
        { lc( $taxonomy->[ $ID ]->{'name'} ) => 1 };
    while ( $taxonomy->[ $ID ]->{'parent'} != $ID ) {
      $ID = $taxonomy->[ $ID ]->{'parent'};
      $this->{'isACache'}->{ lc( $species ) }
          ->{ lc( $taxonomy->[ $ID ]->{'name'} ) } = 1;
    }
  }
  return 0
      if (
        not defined $this->{'isACache'}->{ lc( $species ) }->{ lc( $group ) } );

  return 1;

}

##-------------------------------------------------------------------------##

=begin

=over 4 

=item my $bool = predates( $obj, $species, $ancSpecies1, $ancSpecies2 );

Returns true (1) if $species predates the common ancester of 
ancSpecies1 and ancSpecies2 false (0) if not, and (-1) if either 
$species or $ancSpecies1/2 is not in the Taxonomy database.

=back

=cut

##-------------------------------------------------------------------------##
sub predates {
  my $this        = shift;
  my $species     = shift;
  my $ancSpecies1 = shift;
  my $ancSpecies2 = shift;

  my $taxonomy = $this->{'taxonomy'};

  $species     =~ s/[\<\>,'\)\(]//g;
  $ancSpecies1 =~ s/[\<\>,'\)\(]//g;
  $ancSpecies2 =~ s/[\<\>,'\)\(]//g;

  # Make sure ancestral species are not a synonyms
  my $ancSpecies1ID = -1;
  if ( exists $this->{'synonyms'}->{ lc( $ancSpecies1 ) } ) {
    $ancSpecies1ID = $this->{'synonyms'}->{ lc( $ancSpecies1 ) };
    $ancSpecies1   = lc( $taxonomy->[ $ancSpecies1ID ]->{'name'} );
  }
  else {
    return ( -1 );
  }
  my $ancSpecies2ID = -1;
  if ( exists $this->{'synonyms'}->{ lc( $ancSpecies2 ) } ) {
    $ancSpecies2ID = $this->{'synonyms'}->{ lc( $ancSpecies2 ) };
    $ancSpecies2   = lc( $taxonomy->[ $ancSpecies2ID ]->{'name'} );
  }
  else {
    return ( -1 );
  }

  unless ( exists $this->{'synonyms'}->{ lc( $species ) } ) {
    return ( -1 );
  }

  my $ID = $this->{'synonyms'}->{ lc( $species ) };

  # Find common ancestor to ancSpeces 1&2
  my $commonName = "";
  my $commonID   = -1;
  if ( $this->isA( $ancSpecies1, $ancSpecies2 ) ) {

    # Common ancestor is ancSpecies2
    $commonName = $ancSpecies2;
    $commonID   = $ancSpecies2ID;
  }
  elsif ( $this->isA( $ancSpecies2, $ancSpecies1 ) ) {

    # Common ancestor is ancSpecies1
    $commonName = $ancSpecies1;
    $commonID   = $ancSpecies1ID;
  }
  else {

    # Common ancestor needs to be found
    # Must traverse the tree

    $commonID   = $ancSpecies1ID;
    $commonName = $ancSpecies1;

    #print "Looking for common ancestor to $commonName\n";
    while ( $taxonomy->[ $commonID ]->{'parent'} != $commonID ) {
      $commonID   = $taxonomy->[ $commonID ]->{'parent'};
      $commonName = lc( $taxonomy->[ $commonID ]->{'name'} );

      #print " ..common to: $commonName\n";
      if ( $this->isA( $ancSpecies2, $commonName ) == 1 ) {
        last;
      }
    }

  }
  return ( $this->isA( $commonName, $species ) );
}

##-------------------------------------------------------------------------##

=begin

=over 4 

=item my $bool = getTree( $obj, $rootSpecies, $fullLineageHashRef, 
                          [ $actualSpeciesHashRef ] );

  $rootSpecies:          The species name to use as the root of the tree.
  $fullLineageHashRef:   A hash which contains all the species/groups/taxa
                         which are to be included in the output.  The keys
                         are species ID and the values are the total number
                         elements contained at this level.
  $actualSpeciesHashRef: A hash which contains only the actual species
                         elements.  This optional parameter is used to 
                         put a second number ( the actual number of elements )
                         at one level in the tree.

Returns a string containing the tree using a depth first search and
indented to designates levels in the tree.

=back

=cut

##-------------------------------------------------------------------------##
sub getTree {
  my $this                 = shift;
  my $rootSpecies          = shift;
  my $querySpeciesHashRef  = shift;
  my $actualSpeciesHashRef = shift;

  $rootSpecies =~ s/[\<\>,'\)\(]//g;

  # Make parent->children relation
  my %children = ();
  for ( my $i = 0 ; $i <= $#{ $this->{'taxonomy'} } ; $i++ ) {
    next if ( !defined $this->{'taxonomy'}->[ $i ]->{'name'} );
    next
        if (
        !defined $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $i ]->{'name'} }
        );
    push @{ $children{ $this->{'taxonomy'}->[ $i ]->{'parent'} } }, $i;
  }

  my $speciesID = -1;
  if ( exists $this->{'synonyms'}->{ lc( $rootSpecies ) } ) {
    $speciesID = $this->{'synonyms'}->{ lc( $rootSpecies ) };
  }
  else {
    return ( "Species name \"$rootSpecies\" not defined!" );
  }

  # Perform a depth first search of the tree
  my @T      = ( @{ $children{$speciesID} } );
  my $outStr = "";
  my $level  = 0;
  while ( scalar( @T ) > 0 ) {
    $speciesID = pop @T;
    while ( $speciesID eq "*" ) {
      $speciesID = pop @T;
      $level--;
    }
    last if ( scalar( @T ) < 1 );
    if (
       $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]->{'name'} } >
       0 )
    {
      $outStr .=
            " " x ( $level * 4 )
          . $this->{'taxonomy'}->[ $speciesID ]->{'name'} . " [ "
          . $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]
            ->{'name'} };
      if (
        defined
        $actualSpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]->{'name'} }
          )
      {
        $outStr .= " {"
            . $actualSpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]
              ->{'name'} }
            . "} ]\n";
      }
      else {
        $outStr .= " ]\n";
      }
    }
    if (    exists $children{$speciesID}
         && @{ $children{$speciesID} }
         && scalar( @{ $children{$speciesID} } ) > 0 )
    {
      $level++;
      push @T, "*";
      push @T, ( @{ $children{$speciesID} } );
    }
  }
  return $outStr;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

sub _buildFromStructures {
  my $this             = shift;
  my $taxonomyArrayRef = shift;
  my $synonymHashRef   = shift;
  my $dataFile         = shift;
  my $supplSynonyms    = shift;

  my $synonyms = $this->{'synonyms'} = $synonymHashRef;
  my $taxonomy = $this->{'taxonomy'} = $taxonomyArrayRef;

  # Supplement NCBI synonyms
  foreach my $syn ( keys( %{$supplSynonyms} ) ) {
    if ( exists $synonyms->{ $supplSynonyms->{$syn} } ) {
      $synonyms->{$syn} = $synonyms->{ $supplSynonyms->{$syn} };
    }
    else {
      print "Warning: Species $supplSynonyms->{$syn} does not exist. "
          . "Cannot create supplemental synonym $syn!\n";
    }
  }

  print "Maximum NCBI tax_id = " . $#{$taxonomy} . "\n";
  my $realTotal = 0;
  for ( my $i = 0 ; $i <= $#{$taxonomy} ; $i++ ) {
    $realTotal++
        if (    defined $taxonomy->[ $i ]
             && defined $taxonomy->[ $i ]->{'name'} );
  }
  print "taxonomy.dat total entries = $realTotal\n";
  nstore $this, $dataFile;

}

sub _buildFromNCBIDB {
  my $this          = shift;
  my $namesDmpFile  = shift;
  my $nodesDmpFile  = shift;
  my $dataFile      = shift;
  my $supplSynonyms = shift;
  my $RMDatabase    = shift;

  open TNAMES, "<$namesDmpFile"
      or croak $CLASS
      . "::buildFromNCBIDB: Could not open NCBI names dmp file: "
      . "$namesDmpFile\n";
  my $synonyms = $this->{'synonyms'} = {};
  my $taxonomy = $this->{'taxonomy'} = [];

  while ( <TNAMES> ) {
    s/\t//g;
    my @fields = split /\|/;
    if ( $fields[ 3 ] eq "scientific name" ) {
      my $name = lc( $fields[ 1 ] );
      $name = lc( $fields[ 2 ] ) if ( $fields[ 2 ] =~ /\S+/ );
      $name =~ s/[\<\>,'\)\(]//g;

      if ( defined( $synonyms->{$name} ) ) {
        print "Warning name $name is already defined!\n";
      }
      $synonyms->{$name} = $fields[ 0 ];
      if ( not defined $taxonomy->[ $fields[ 0 ] ] ) {
        $taxonomy->[ $fields[ 0 ] ] = {
                                        "id"   => $fields[ 0 ],
                                        "name" => $name
        };
      }
      else {
        print $CLASS
            . "::buildFromNCBIDB: Warning taxonomy record "
            . "for $name is already defined!\n";
      }
    }
    elsif (    $fields[ 3 ] eq "synonym"
            || $fields[ 3 ] eq "common name"
            || $fields[ 3 ] eq "genbank common name"
            || $fields[ 3 ] eq "genbank synonym" )
    {
      my $name = lc( $fields[ 1 ] );
      $name = lc( $fields[ 2 ] ) if ( $fields[ 2 ] =~ /\S+/ );
      if ( defined( $synonyms->{$name} ) ) {
        print "Warning synonym $name is already defined!\n";
      }
      $synonyms->{$name} = $fields[ 0 ];
    }
  }
  close TNAMES;

  # Supplement NCBI synonyms
  foreach my $syn ( keys( %{$supplSynonyms} ) ) {
    if ( exists $synonyms->{ $supplSynonyms->{$syn} } ) {
      $synonyms->{$syn} = $synonyms->{ $supplSynonyms->{$syn} };
    }
    else {
      print "Warning: Species $supplSynonyms->{$syn} does not exist. "
          . "Cannot create supplemental synonym $syn!\n";
    }
  }

  open TNODES, "<$nodesDmpFile"
      or croak $CLASS
      . "::buildFromNCBIDB: Could not open NCBI nodes dmp file: "
      . "$nodesDmpFile\n";
  while ( <TNODES> ) {
    s/\t//g;
    my @fields = split /\|/;
    if ( defined $taxonomy->[ $fields[ 0 ] ] ) {
      my $recRef = $taxonomy->[ $fields[ 0 ] ];
      $recRef->{'parent'} = $fields[ 1 ];
      $recRef->{'rank'}   = $fields[ 2 ];
    }
  }
  close TNODES;

  if ( defined $RMDatabase ) {
    my %RMSpecies = ();
    my $db        = EMBL->new( fileName => $RMDatabase );
    my $seqCount  = $db->getRecordCount();
    for ( my $i = 0 ; $i < $seqCount ; $i++ ) {
      my $record = $db->getRecord( $i );
      foreach my $name ( $record->getRMSpeciesArray() ) {
        $name =~ s/_/ /g;
        $name = lc( $name );
        next if ( $name =~ /root/ );
        $RMSpecies{$name} = 1;
      }
    }
    undef $db;

    # Create a children list
    my %children = ();
    print "  - Reversing the lookup tree\n";
    foreach my $tax ( @{$taxonomy} ) {
      my $id  = $tax->{'id'};
      my $pid = $tax->{'parent'};
      if ( defined $pid ) {
        $children{$pid}->{"$id"} = 1;
      }
    }

    # Add all children
    my @childlist           = ();
    my $unresolvedRMSpecies = 0;
    print "  - Adding children of RM's species\n";
    foreach my $name ( keys( %RMSpecies ) ) {
      print "    -- adding $name\n";
      if ( defined $synonyms->{$name} ) {
        my $id = $synonyms->{$name};
        push @childlist, $id;
      }
      else {
        print $CLASS
            . "::_buildFromNCBIDB: Warning $name in the "
            . "RepeatMasker database is not known to the NCBI "
            . "Taxonomy database!\n";
        $unresolvedRMSpecies++;
      }
    }
    croak $CLASS
        . "::_buildFromNCBIDB: There were unresolved "
        . "RepeatMasker species cannot continue building "
        . "taxonomy database until these are fixed.\n"
        if ( $unresolvedRMSpecies );

    # Adding all children of the seeeded list
    my %goodTaxIDs = ();

    while ( @childlist ) {
      my $id = pop @childlist;
      next if ( defined $goodTaxIDs{$id} );
      $goodTaxIDs{$id} = 1;
      push @childlist, keys( %{ $children{$id} } );
    }
    undef %children;
    undef @childlist;

    print "  - Adding parents of RM's species\n";
    foreach my $species ( keys( %RMSpecies ) ) {
      my $id = $synonyms->{$species};
      while ( $taxonomy->[ $id ]->{'parent'} != $id ) {
        $id = $taxonomy->[ $id ]->{'parent'};
        $goodTaxIDs{$id} = 1;
      }
    }
    undef %RMSpecies;

    foreach my $name ( keys( %{$synonyms} ) ) {
      if ( not defined $goodTaxIDs{ $synonyms->{$name} } ) {
        $taxonomy->[ $synonyms->{$name} ] = undef;
        delete $synonyms->{$name};
      }
    }
    undef %goodTaxIDs;

  }

  #   id, parent_id, name
  #   id, ( synonyms )
  #   @taxomony->[id]->{'id' => #,
  #                     'parent' => #,
  #                     'name' => "",
  #                     'rank' => };
  #   %synonyms->{name} = id#

  print "Total entries = " . $#{$taxonomy} . "\n";
  my $realTotal = 0;
  for ( my $i = 0 ; $i <= $#{$taxonomy} ; $i++ ) {
    $realTotal++
        if (    defined $taxonomy->[ $i ]
             && defined $taxonomy->[ $i ]->{'name'} );
  }
  print "Actual entries = $realTotal\n";
  print "Dumper: " . Dumper( $taxonomy ) . "\n";
  nstore $this, $dataFile;

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

