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

Taxonomy

=head1 SYNOPSIS

use Taxonomy;

Usage: 

my $taxDB = Taxonomy->new( famdbfile=>"Dfam.h5" );

if ( $taxDB->isA( "Mouse", "Mammalia" ) ) {
  print "A Mouse is a Mammal!\n";
}else {
  print "Error...a Mouse should be Mammal!\n";
}

=head1 DESCRIPTION

A general object for querying the NCBI Taxonomy database.
This is now a wrapper around functionality provided by the
FamDB format.

The NCBI Taxonomy database is filtered for use solely with Dfam or the
RepeatMasker library.  Names in the taxonomy database are removed
if they are not related ( child or ancestor ) to a species contained 
in the RepeatMasker database.

=head1 INSTANCE METHODS

=cut 

package Taxonomy;
use strict;
use FindBin;
use Data::Dumper;
use FastaDB;
use EMBL;
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
my $CLASS   = "Taxonomy";

my $FAMDB = "$FindBin::Bin/famdb.py";

##-------------------------------------------------------------------------##

=begin

=over 4

=item my $instance = Taxonomy->new( famdbfile=>"filename" );

Construct a new Taxonomy object.  Use the FamDB file
specified for queries.

=back

=cut

my %supplementalSynonyms = (
                             # additional common names / plural vs singular,
                             # not in the NCBI database
                             'carnivore'   => "carnivora",
                             'artiodactyl' => "artiodactyla",
                             'puffer'      => "takifugu",
                             'fruitfly'    => "drosophila_flies_genus",
                             'cionaint'    => "ciona intestinalis",
                             'cionasav'    => "ciona savignyi",
                             'diatom'      => "bacillariophyta",
                             'rat'         => "rattus",
                             'mammal'      => "mammalia",

                             # disambiguate some species to specific taxa
                             'anopheles'   => "anopheles_genus",
                             'drosophila'  => "drosophila_flies_genus",
                             'worm'        => "caenorhabditis",
                             'mustard'     => "arabidopsis",
                             'corn'        => "zea",
                             'elegans'     => "caenorhabditis elegans",
                             'maize'       => "zea",
                             "mus"         => "mus musculus",
                             "mouse"       => "mus musculus",
);

##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  my $this = {};

  if ( defined $nameValuePairs{'famdbfile'}
          && -s $nameValuePairs{'famdbfile'} )
  {

    # store the database filename to use later
    $this = {
      famdbfile => $nameValuePairs{'famdbfile'},
      isACache => {},
    };

    # Bless this hash in the name of the father, the son...
    bless $this, $class;

  }
  else {
    croak $CLASS. "::new() needs a path for a famdb file!\n";
  }

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

  $species =~ s/^\s+(\S.*)$/$1/;
  $species =~ s/(.*\S)\s+$/$1/;
  $species = lc($species);
  $species = $supplementalSynonyms{$species} if exists $supplementalSynonyms{$species} ;

  my @lineage = ();

  my $result = $this->_invokeFamDB([ "lineage", "--ancestors", "--format=semicolon", $species ]);

  if ( $result =~ /(\d+):\s*(.*)\s*\[(\d+)\]/ ) {
    my ( $taxId, $lineage, $count ) = ( $1, $2, $3 );
    $lineage =~ s/^\s*|\s*$//g;
    @lineage = split ';', $lineage;
  }

  return @lineage;
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

  my $result = $this->_invokeFamDB([ "lineage", "--format=semicolon", $species ]);

  $species =~ s/^\s+(\S.*)$/$1/;
  $species =~ s/(.*\S)\s+$/$1/;
  $species = lc($species);
  $species = $supplementalSynonyms{$species} if exists $supplementalSynonyms{$species} ;

  if ( $result =~ /(\d+):\s*(.*)\s*\[(\d+)\]/ ) {
    my $lineage = $2;
    $lineage =~ s/^\s*|\s*$//g;
    my @lineage = split ';', $lineage;
    my $species = $lineage[$#lineage];

    $species = _sanitizeName( $species );

    return lc($species);
  }
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

  if ( ! exists $this->{'isACache'}->{ $species } ) {
    $this->{'isACache'}->{ $species } = {};
  }

  if ( ! exists $this->{'isACache'}->{ $species }->{ $group } ) {
    my @speciesLineage = $this->getLineage( $species );
    my @groupLineage = $this->getLineage( $group );
    my $groupName = $groupLineage[$#groupLineage];

    my $result;

    if ( $#speciesLineage && $#groupLineage ) {
      foreach my $ancestor ( @speciesLineage ) {
        if ( lc($ancestor) eq lc($groupName) ) {
          $result = 1;
          last;
        }
      }
      $result = 0 if not defined $result;
    } else {
      $result = -1;
    }

    $this->{'isACache'}->{ $species }->{ $group } = $result;
  }

  return $this->{'isACache'}->{ $species }->{ $group };
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

  my @lineage1 = $this->getLineage( $ancSpecies1 );
  my @lineage2 = $this->getLineage( $ancSpecies2 );

  if ($#lineage1 < 1 || $#lineage2 < 1) {
    return -1;
  }

  my $commonAnc = undef;

  my $i = 0;
  while ( $i < $#lineage1 && $i < $#lineage2 ) {
    last if !( $lineage1[$i] eq $lineage2[$i] );
    $commonAnc = $lineage1[$i];
  }

  return $this->isA( $commonAnc, $species );
}
##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

sub _sanitizeName {
    my $name = shift;
    $name =~ s/[\<\>,'\)\(]//g;
    return $name;
}

# Invokes famdb, returning undef if there is any problem
# with the species name (i.e. missing, ambiguous).
#
# $this->_invokeFamDB ( [ "lineage", "--ancestors", $taxId ] );
# $this->_invokeFamDB ( [ "names", $speciesName ] );
sub _invokeFamDB {
  my $this = shift;
  my $args = shift;

  my $dbfile = $this->{famdbfile};

  my $args_quoted = "";
  for my $arg (@{$args}) {
    my $argq = $arg;
    $argq =~ s/'/'"'"'/g;
    $args_quoted .= " '$argq'";
  }

  my $result = `$FAMDB --file $dbfile $args_quoted 2>&1`;

  if (    $result =~ /^\s*no results/i
       || $result =~ /^\s*no species/i
       || $result =~ /^\s*ambiguous/i ) {
    return undef;
  } else {
    return $result;
  }
}

## DEPRECATED / DELETED methods.
## These used to be implemented on top of the hash-based structures,
## and would need to be reimplemented on top of / in famdb but
## are not used by any other code.

###-------------------------------------------------------------------------##
#
#=begin
#
#=over 4
#
#=item my @children = getChildren( $obj, $clade );
#
#Return an array of species names which are direct children
#of the clade specified.  A null list is returned if 
#clade is not found in the database.
#
#=back
#
#=cut
#
###-------------------------------------------------------------------------##
#sub getChildren {
#  my $this    = shift;
#  my $species = shift;
#
#  $species =~ s/[\<\>,'\)\(]//g;
#  $species = lc( $species );
#
#  my @children = ();
#  if ( defined $this->{'synonyms'}->{$species} ) {
#    my $id = $this->{'synonyms'}->{$species};
#
#    # Gather your children around you
#    for ( my $i = 0 ; $i <= $#{ $this->{'taxonomy'} } ; $i++ ) {
#      next if ( !defined $this->{'taxonomy'}->[ $i ]->{'name'} );
#      if ( $this->{'taxonomy'}->[ $i ]->{'parent'} == $id ) {
#        push @children, $this->{'taxonomy'}->[ $i ]->{'name'};
#      }
#    }
#  }
#
#  return @children;
#}
#
#
###-------------------------------------------------------------------------##
#
#=begin
#
#=over 4
#
#=item my @species = getSimilarSoundingSpecies( $obj, $species );
#
#Use the soundex algorithm to find similar sounding species in the
#Taxonomy database.  
#
#=back
#
#=cut
#
###-------------------------------------------------------------------------##
#sub getSimilarSoundingSpecies {
#  my $this     = shift;
#  my $species  = shift;
#  my $maxCount = shift;
#
#  my $count = 1;
#  $species = _sanitizeName ( $species );
#  $species =~ s/^\s+(\S.*)$/$1/;
#  $species =~ s/(.*\S)\s+$/$1/;
#  my $speciesCode    = soundex( lc( $species ) );
#  my @codes          = soundex( keys( %{ $this->{'synonyms'} } ) );
#  my @similarSpecies = ();
#  for ( my $i = 0 ; $i <= $#codes ; $i++ ) {
#
#    if ( $codes[ $i ] eq $speciesCode ) {
#      push @similarSpecies, ( keys( %{ $this->{'synonyms'} } ) )[ $i ];
#      $count++;
#      last if ( defined $maxCount && $count > $maxCount );
#    }
#  }
#  return @similarSpecies;
#}
#
###-------------------------------------------------------------------------##
#
#=begin
#
#=over 4 
#
#=item my $bool = getTree( $obj, $rootSpecies, $fullLineageHashRef, 
#                          [ $actualSpeciesHashRef ] );
#
#  $rootSpecies:          The species name to use as the root of the tree.
#  $fullLineageHashRef:   A hash which contains all the species/groups/taxa
#                         which are to be included in the output.  The keys
#                         are species ID and the values are the total number
#                         elements contained at this level.
#  $actualSpeciesHashRef: A hash which contains only the actual species
#                         elements.  This optional parameter is used to 
#                         put a second number ( the actual number of elements )
#                         at one level in the tree.
#
#Returns a string containing the tree using a depth first search and
#indented to designates levels in the tree.
#
#=back
#
#=cut
#
###-------------------------------------------------------------------------##
#sub getTree {
#  my $this                 = shift;
#  my $rootSpecies          = shift;
#  my $querySpeciesHashRef  = shift;
#  my $actualSpeciesHashRef = shift;
#
#  $rootSpecies =~ s/[\<\>,'\)\(]//g;
#
#  # Make parent->children relation
#  my %children = ();
#  for ( my $i = 0 ; $i <= $#{ $this->{'taxonomy'} } ; $i++ ) {
#    next if ( !defined $this->{'taxonomy'}->[ $i ]->{'name'} );
#    next
#        if (
#        !defined $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $i ]->{'name'} }
#        );
#    push @{ $children{ $this->{'taxonomy'}->[ $i ]->{'parent'} } }, $i;
#  }
#
#  my $speciesID = -1;
#  if ( exists $this->{'synonyms'}->{ lc( $rootSpecies ) } ) {
#    $speciesID = $this->{'synonyms'}->{ lc( $rootSpecies ) };
#  }
#  else {
#    return ( "Species name \"$rootSpecies\" not defined!" );
#  }
#
#  # Perform a depth first search of the tree
#  my @T      = ( @{ $children{$speciesID} } );
#  my $outStr = "";
#  my $level  = 0;
#  while ( scalar( @T ) > 0 ) {
#    $speciesID = pop @T;
#    while ( $speciesID eq "*" ) {
#      $speciesID = pop @T;
#      $level--;
#    }
#    last if ( scalar( @T ) < 1 );
#    if (
#       $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]->{'name'} } >
#       0 )
#    {
#      $outStr .=
#            " " x ( $level * 4 )
#          . $this->{'taxonomy'}->[ $speciesID ]->{'name'} . " [ "
#          . $querySpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]
#            ->{'name'} };
#      if (
#        defined
#        $actualSpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]->{'name'} }
#          )
#      {
#        $outStr .= " {"
#            . $actualSpeciesHashRef->{ $this->{'taxonomy'}->[ $speciesID ]
#              ->{'name'} }
#            . "} ]\n";
#      }
#      else {
#        $outStr .= " ]\n";
#      }
#    }
#    if (    exists $children{$speciesID}
#         && @{ $children{$speciesID} }
#         && scalar( @{ $children{$speciesID} } ) > 0 )
#    {
#      $level++;
#      push @T, "*";
#      push @T, ( @{ $children{$speciesID} } );
#    }
#  }
#  return $outStr;
#}
#
#sub _buildFromStructures {
#  my $this             = shift;
#  my $taxonomyArrayRef = shift;
#  my $synonymHashRef   = shift;
#  my $dataFile         = shift;
#  my $supplSynonyms    = shift;
#
#  my $synonyms = $this->{'synonyms'} = $synonymHashRef;
#  my $taxonomy = $this->{'taxonomy'} = $taxonomyArrayRef;
#
#  # Supplement NCBI synonyms
#  foreach my $syn ( keys( %{$supplSynonyms} ) ) {
#    if ( exists $synonyms->{ $supplSynonyms->{$syn} } ) {
#      $synonyms->{$syn} = $synonyms->{ $supplSynonyms->{$syn} };
#    }
#    else {
#      print "Warning: Species $supplSynonyms->{$syn} does not exist. "
#          . "Cannot create supplemental synonym $syn!\n";
#    }
#  }
#
#  print "Maximum NCBI tax_id = " . $#{$taxonomy} . "\n";
#  my $realTotal = 0;
#  for ( my $i = 0 ; $i <= $#{$taxonomy} ; $i++ ) {
#    $realTotal++
#        if (    defined $taxonomy->[ $i ]
#             && defined $taxonomy->[ $i ]->{'name'} );
#  }
#  print "taxonomy.dat total entries = $realTotal\n";
#  nstore $this, $dataFile;
#
#}
#
#sub _buildFromNCBIDB {
#  my $this          = shift;
#  my $namesDmpFile  = shift;
#  my $nodesDmpFile  = shift;
#  my $dataFile      = shift;
#  my $supplSynonyms = shift;
#  my $RMDatabase    = shift;
#
#  open TNAMES, "<$namesDmpFile"
#      or croak $CLASS
#      . "::buildFromNCBIDB: Could not open NCBI names dmp file: "
#      . "$namesDmpFile\n";
#  my $synonyms = $this->{'synonyms'} = {};
#  my $taxonomy = $this->{'taxonomy'} = [];
#
#  while ( <TNAMES> ) {
#    s/\t//g;
#    my @fields = split /\|/;
#    if ( $fields[ 3 ] eq "scientific name" ) {
#      my $name = lc( $fields[ 1 ] );
#      $name = lc( $fields[ 2 ] ) if ( $fields[ 2 ] =~ /\S+/ );
#      $name =~ s/[\<\>,'\)\(]//g;
#
#      if ( defined( $synonyms->{$name} ) ) {
#        print "Warning name $name is already defined!\n";
#      }
#      $synonyms->{$name} = $fields[ 0 ];
#      if ( not defined $taxonomy->[ $fields[ 0 ] ] ) {
#        $taxonomy->[ $fields[ 0 ] ] = {
#                                        "id"   => $fields[ 0 ],
#                                        "name" => $name
#        };
#      }
#      else {
#        print $CLASS
#            . "::buildFromNCBIDB: Warning taxonomy record "
#            . "for $name is already defined!\n";
#      }
#    }
#    elsif (    $fields[ 3 ] eq "synonym"
#            || $fields[ 3 ] eq "common name"
#            || $fields[ 3 ] eq "genbank common name"
#            || $fields[ 3 ] eq "genbank synonym" )
#    {
#      my $name = lc( $fields[ 1 ] );
#      $name = lc( $fields[ 2 ] ) if ( $fields[ 2 ] =~ /\S+/ );
#      if ( defined( $synonyms->{$name} ) ) {
#        print "Warning synonym $name is already defined!\n";
#      }
#      $synonyms->{$name} = $fields[ 0 ];
#    }
#  }
#  close TNAMES;
#
#  # Supplement NCBI synonyms
#  foreach my $syn ( keys( %{$supplSynonyms} ) ) {
#    if ( exists $synonyms->{ $supplSynonyms->{$syn} } ) {
#      $synonyms->{$syn} = $synonyms->{ $supplSynonyms->{$syn} };
#    }
#    else {
#      print "Warning: Species $supplSynonyms->{$syn} does not exist. "
#          . "Cannot create supplemental synonym $syn!\n";
#    }
#  }
#
#  open TNODES, "<$nodesDmpFile"
#      or croak $CLASS
#      . "::buildFromNCBIDB: Could not open NCBI nodes dmp file: "
#      . "$nodesDmpFile\n";
#  while ( <TNODES> ) {
#    s/\t//g;
#    my @fields = split /\|/;
#    if ( defined $taxonomy->[ $fields[ 0 ] ] ) {
#      my $recRef = $taxonomy->[ $fields[ 0 ] ];
#      $recRef->{'parent'} = $fields[ 1 ];
#      $recRef->{'rank'}   = $fields[ 2 ];
#    }
#  }
#  close TNODES;
#
#  if ( defined $RMDatabase ) {
#    my %RMSpecies = ();
#    my $db        = EMBL->new( fileName => $RMDatabase );
#    my $seqCount  = $db->getRecordCount();
#    for ( my $i = 0 ; $i < $seqCount ; $i++ ) {
#      my $record = $db->getRecord( $i );
#      foreach my $name ( $record->getRMSpeciesArray() ) {
#        $name =~ s/_/ /g;
#        $name = lc( $name );
#        next if ( $name =~ /root/ );
#        $RMSpecies{$name} = 1;
#      }
#    }
#    undef $db;
#
#    # Create a children list
#    my %children = ();
#    print "  - Reversing the lookup tree\n";
#    foreach my $tax ( @{$taxonomy} ) {
#      my $id  = $tax->{'id'};
#      my $pid = $tax->{'parent'};
#      if ( defined $pid ) {
#        $children{$pid}->{"$id"} = 1;
#      }
#    }
#
#    # Add all children
#    my @childlist           = ();
#    my $unresolvedRMSpecies = 0;
#    print "  - Adding children of RM's species\n";
#    foreach my $name ( keys( %RMSpecies ) ) {
#      print "    -- adding $name\n";
#      if ( defined $synonyms->{$name} ) {
#        my $id = $synonyms->{$name};
#        push @childlist, $id;
#      }
#      else {
#        print $CLASS
#            . "::_buildFromNCBIDB: Warning $name in the "
#            . "RepeatMasker database is not known to the NCBI "
#            . "Taxonomy database!\n";
#        $unresolvedRMSpecies++;
#      }
#    }
#    croak $CLASS
#        . "::_buildFromNCBIDB: There were unresolved "
#        . "RepeatMasker species cannot continue building "
#        . "taxonomy database until these are fixed.\n"
#        if ( $unresolvedRMSpecies );
#
#    # Adding all children of the seeeded list
#    my %goodTaxIDs = ();
#
#    while ( @childlist ) {
#      my $id = pop @childlist;
#      next if ( defined $goodTaxIDs{$id} );
#      $goodTaxIDs{$id} = 1;
#      push @childlist, keys( %{ $children{$id} } );
#    }
#    undef %children;
#    undef @childlist;
#
#    print "  - Adding parents of RM's species\n";
#    foreach my $species ( keys( %RMSpecies ) ) {
#      my $id = $synonyms->{$species};
#      while ( $taxonomy->[ $id ]->{'parent'} != $id ) {
#        $id = $taxonomy->[ $id ]->{'parent'};
#        $goodTaxIDs{$id} = 1;
#      }
#    }
#    undef %RMSpecies;
#
#    foreach my $name ( keys( %{$synonyms} ) ) {
#      if ( not defined $goodTaxIDs{ $synonyms->{$name} } ) {
#        $taxonomy->[ $synonyms->{$name} ] = undef;
#        delete $synonyms->{$name};
#      }
#    }
#    undef %goodTaxIDs;
#
#  }
#
#  #   id, parent_id, name
#  #   id, ( synonyms )
#  #   @taxomony->[id]->{'id' => #,
#  #                     'parent' => #,
#  #                     'name' => "",
#  #                     'rank' => };
#  #   %synonyms->{name} = id#
#
#  print "Total entries = " . $#{$taxonomy} . "\n";
#  my $realTotal = 0;
#  for ( my $i = 0 ; $i <= $#{$taxonomy} ; $i++ ) {
#    $realTotal++
#        if (    defined $taxonomy->[ $i ]
#             && defined $taxonomy->[ $i ]->{'name'} );
#  }
#  print "Actual entries = $realTotal\n";
#  print "Dumper: " . Dumper( $taxonomy ) . "\n";
#  nstore $this, $dataFile;
#
#}
#
###-------------------------------------------------------------------------##
### Serialization & Debug Routines
###-------------------------------------------------------------------------##
#
###-------------------------------------------------------------------------##
### Use: my $string = toString([$this]);
###
###      $this         : Normally passed implicitly
###
###  Returns
###
###      Uses the Data::Dumper to create a printable reprentation
###      of a data structure.  In this case the object data itself.
###
###-------------------------------------------------------------------------##
#sub toString {
#  my $this = shift;
#  my $data_dumper = new Data::Dumper( [ $this ] );
#  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
#  return $data_dumper->Dump();
#}
#
###-------------------------------------------------------------------------##
### Use: my serializeOUT( $filename );
###
###	  $filename	: A filename to be created
###
###  Returns
###
###	Uses the Data::Dumper module to save out the data
###	structure as a text file.  This text file can be
###	read back into an object of this type.
###
###-------------------------------------------------------------------------##
#sub serializeOUT {
#  my $this     = shift;
#  my $fileName = shift;
#
#  my $data_dumper = new Data::Dumper( [ $this ] );
#  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
#  open OUT, ">$fileName";
#  print OUT $data_dumper->Dump();
#  close OUT;
#}
#
###-------------------------------------------------------------------------##
### Use: my serializeIN( $filename );
###
###	$filename	: A filename containing a serialized object
###
###  Returns
###
###	Uses the Data::Dumper module to read in data
###	from a serialized PERL object or data structure.
###
###-------------------------------------------------------------------------##
#sub serializeIN {
#  my $this         = shift;
#  my $fileName     = shift;
#  my $fileContents = "";
#  my $oldSep       = $/;
#  undef $/;
#  my $in;
#  open $in, "$fileName";
#  $fileContents = <$in>;
#  $/            = $oldSep;
#  close $in;
#  return eval( $fileContents );
#}

1;

