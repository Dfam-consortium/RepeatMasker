#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) queryRepeatDatabase.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to assist in querying the monolithic
##      RM repeat sequence database.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2005 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
#
# ChangeLog
#
#     $Log: queryRepeatDatabase.pl,v $
#     Revision 1.100  2017/02/01 21:01:57  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

queryRepeatDatabase.pl - Query the RepeatMasker repeat database.

=head1 SYNOPSIS

  queryRepeatDatabase.pl [-version] [-species <species> |
                                     -stage <stage num> |
                                     -class <class> |
                                     -id <id>]
                                    [-stat]
                                    [-tree]
                                    [-clade]

=head1 DESCRIPTION

  A utility script to query the RepeatMasker repeat database.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=over 4

=item -species "species name"

The full name ( case insensitive ) of the species you would like
to search for in the database.  This will return all the repeats
which would be used in a RepeatMasker search against this species.
This includes repeats contained in the clade given by "species name"
and ancestral repeats of "species name".  Lastly ubiquitous 
sequences such as RNAs and simple repeats are also included.

=item -clade 

This will modify the default behaviour of the species option and
return only the repeats which are specific to your species or any
of it descendents.  This is useful for identifying how rich
the database of repeats is for a given species/clade.  

=item -stage <stage num>

The number of the RepeatMasker stage for which you would like
repeats.  In the past these stages were individual libraries
with the following general names:

  Stage          Library
  -----          -------
   0             species.lib
  10             is.lib
  15             rodspec.lib
  20             humspec.lib
  25             simple.lib
  30             at.lib
  35             sinecutlib
  40             shortcutlib
  45             cutlib
  50             shortlib
  55             longlib
  60             mirs.lib
  65             mir.lib
  70             retrovirus.lib
  75             l1.lib

=item -class <class>

Retrieve all elements of a particular class.  For example:

  DNA
  SINE
  LINE
  LTR
  Other
  RC
  Satellite
  tRNA
  Simple_repeat
  Unknown
  snRNA

=item -id <id>

Retrieve only a single id from the database.

=item -stat

Returns statistics on the sequences

=item -tree

Prints the taxonomy tree for all species in the database.

=back

=head1 SEE ALSO

ReapeatMasker

=head1 COPYRIGHT

Copyright 2005-2011 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Getopt::Long;
use Data::Dumper;
use FastaDB;
use RmEMBL;
use DFAM;
use Taxonomy;
use File::Basename;
use LibraryUtils;

#
# Version
#
my $Version = 0.1;
my $DEBUG   = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',     # print out the version and exit
                    '-dfam',
                    '-embl',
                    '-fasta',
                    '-species=s',
                    '-clade',
                    '-stage=i',
                    '-class=s',
                    '-quiet',
                    '-id=s',
                    '-tree',
                    '-stat'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit( 1 );
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

my $fileFormat = "Unknown";
my $RMLib;
my $libVersion;
if ( exists $options{'dfam'} ) {

  # Dfam file format
  $RMLib      = "$FindBin::Bin/../Libraries/Dfam.hmm";
  $fileFormat = "dfam";
  my @values =
      LibraryUtils::validateLibraries( "$FindBin::Bin/../Libraries", "HMM" );
  $libVersion = $values[ 2 ];
}
elsif ( exists $options{'fasta'} ) {

  # Old legacy library format
  $RMLib      = "$FindBin::Bin/../Libraries/RepeatMasker.lib";
  $fileFormat = "fasta";
}
else {

  # RepeatMasker Combined Library format
  $RMLib      = "$FindBin::Bin/../Libraries/RepeatMaskerLib.embl";
  $fileFormat = "embl";
  my @values = LibraryUtils::validateLibraries( "$FindBin::Bin/../Libraries",
                                                "CONSENSUS" );
  $libVersion = $values[ 2 ];
}

if ( !-f $RMLib ) {
  print "Error!  Could not find the RepeatMasker library: $RMLib!\n";
  exit;
}

unless ( $options{'quiet'} ) {
  print STDERR "queryRepeatDatabase\n";
  print STDERR "===================\n";
  print STDERR "RepeatMasker Database: " . basename( $RMLib ) . "\n";
  print STDERR "Version: $libVersion\n";
}

# Open up the taxonomy database
my $taxDB = "$FindBin::Bin/../Libraries/taxonomy.dat";
my $tax   = Taxonomy->new( taxonomyDataFile => $taxDB );

if ( $options{'tree'} ) {
  my $db = RmEMBL->new( fileName => $RMLib );
  &_displayLibraryTaxonomy( library_ref => $db, taxonomy_ref => $tax );
  exit;
}

my $species         = $options{'species'};
my $realSpeciesName = "";
if ( $species ne "" ) {
  if ( ( $realSpeciesName = $tax->isSpecies( $species ) ) eq "" ) {
    print STDERR "Species " . $species
        . " is not in the database. "
        . "Here is a list of possible similar\nsounding substitutes:\n";
    foreach my $species ( $tax->getSimilarSoundingSpecies( $species ) ) {
      print STDERR "  $species\n";
    }
    exit;
  }
}
unless ( $options{'quiet'} ) {
  if ( $options{'clade'} ) {
    print STDERR "Clade: $species ( $realSpeciesName )\n";
  }
  elsif ( $options{'species'} ) {
    print STDERR "Species: $species ( $realSpeciesName )\n";
  }
}

my $stagePattern = "";
$stagePattern = $options{'stage'} if ( defined $options{'stage'} );
my $idPattern = "";
$idPattern = $options{'id'} if ( defined $options{'id'} );
my $classPattern = "";
$classPattern = $options{'class'} if ( defined $options{'class'} );

my $db       = undef;
my $seqCount = 0;
my @ids      = ();
if ( $fileFormat eq "fasta" ) {
  $db       = FastaDB->new( fileName => $RMLib, openMode => SeqDBI::ReadOnly );
  $seqCount = scalar( $db->getIDs() ) + 1;
  @ids      = $db->getIDs();
}
elsif ( $fileFormat eq "embl" ) {
  $db = RmEMBL->new( fileName => $RMLib );
  $seqCount = $db->getRecordCount();
}
elsif ( $fileFormat eq "dfam" ) {
  $db = DFAM->new( fileName => $RMLib );
  $seqCount = $db->getRecordCount();
}

my $totalLength       = 0;
my $totalNumber       = 0;
my $cladeCnt          = 0;
my $equalCnt          = 0;
my $cladeLen          = 0;
my $ancestCnt         = 0;
my $ancestLen         = 0;
my $equalLen          = 0;
my %uniqSpeciesLabels = ();
for ( my $i = 0 ; $i < $seqCount ; $i++ ) {

  # Look for species match
  my $isDescendant = 0;
  my $isAncestor   = 0;
  my $isEqual      = 0;

  if ( $realSpeciesName ne "" ) {
    my $match = 0;
    if ( $db->isa( "FastaDB" ) ) {
      my $descLine = $db->getDescription( $ids[ $i ] );
      while ( $descLine =~ /@([\w\.]+)/ig ) {
        my $name = $1;
        $name =~ s/_/ /g;

       #        $isDescendant = $tax->isA( $name,            $realSpeciesName );
       #        $isAncestor   = $tax->isA( $realSpeciesName, $name );
       #        if ( $options{'clade'} ) {
       #          $cladeCnt++ if ( $isDescendant == 1 );
       #          $match = 1 if ( $isDescendant == 1 );
       #          last;
       #        }
       #        else {
       #          $cladeCnt++  if ( $isDescendant == 1  && !$isAncestor == 1 );
       #          $cladeCnt++  if ( $isDescendant == 1  && $isAncestor == 1 );
       #          $ancestCnt++ if ( !$isDescendant == 1 && $isAncestor == 1 );
       #          $match = 1 if ( $isDescendant == 1 || $isAncestor == 1 );
       #          last;
       #        }
        if ( $tax->isSpecies( $name ) eq $realSpeciesName ) {
          $isEqual = 1;
        }
        elsif ( $tax->isA( $name, $realSpeciesName ) == 1 ) {
          $isDescendant = 1;
        }
        elsif ( $tax->isA( $realSpeciesName, $name ) == 1 ) {
          $isAncestor = 1;
        }

        if ( $isEqual ) {
          $equalCnt++;
          $match = 1;
          last;
        }
        elsif ( $isDescendant ) {
          $cladeCnt++;
          $match = 1;
          last;
        }
        elsif ( ( !$options{'clade'} )
                && $isAncestor )
        {
          $ancestCnt++;
          $match = 1;
          last;
        }

      }
    }
    else {
      my $record = $db->getRecord( $i );
      foreach my $name ( $record->getRMSpeciesArray() ) {
        $name =~ s/_/ /g;

        if ( $tax->isSpecies( $name ) eq $realSpeciesName ) {
          $isEqual = 1;
        }
        elsif ( $tax->isA( $name, $realSpeciesName ) == 1 ) {
          $isDescendant = 1;
        }
        elsif ( $tax->isA( $realSpeciesName, $name ) == 1 ) {
          $isAncestor = 1;
        }

        if ( $isEqual ) {
          $equalCnt++;
          $match = 1;
          last;
        }
        elsif ( $isDescendant ) {
          $cladeCnt++;
          $match = 1;
          last;
        }
        elsif ( ( !$options{'clade'} )
                && $isAncestor )
        {
          $ancestCnt++;
          $match = 1;
          last;
        }
      }
    }
    next if ( $match == 0 );
  }    # End - Look for species match

  # Look for stages
  if ( $stagePattern ne "" ) {
    my $match = 0;
    if ( $db->isa( "FastaDB" ) ) {
      my $descLine = $db->getDescription( $ids[ $i ] );
      if ( ( $stagePattern == 0 && $descLine =~ /\[S:\]/ )
           || $descLine =~ /\[S:(\d+,)*$stagePattern(,\d+)*\]/ )
      {
        $match = 1;
      }
    }
    else {
      my $record = $db->getRecord( $i );
      my @stages = $record->getRMSearchStagesArray();
      foreach my $stage ( @stages ) {
        if ( $stage eq $stagePattern ) {
          $match = 1;
          last;
        }
      }
    }
    next if ( $match == 0 );
  }    # End - Look for stage match

  # Look for id match
  if ( $idPattern ne "" ) {
    my $match = 0;
    if ( $db->isa( "FastaDB" ) ) {
      if ( $ids[ $i ] =~ /$idPattern/i ) {
        $match = 1;
      }
    }
    else {
      my $record = $db->getRecord( $i );
      if ( $record->getId() =~ /$idPattern/i ) {
        $match = 1;
      }
    }
    next if ( $match == 0 );
  }    # End - Look for id match

  # Look for class match
  if ( $classPattern ne "" ) {
    my $match = 0;
    if ( $db->isa( "FastaDB" ) ) {
      if ( $ids[ $i ] =~ /\#$classPattern/i ) {
        $match = 1;
      }
    }
    else {
      my $record = $db->getRecord( $i );
      if ( $record->getRMType() =~ /$classPattern/i ) {
        $match = 1;
      }
    }
    next if ( $match == 0 );
  }    # End - Look for class match

  if ( $db->isa( "FastaDB" ) ) {
    my $seq = $db->getSequence( $ids[ $i ] );
    if ( defined $options{'stat'} ) {
      my $seqLen = length( $seq );
      print ">" . $ids[ $i ] . "  Length = $seqLen bp\n";
      $totalLength += $seqLen;
      if ( $isAncestor ) {
        $ancestLen += $seqLen;
      }
      elsif ( $isDescendant ) {
        $cladeLen += $seqLen;
      }
      elsif ( $isEqual ) {
        $equalLen += $seqLen;
      }
      ++$totalNumber;
      $uniqSpeciesLabels{$species}++;
    }
    else {
      my $descLine = $db->getDescription( $ids[ $i ] );
      print ">" . $ids[ $i ];
      if ( $descLine ne "" ) {
        print " " . $descLine;
      }
      print "\n";
      $seq =~ s/(\S{50})/$1\n/g;
      $seq .= "\n"
          unless ( $seq =~ /.*\n+$/s );
      print $seq;
    }
  }
  elsif ( $fileFormat eq "embl" ) {
    my $record = $db->getRecord( $i );
    my $seq    = $record->getSequence();
    if ( defined $options{'stat'} ) {
      my $seqLen = length( $seq );
      print ">"
          . $record->getId() . "#"
          . $record->getRMType() . "/"
          . $record->getRMSubType()
          . "  Length = "
          . "$seqLen bp";

      foreach my $name ( $record->getRMSpeciesArray() ) {
        print " \@$name";
        $uniqSpeciesLabels{$name}++;
      }
      print "\n";

      $totalLength += $seqLen;

      if ( $isAncestor ) {
        $ancestLen += $seqLen;
      }
      elsif ( $isDescendant ) {
        $cladeLen += $seqLen;
      }
      elsif ( $isEqual ) {
        $equalLen += $seqLen;
      }
      ++$totalNumber;
    }
    else {
      print ">" . $record->getId() . "#" . $record->getRMType();
      if ( $record->getRMSubType() ne "" ) {
        print "/" . $record->getRMSubType();
      }
      print " ";
      print $record->getDescription();
      print "\n";
      $seq =~ s/(\S{50})/$1\n/g;
      $seq .= "\n"
          unless ( $seq =~ /.*\n+$/s );
      print $seq;
    }
  }
  elsif ( $fileFormat eq "dfam" ) {
    my $record = $db->getRecord( $i );
    if ( defined $options{'stat'} ) {
      my $seqLen = 0;
      print ">"
          . $record->getId() . "#"
          . $record->getRMType() . "/"
          . $record->getRMSubType()
          . "  Length = "
          . "$seqLen bp";

      foreach my $name ( $record->getRMSpeciesArray() ) {
        print " \@$name";
      }
      print "\n";

      $totalLength += $seqLen;
      $ancestLen += $seqLen if ( $isAncestor == 1 && !$isDescendant == 1 );
      $cladeLen += $seqLen if ( $isDescendant == 1 );
      ++$totalNumber;
    }
    else {
      print "" . $record->getRecordLines();
    }
  }
}
print "\n";

if ( defined $options{'stat'} ) {
  if ( defined $options{'species'} ) {
    if ( !$options{'clade'} ) {
      print STDERR "$ancestCnt ancestral and ubiquitous sequence(s) with a "
          . "total length of "
          . $ancestLen . " bp\n";
    }
    print STDERR
"$equalCnt $realSpeciesName specific repeats with a total length of $equalLen bp\n";
    print STDERR
        "$cladeCnt lineage specific sequence(s) with a total length of "
        . $cladeLen . " bp\n";
    print STDERR "-" x ( 80 ) . "\n";
  }
  else {
    print STDERR "$totalNumber sequence(s) with a total length of "
        . $totalLength . " bp\n";
    print STDERR ""
        . scalar( keys( %uniqSpeciesLabels ) )
        . " species labels defined\n";
  }
}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my _displayLibraryTaxonomy( library_ref => value,
##                                  taxonomy_ref => value );
##
##      library_ref       : A reference to a RmEMBL object
##      taxonomy_ref      : A reference to a Taxonomy object
##
##  Returns
##     Prints the species tree for all repeat species in the
##     RepeatMasker library.
##
##-------------------------------------------------------------------------##
sub _displayLibraryTaxonomy {
  my %parameters = @_;

  print ""
      . ( &caller( 0 ) )[ 0 ] . "::"
      . ( &caller( 0 ) )[ 3 ] . "( "
      . @{ [ %parameters ] }
      . "): Called\n"
      if ( $DEBUG );

  my $libRef = $parameters{'library_ref'};
  my $taxRef = $parameters{'taxonomy_ref'};

  my $seqCount      = $libRef->getRecordCount();
  my %uniqSpecies   = ();
  my %actualSpecies = ();
  for ( my $i = 0 ; $i < $seqCount ; $i++ ) {
    my $record = $libRef->getRecord( $i );
    foreach my $name ( $record->getRMSpeciesArray() ) {
      $name =~ s/_/ /g;
      $name = $tax->isSpecies( $name );
      $uniqSpecies{$name}++;
      $actualSpecies{$name}++;
      my @ancSpecs = $taxRef->getLineage( $name );
      foreach my $ancSpec ( @ancSpecs ) {
        next if ( $ancSpec eq $name );
        $uniqSpecies{$ancSpec}++;
      }
    }
  }
  print "" . $taxRef->getTree( "root", \%uniqSpecies, \%actualSpecies ) . "\n";

}

1;

