#!/u1/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) LibraryUtils.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Module to assist in the validation and updating of
##      RepeatMasker's growing set of library files.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2017-2017 Developed by
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
#     $Log: LibraryUtils.pm,v $
#     Revision 1.6  2017/02/01 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

LibraryUtils.pm - Validate and update RepeatMasker libaries

=head1 SYNOPSIS

use LibraryUtils;

Usage:

LibraryUtils::validate()

=head1 DESCRIPTION

  A set of subroutines to assist RepeatMasker with managing the installation
  of multiple repeat libraries.

The options are:

=head1 INSTANCE METHODS

=cut

#
# Module Dependence
#
package LibraryUtils;
use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Data::Dumper;
use RmEMBL;
use DFAM;
use File::Basename;

my $DEBUG = 0;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

##-------------------------------------------------------------------------##

=begin

=over 4

=item my $versionString = getLibraryVersionStr( $libFile );

Return the version from the header of the given library file.

=back

=cut

##-------------------------------------------------------------------------##
sub getLibraryVersionStr {
  my $libFile = shift;

  open INVER, "<$libFile"
      or die "getLibraryVersion(): Could not open file $libFile";
  my $releaseStr  = "";
  my $searchLimit = 50;    # 50 lines max
  while ( <INVER> ) {

    # CC   RELEASE 20170125;
    # #   Release: Dfam_2.0
    if (    /^#\s+Release:\s+(Dfam_\S+).*/
         || /^(?:CC|##)\s+RELEASE\s+(\S+);.*/ )
    {
      $releaseStr = $1;
      last;
    }
    last if ( $searchLimit-- == 0 );
  }
  close INVER;

  return $releaseStr;
}

sub getCombinedLibrarySources {
  my $libFile = shift;

  open INVER, "<$libFile"
      or die "getLibraryVersion(): Could not open file $libFile";
  my %sources     = ();
  my $searchLimit = 50;    # 50 lines max
  my $inSources   = 0;
  while ( <INVER> ) {
    $inSources = 1 if ( /^CC\s+Sources:.*/ );

    # CC   Dfam_Consensus RELEASE 20170125;
    # CC   RepBase RELEASE 20170125;
    if ( /^CC\s+(\S+)\s+RELEASE\s+(\S+);.*/ ) {
      $sources{$1} = $2;
    }
    last if ( $searchLimit-- == 0 || ( $inSources && /^CC\s*$/ ) );
  }
  close INVER;

  return \%sources;
}

sub validateLibraries {

  # Interrogate the Libraries subdirectory
  my $libDir  = shift;
  my $libType = shift;

# Old filestruce ( pre-2017 )
#   /Libraries
#     RepeatMaskerLib.embl         : RepBase RepeatMasker Edition
#                                    or
#                                    Minimal RepeatMasker Library containing only
#                                    non-RepBase repeats
#     Dfam.hmm                     : Dfam library ( no processing necessary )
#
# File Structure
#    /Libraries
#         RepeatMaskerLib.embl     : The combined libraries for RepeatMasker
#         DfamConsensus.embl       : The Dfam Consensus library
#                                    ( which now contains the Minimal RM Library )
#         RMRBMeta.embl            : The RepeatMasker Metadata for Repbase
#         RMRBSeqs.embl            : RepBase Library ( RepeatMasker Edition ) from GIRI
#         Dfam.hmm                 : Dfam library ( no processing necessary )
#

  my $mainLibrary          = "RepeatMaskerLib.embl";
  my $dfamConsensusLibrary = "DfamConsensus.embl";
  my $RBRMSeqLibrary       = "RMRBSeqs.embl";
  my $RBRMMetaLibrary      = "RMRBMeta.embl";
  my $dfamLibrary          = "Dfam.hmm";

  my %dbAlias = (
                  'Dfam_Consensus' => 'dc',
                  'RepBase'        => 'rb'
  );

  my $rmLibraryVersionKey  = "";
  my $isLibraryCombined    = 0;
  my $rmLibraryDescription = "";

  if ( $libType eq "HMM" ) {
    if ( -s "$libDir/$dfamLibrary" ) {
      my $dfamVersion = getLibraryVersionStr( "$libDir/$dfamLibrary" );
      $rmLibraryVersionKey  = $dfamVersion;
      $rmLibraryDescription = "Dfam database version $dfamVersion";
    }
    else {
      print
          "\n\nThe Dfam Profile HMM library ( $libDir/$dfamLibrary ) was not\n"
          . "found.  Please download and install the latest version from\n"
          . "http://www.dfam.org.\n\n";
      die;
    }
  }
  else {

    # Do we have a main consensus library file?
    if ( -s "$libDir/$mainLibrary" ) {

      # Which structure do we have?
      my $libSources = getCombinedLibrarySources( "$libDir/$mainLibrary" );
      if ( defined $libSources && keys( %{$libSources} ) ) {

        # Main library file is a combined library ( new structure )
        $isLibraryCombined = 1;

        # Validate RepBase RepeatMasker Edition
        my $mustRebuild = 0;
        if ( -s "$libDir/$RBRMSeqLibrary" ) {

          # Validate that we have included RepBase
          my $rbSeqVersion = getLibraryVersionStr( "$libDir/$RBRMSeqLibrary" );
          if ( $rbSeqVersion ne $libSources->{'RepBase'} ) {

            #
            my $rmRbMetaVersion;
            if ( -s "$libDir/$RBRMMetaLibrary" ) {
              $rmRbMetaVersion =
                  getLibraryVersionStr( "$libDir/$RBRMMetaLibrary" );
              $rmRbMetaVersion = undef if ( $rmRbMetaVersion ne $rbSeqVersion );
            }
            if ( !$rmRbMetaVersion ) {
              print
                  "\n\nThe Repbase RepeatMasker Edition database has changed\n"
                  . "( RELEASE = $rbSeqVersion ), however the corresponding\n"
                  . "metadata library file ( $libDir/$RBRMMetaLibrary ) is missing or\n"
                  . "out of date.  Please obtain the metadata library from:\n"
                  . "http:/www.repeatmasker.org/libraries/RepeatMaskerMetaData-$rbSeqVersion.tar.gz\n"
                  . "Once this file is in placed in $libDir rerun RepeatMasker to continue.\n\n";
              die;
            }
            print
"RepBase RepeatMasker Edition database changed ( RELEASE = $rbSeqVersion ).\n";
            $mustRebuild = 1;
          }
        }

        # Validate Dfam_consensus
        if ( -s "$libDir/$dfamConsensusLibrary" ) {

          # Validate that we have included Dfam_consensus
          my $dfamConsVersion =
              getLibraryVersionStr( "$libDir/$dfamConsensusLibrary" );
          if ( $dfamConsVersion ne $libSources->{'Dfam_Consensus'} ) {
            print
"Dfam_consensus database changed ( RELEASE = $dfamConsVersion ).\n";
            $mustRebuild = 1;
          }
        }
        if ( $mustRebuild ) {
          rebuildMainLibrary( $libDir );
          $libSources = getCombinedLibrarySources( "$libDir/$mainLibrary" );
        }

        $rmLibraryVersionKey = join( "-",
                                     map { $dbAlias{$_} . $libSources->{$_} }
                                         sort keys( %{$libSources} ) );
        $rmLibraryDescription = "RepeatMasker Combined Database: "
            . join( ", ",
                    map { $_ . "-" . $libSources->{$_} }
                        sort keys( %{$libSources} ) );
      }
      else {

        # LEGACY SUPPORT
        # Main library is either a minimal library or a older
        # RepBase RepeatMasker Edition file
        my $isMinimum    = 0;
        my $rmlibVersion = getLibraryVersionStr( "$libDir/$mainLibrary" );
        if ( $rmlibVersion =~ /(\d+)-min/ ) {
          $rmlibVersion = $1;
          $isMinimum    = 1;
        }
        if ( $rmlibVersion && $rmlibVersion <= 20160829 ) {

          # Last valid version for legacy structure
          # Check to see if we have newer files in the directory and
          # warn that they are not being used.
          my @extraLibs;
          push @extraLibs, $dfamConsensusLibrary
              if ( -s "$libDir/$dfamConsensusLibrary" );
          push @extraLibs, $RBRMSeqLibrary if ( -s "$libDir/$RBRMSeqLibrary" );
          if ( @extraLibs ) {
            print "\n\nNewer libraries exist in $libDir: "
                . join( ", ", @extraLibs ) . "\n"
                . "but are not configured for use by RepeatMasker.  To enable them, remove\n"
                . "the $libDir/$mainLibrary file and rerun RepeatMasker to rebuild it.\n\n";
            die;
          }
          print
              "Legacy format: rmlibVersion = $rmlibVersion.  Ok to continue\n";
        }
        else {
          die "Legacy file format for $libDir/$mainLibrary yet\n"
              . "version ( $rmlibVersion ) is not valid.";
        }
        $rmLibraryVersionKey  = $rmlibVersion;
        $rmLibraryDescription =
            "RepBase/RepeatMasker database version $rmlibVersion";
        $rmLibraryDescription .= "-min" if ( $isMinimum );

    # Used to be
    # "RepBase Update $rmLibraryVersion, RM database version $rmLibraryVersion";
      }
    }
    else {

      # We don't have a main library file
      if ( -s "$libDir/$dfamConsensusLibrary" || -s "$libDir/$RBRMSeqLibrary" )
      {
        $isLibraryCombined = 1;
        rebuildMainLibrary( $libDir );
        my $libSources = getCombinedLibrarySources( "$libDir/$mainLibrary" );
        $rmLibraryVersionKey = join( "-",
                                     map { $dbAlias{$_} . $libSources->{$_} }
                                         sort keys( %{$libSources} ) );
        $rmLibraryDescription = "RepeatMasker Combined Database: "
            . join( ", ",
                    map { $_ . "-" . $libSources->{$_} }
                        sort keys( %{$libSources} ) );

      }
      else {
        print
            "\n\nNo repeat libraries found!  At a minimum the Dfam_consensus\n"
            . "is required to run.  Please download and install the latest \n"
            . "Dfam_consensus.  It is highly recommended that you also install the\n"
            . "latest RepBase RepeatMasker Edition library obtainable from GIRI.\n"
            . "General instructions can be found here: http://www.repeatmasker.org\n\n";
        die;
      }
    }
  }
  return ( $isLibraryCombined, $rmLibraryVersionKey, $rmLibraryDescription );
}

sub rebuildMainLibrary {
  my $libDir = shift;

  my $mainLibrary          = "RepeatMaskerLib.embl";
  my $dfamConsensusLibrary = "DfamConsensus.embl";
  my $RBRMSeqLibrary       = "RMRBSeqs.embl";
  my $RBRMMetaLibrary      = "RMRBMeta.embl";

  print "Rebuilding $mainLibrary library\n";

  # Backup old library ( only one backup kept )
  unlink( "$libDir/$mainLibrary.old" )
      if ( -s "$libDir/$mainLibrary.old" );
  rename( "$libDir/$mainLibrary", "$libDir/$mainLibrary.old" )
      if ( -s "$libDir/$mainLibrary" );

  my $headerSources = "";

  my $combinedDb = new RmEMBL();

  if ( -s "$libDir/$dfamConsensusLibrary" ) {
    my $savBuf = $|;
    $| = 1;
    print "    Reading Dfam_consensus database...";
    $| = $savBuf;
    my $db = RmEMBL->new( fileName => "$libDir/$dfamConsensusLibrary" );
    print "\r  - Read in "
        . $db->size()
        . " sequences from $libDir/$dfamConsensusLibrary\n";
    $combinedDb->addAll( $db );
    my $libVersion = getLibraryVersionStr( "$libDir/$dfamConsensusLibrary" );
    $headerSources .=
"CC    Dfam_Consensus RELEASE $libVersion;                            *\n";
  }

  if ( -s "$libDir/$RBRMSeqLibrary" ) {
    my $savBuf = $|;
    $| = 1;
    print "    Reading RepBase RepeatMasker Edition database...";
    $| = $savBuf;
    my $seqs = RmEMBL->new( fileName => "$libDir/$RBRMSeqLibrary" );
    my %seqId = ();
    for ( my $i = 0 ; $i < $seqs->size() ; $i++ ) {
      my $rec = $seqs->get( $i );
      $seqId{ $rec->getId() } = $rec;
    }
    print "\r  - Read in "
        . $seqs->size()
        . " sequences from $libDir/$RBRMSeqLibrary\n";
    undef $seqs;

    my $libVersion = getLibraryVersionStr( "$libDir/$RBRMSeqLibrary" );
    $headerSources .=
"CC    RepBase RELEASE $libVersion;                                   *";

    my $savBuf = $|;
    $| = 1;
    print "    Reading metadata database...";
    $| = $savBuf;
    my $meta = RmEMBL->new( fileName => "$libDir/$RBRMMetaLibrary" );
    for ( my $i = 0 ; $i < $meta->size() ; $i++ ) {
      my $mRec = $meta->get( $i );
      my $sRec = $seqId{ $mRec->getId() };
      if ( $mRec && $sRec ) {
        $mRec->setLength( $sRec->getLength() );
        $mRec->setSequence( $sRec->getSequence() );
        $mRec->setComposition( 'A', $sRec->getCompositionElement( 'A' ) );
        $mRec->setComposition( 'C', $sRec->getCompositionElement( 'C' ) );
        $mRec->setComposition( 'G', $sRec->getCompositionElement( 'G' ) );
        $mRec->setComposition( 'T', $sRec->getCompositionElement( 'T' ) );
        $mRec->setComposition( 'other',
                               $sRec->getCompositionElement( 'other' ) );
        $mRec->pushComments( "Source: RepBase RepeatMasker Edition\n" );
      }
      else {
        print "Error! Could not find " . $mRec->getId() . "\n";
      }
    }
    print "\r  - Read in "
        . $meta->size()
        . " annotations from $libDir/$RBRMMetaLibrary\n";
    $combinedDb->addAll( $meta );
  }

  my $savBuf = $|;
  $| = 1;
  print "  Saving $mainLibrary library...";
  $| = $savBuf;
  my $headerStr =
      "CC ****************************************************************
CC                                                                *
CC   RepeatMasker Combined Library                                *
CC    This is a merged file of external library sources.          *
CC    See the original libraries for detailed copyright           *
CC    and licensing restrictions.                                 *
CC                                                                *
CC   Sources:                                                     *
$headerSources
CC                                                                *
CC   RepeatMasker software and database development and           *
CC   maintenance are currently funded by an NIH/NHGRI             *
CC   R01 grant HG02939-01 to Arian Smit.  RepBase Update          *
CC   development and maintenance are funded by NIH/NLM grant      *
CC   No.2P41LM006252-07A1 to Jerzy Jurka.                         *
CC                                                                *
CC ****************************************************************";

  $combinedDb->writeEMBLFile( "$libDir/$mainLibrary", $headerStr );
  print "\r$mainLibrary: " . $combinedDb->size() . " total sequences.\n";

}

1;
