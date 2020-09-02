#!/u1/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) addRepBase.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Jeb Rosen          jeb.rosen@isbscience.org
##  Description:
##      Script to merge the RepBase RepeatMasker Edition
##      libraries into the main RepeatMasker library.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2017-2020 Developed by
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
#     $Log$
#
###############################################################################
#
# To Do:
#

=head1 NAME

addRepBase.pl - Merge RepBase RepeatMasker Edition into the RepeatMasker library

=head1 SYNOPSIS

Usage:

    addRepBase.pl

=head1 DESCRIPTION

  Merges RepBase RepeatMasker Edition (RMRBMeta.embl and RMRBSeqs.embl) with
  Dfam.h5, and makes it the main RepeatMasker library (RepeatMaskerLib.h5)

The options are:

=over 4

=item -libdir [path_to_library_directory]

Use an alternate library directory to for the primary repeat libraries.
These include the Dfam.h5 and the RBRM (Repbase RepeatMasker Edition )
distribution files. Normally these are stored/updated in the "Libraries/"
subdirectory of the main program which RepeatMasker will use by default
This parameter should only be used when it's not possible to keep the
libraries in the same place as the program.

=back


=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use EMBL;
use File::Copy;

my @getopt_args = ( '-libdir=s' );

my %options = ();
unless ( &GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0\n";
  exec "pod2text $0";
  exit( 1 );
}

# Allow the user to override the default library directory ( currently only by the command line )
my $LIBDIR = "$FindBin::Bin/Libraries";
if ( $options{'libdir'} ) {
  $LIBDIR = $options{'libdir'};
  if ( ! -d $LIBDIR ) {
    die "The specified library directory $options{'libdir'} does not exist!\n";
  }
}


my $mainLibrary      = "RepeatMaskerLib.h5";
my $dfamLibrary      = "Dfam.h5";
my $RMRBLibrary      = "RMRB.embl";

print "Rebuilding $mainLibrary master library\n";

my $versionWarned = 0;
if ( ! -s "$LIBDIR/$RMRBLibrary") {
  combineRMRBMetaWithSeqs( $LIBDIR );
}

my $libVersion = getLibraryVersionStr( "$LIBDIR/$RMRBLibrary" );
if ( !$versionWarned && $libVersion ne "20181026" ) {
  $versionWarned = 1;
  print "WARNING: RepBase RepeatMasker edition older than 20181026 is untested with this version of RepeatMasker.\n";
}

my $savBuf = $|;
$| = 1;
print "  Merging Dfam + RepBase into $mainLibrary library...";
$| = $savBuf;

copy( "$LIBDIR/$dfamLibrary", "$LIBDIR/$mainLibrary.writing" )
  or die "Failed to copy $LIBDIR/$dfamLibrary to $LIBDIR/$mainLibrary.writing.\n" .
         "Is the source file missing, or is your user missing write permissions to the directory?\n";

my $REPEATMASKER_DIR = "$FindBin::Bin";
my $FAMDB = "$REPEATMASKER_DIR/famdb.py";

system("$FAMDB -i $LIBDIR/$mainLibrary.writing append $LIBDIR/$RMRBLibrary --name 'Dfam withRBRM' --description 'RBRM - RepBase RepeatMasker Edition - version $libVersion'");
my $status = $?;
if ( $status ) {
  die "Failed to append $LIBDIR/$RMRBLibrary to $LIBDIR/$mainLibrary.writing.\n" .
      "Process exited with code " . ( $status >> 8 ) . ".\n";
}

# Backup old library ( only one backup kept )
if ( -s "$LIBDIR/$mainLibrary" ) {
  unlink( "$LIBDIR/$mainLibrary.old" )
    if ( -s "$LIBDIR/$mainLibrary.old" );
  rename( "$LIBDIR/$mainLibrary", "$LIBDIR/$mainLibrary.old" )
}

# rename temporary file
rename( "$LIBDIR/$mainLibrary.writing", "$LIBDIR/$mainLibrary" );

print "\r\n\n";

system("$FAMDB -i $LIBDIR/$mainLibrary info");

# my $versionString = getLibraryVersionStr( $libFile );
# Return the version from the header of the given library file.
sub getLibraryVersionStr {
  my $libFile = shift;

  open INVER, "<$libFile"
      or die "getLibraryVersion(): Could not open file $libFile";
  my $releaseStr  = "";
  my $searchLimit = 50;    # 50 lines max
  while ( <INVER> ) {

    # RepBase RepeatMasker Edition
    #    CC   RELEASE 20170125;
    if ( /^(?:CC|##)\s+RELEASE\s+(\S+);.*/ )
    {
      $releaseStr = $1;
      last;
    }
    last if ( $searchLimit-- == 0 );
  }
  close INVER;

  return $releaseStr;
}

sub combineRMRBMetaWithSeqs {
  my $LIBDIR = shift;

  my $RMRBSeqLibrary   = "RMRBSeqs.embl";
  my $RMRBMetaLibrary  = "RMRBMeta.embl";
  my $RMRBLibrary      = "RMRB.embl";

  my $RMRBDb     = new EMBL();

  if ( -s "$LIBDIR/$RMRBSeqLibrary" ) {
    my $savBuf = $|;
    $| = 1;
    print "    Reading RepBase RepeatMasker Edition database...";
    my $libVersion = getLibraryVersionStr( "$LIBDIR/$RMRBSeqLibrary" );
    if ( $libVersion ne "20181026" ) {
      $versionWarned = 1;
      print "\nWARNING: RepBase RepeatMasker edition older than 20181026 is untested with this version of RepeatMasker.\n";
    }
    $| = $savBuf;

    my $seqs = EMBL->new( fileName => "$LIBDIR/$RMRBSeqLibrary" );
    my %seqId = ();
    for ( my $i = 0 ; $i < $seqs->size() ; $i++ ) {
      my $rec = $seqs->get( $i );
      $seqId{ $rec->getId() } = $rec;
    }
    print "\r  - Read in "
        . $seqs->size()
        . " sequences from $LIBDIR/$RMRBSeqLibrary\n";
    undef $seqs;

    my $savBuf = $|;
    $| = 1;
    print "    Reading metadata database...";
    $| = $savBuf;
    my $meta = EMBL->new( fileName => "$LIBDIR/$RMRBMetaLibrary" );
    my $numMerged = 0;
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

        # Save combined RMRB.embl file
        $RMRBDb->add( $mRec );
      }
      else {
        print "Error! Could not find " . $mRec->getId() . "\n";
      }
    }

    # Save out complete RepBase RepeatMasker Edition
    my $headerStr =
        "CC ****************************************************************
CC                                                                *
CC   RepBase RepeatMasker Edition                                 *
CC    Please refer to GIRI (https://www.girinst.org/) for         *
CC    detailed copyright and licensing restrictions.              *
CC                                                                *
CC    RELEASE $libVersion;                                   *
CC
CC   RepeatMasker software, and maintenance are currently         *
CC   funded by an NIH/NHGRI R01 grant HG02939-01 to Arian Smit.   *
CC                                                                *
CC ****************************************************************";
    $RMRBDb->writeEMBLFile( "$LIBDIR/$RMRBLibrary", $headerStr );

    print "\r  - Read in "
        . $meta->size()
        . " annotations from $LIBDIR/$RMRBMetaLibrary\n";
  } else {
    die "Error! Could not find the RepBase RepeatMasker Edition sequences file,\n" .
        "$RMRBSeqLibrary. Please download RepBase RepeatMasker Edition and\n" .
        "extract it to the Libraries directory.\n";
  }
}
