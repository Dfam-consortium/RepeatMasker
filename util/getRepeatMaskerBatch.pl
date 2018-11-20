#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) getRepeatMaskerBatch.pl
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      Using RepeatMasker's batching mechanism generate a particular
##      numbered batch from the input file.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2002-2011 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
#  ChangeLog:
#
#    $Log$
#
###############################################################################
#
# To Do:
#
#

=head1 NAME

getRepeatMaskerBatch.pl - Retrieve a RM batch given a fasta file

=head1 SYNOPSIS

  getRepeatMaskerBatch.pl [-options] <seqfiles(s) in fasta format>

=head1 DESCRIPTION

The options are:

=over 4

=item -frag [number] 

Maximum sequence length masked without fragmenting (default 60000,
300000 for DeCypher)

=item -batch [number]  

Retrieve only batch #number.

=back

=head1 SEE ALSO

=over 4

RepeatMasker

=back

=head1 COPYRIGHT

Copyright 2007-2011 Arian Smit, Institute for Systems Biology

=head1 AUTHORS

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use Getopt::Long;
use POSIX qw(:sys_wait_h);
use File::Copy;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd;

# RepeatMasker Libraries
use FastaDB;
use SimpleBatcher;

# A bug in 5.8 produces way too many warnings
if ( $] && $] >= 5.008003 ) {
  use warnings;
}

# Debugging flag
my $DEBUG = 0;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co RepeatMasker" the
#  $Name:  $ macro will be blank and thus we use
#  this as the development version.
#
my $CVSTag = '$Name:  $';
my $version;
if ( $CVSTag =~ /\$\s*Name:\s*open-(\S+)\s*\$/ ) {
  $version = $1;
  $version =~ s/-/./g;
  $version = "open-$version";
}
else {
  $version =
'development-$Id: getRepeatMaskerBatch.pl,v 1.16 2017/02/01 21:01:57 rhubley Exp $';
}

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @opts = qw( frag=s batch=s );

#
# Get the supplied command line options, and set flags
#
my %options = ();
unless ( &GetOptions( \%options, @opts ) ) {
  exec "pod2text $0";
  exit( 0 );
}

# Print the internal POD documentation if something is missing
if ( $#ARGV == -1 && !$options{'help'} ) {
  print "No query sequence file indicated\n\n";

  # This is a nifty trick so we don't have to have
  # a duplicate "USAGE()" subroutine.  Instead we
  # just recycle our POD docs.  See PERL POD for more
  # details.
  exec "pod2text $0";
  die;
}

##
## RepeatMasker Batching Parameters
##
# Most computers can handle the memory/cpu requirements of
# running 60kb chunks.  Currently our batching mechanism will
# not breakup anything smaller than this.
my $fragmentSize = 60000;

#
# Selection of a batch overlap length can have large impacts on the program.
# The overlap boundaries are places where edge effects produce partial
# overlapping annotations.  Also matrix differences in flanking batches
# can cause the same repeat to have different divergence, score and alignment
# length characteristics.
#
my $overlapLen = 2000;

#
# User supplied fragment length
#
if ( defined $options{'frag'} ) {
  if ( $options{'frag'} < ( 2 * $overlapLen ) ) {
    warn "RepeatMasker: You may not use a fragment size (-frag "
        . "$options{'frag'} ) which is less than 2 times the overlap "
        . "length (overlapLen = $overlapLen).  Defaulting to $fragmentSize\n";
  }
  else {
    $fragmentSize = $options{'frag'};
  }
}

#
# Parse filenames
#
foreach my $file ( @ARGV ) {
  if ( $file =~ /\s/ ) {
    die "RepeatMasker can not handle filenames with spaces "
        . "like the file \"$file\"\n";
  }
  elsif ( $file =~ /([\`\!\$\^\&\*\(\)\{\}\[\]\|\\\;\"\'\<\>\?])/ ) {
    die "RepeatMasker can not handle filenames with the special "
        . "character \"$1\" as in the file \"$file\"\n";
  }
}

#
# Main loop
#
FILECYCLE:
foreach my $file ( @ARGV ) {

  unless ( -r $file ) {
    print "cannot read file $file in " . cwd() . "\n";
    next;
  }

  unless ( -s $file ) {
    print "File $file appears to be empty.\n";
    next;
  }

  my $compressed = "";
  if ( $file =~ /\.gz$/ ) {
    unless ( `gunzip $file 2>&1` ) {

      # Name $file only changes if gunzip did not complain
      # (file may end with .gz but not be zipped
      $file =~ s/\.gz$//;
      $compressed = "zipped";
    }
  }
  elsif ( $file =~ /\.Z$/ ) {
    unless ( `uncompress $file 2>&1` ) {
      $file =~ s/\.Z$//;
      $compressed = "Zed";
    }
  }

  # With one file $#ARGV == 0
  print "\nanalyzing file $file\n" if ( $#ARGV >= 0 );

  ## Create a batcher object.  Upon construction the
  ## object will survey the fasta file, check for syntax
  ## errors, and create a byte index for all parseable sequences
  ## in the file.  Copy the file so we don't mess with the
  ## original.
  my $db = FastaDB->new(
                         fileName    => $file,
                         openMode    => SeqDBI::ReadWrite,
                         maxIDLength => 50
  );
  my $batcher = SimpleBatcher->new( $db, $fragmentSize, $overlapLen );

  ## How many sequences are in the fasta file?
  my $totseqcnt = $db->getSeqCount();

  ## Don't process this file unless we have some
  ## unambiguous (ACGT) sequence.
  my $sublength = 0;
  unless ( $sublength = $db->getSubtLength() ) {
    print "INFO: File $file contains only ambiguous sequence. Skipping.";
    next FILECYCLE;
  }

  my $totGClevel = 100 * $db->getGCLength() / $sublength;
  $totGClevel = sprintf "%4.2f", $totGClevel;

  my $totseqlen = $db->getSeqLength();

  my $batchCount = $batcher->getBatchCount();
  my $numberChildren = $options{'parallel'} ? $options{'parallel'} : 1;
  $numberChildren = $batchCount if ( $batchCount < $numberChildren );

  my $badForkCount = 0;     # A failsafe for in case the fork goes badly
  my $badForkMax   = 20;    # The number of bad forks ( in a row ) before exit
  my $retryLimit   = 2;     # Attempt to run a batch 2 times before failing it
  my %children     = ();    # A hash of children PIDs with their batch nums
  my $child_id     = 0;     # The PID of returned by fork
  my %batchStatus  = ();    # A hash which holds the retry count for
                            #  each batch number.

  #
  # Job processing loop
  #   Process all batches stored in the batchStatus hash.  If
  #   for some reason a job fails the entry in $batchStatus is
  #   incremented.  We continue to re-run this batch until the
  #   $retryLimit is reached.
  #
  if ( $options{'batch'} ) {
    die "Batch number outside range ( 1 - $batchCount ): $options{'batch'}\n"
        if ( $options{'batch'} > $batchCount || $options{'batch'} < 1 );

    my $k = $options{'batch'};

    my $batchSeqFile = $file . "_batch-" . $k . ".fa";

    ## Get batch parameters
    my $seq_cnt = $batcher->getBatchSeqCount( $k );
    my $seqlen  = $batcher->getBatchSeqLength( $k );
    my $frac_GC = $batcher->getBatchAverageGC( $k );
    print "Creating Batch " . $k
        . " seq count = $seq_cnt "
        . "len = $seqlen, average GC = $frac_GC\n";

    ## Create the batch file
    $batcher->writeBatchFile( $k, $batchSeqFile );

  }
  else {
    for ( my $k = 1 ; $k <= $batchCount ; $k++ ) {
      my $batchSeqFile = $file . "_batch-" . $k . ".fa";

      ## Get batch parameters
      my $seq_cnt = $batcher->getBatchSeqCount( $k );
      my $seqlen  = $batcher->getBatchSeqLength( $k );
      my $frac_GC = $batcher->getBatchAverageGC( $k );
      print "Creating Batch " . $k
          . " seq count = $seq_cnt "
          . "len = $seqlen, average GC = $frac_GC\n";

      ## Create the batch file
      $batcher->writeBatchFile( $k, $batchSeqFile );
    }    # End for ( my $k = 0 ; $k < $numberToStart...
  }
}    # END FILECYCLE

# We are soooo done
exit( 0 );

