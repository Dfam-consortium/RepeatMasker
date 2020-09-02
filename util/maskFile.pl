#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) maskRepeats.pl
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
#     $Log: maskFile.pl,v $
#     Revision 1.39  2017/02/01 21:01:57  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

maskRepeats.pl - A utility script to mask a sequence file given a *.out file

=head1 SYNOPSIS

  maskRepeats.pl    -fasta <filename.fa> -annotations <filename.out> 
                    [-softmask]
                    [ -minDivergence = <#> ] [ -maxDivergence = <#> ]

=head1 DESCRIPTION


The options are:

=over 4

=back

=head1 SEE ALSO

ReapeatMasker

=head1 COPYRIGHT

Copyright 2006 Robert Hubley, Institute for Systems Biology

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
use CrossmatchSearchEngine;
use FastaDB;
use File::Basename;

my $DEBUG = 0;

#
# Version
#
my $Version = 0.1;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',           # print out the version and exit
                    '-fasta=s',
                    '-softmask',
                    '-annotations=s',
                    '-minDivergence=s',
                    '-maxDivergence=s'
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

usage()
    if ( !( exists $options{'fasta'} && exists $options{'annotations'} ) );

my $db;
if ( -e $options{'fasta'} ) {
  $db = FastaDB->new(
                      fileName    => $options{'fasta'},
                      openMode    => SeqDBI::ReadOnly,
                      maxIDLength => 80
  );
}
else {
  die "Could not find file " . $options{'fasta'} . "\n";
}

if ( $options{'softmask'} ) {
  &maskSequence(
               maskFormat     => 'xsmall',
               seqDB          => $db,
               annotationFile => $options{'annotations'},
               outputFile     => $options{'fasta'} . ".masked",
               minDivergence  => $options{'minDivergence'},
               maxDivergence  => $options{'maxDivergence'}
  );
}else {
  &maskSequence(
               seqDB          => $db,
               annotationFile => $options{'annotations'},
               outputFile     => $options{'fasta'} . ".masked",
               minDivergence  => $options{'minDivergence'},
               maxDivergence  => $options{'maxDivergence'}
  );
}

##-------------------------------------------------------------------------##
## Use:  my ( $seq_cnt, $totalSeqLen, $nonMaskedSeqLen, $totGCLevel,
##             $totBPMasked ) =
##                    &maskSequence ( $seqDB, $annotationFile,
##                                    $outputFile );
##  Returns
##
##     $seq_cnt:          The number of sequences in the FASTA file.
##     $totalSeqLen:      The absoulte length of all sequences combined.
##     $nonMaskedSeqLen:  Length of sequence (excluding runs of >20 N's
##                         and X's) of the pre-masked sequence.
##     $totGCLevel:       The GC content of the original sequence.
##     $totBPMasked:      The total bp we masked
##
##-------------------------------------------------------------------------##
sub maskSequence {
  my %parameters = @_;

  my $maskFormat     = $parameters{'maskFormat'};
  my $seqDB          = $parameters{'seqDB'};
  my $annotationFile = $parameters{'annotationFile'};
  my $outputFile     = $parameters{'outputFile'};
  my $minDivergence  = $parameters{'minDivergence'};
  my $maxDivergence  = $parameters{'maxDivergence'};

  print "maskSequence()\n" if ( $DEBUG );

  my %annots = ();

  #
  # Open up a search results object
  #
  my $searchResults =
      CrossmatchSearchEngine::parseOutput( searchOutput => $annotationFile );

  #
  # Read in annotations and throw away the rest
  #
  my $prevResult;
  for ( my $i = 0 ; $i < $searchResults->size() ; $i++ ) {
    my $result = $searchResults->get( $i );
    my $start  = $result->getQueryStart();
    my $end    = $result->getQueryEnd();
    if (    defined $prevResult
         && $prevResult->getQueryName() eq $result->getQueryName()
         && $prevResult->getQueryEnd() >= $start )
    {
      next if ( $prevResult->getQueryEnd() >= $end );
      $start = $prevResult->getQueryEnd() + 1;
    }
    next
        if ( defined $minDivergence
             && $result->getPctDiverge() < $minDivergence );
    next
        if ( defined $maxDivergence
             && $result->getPctDiverge() > $maxDivergence );
    push @{ $annots{ $result->getQueryName() } },
        {
          'begin' => $start,
          'end'   => $end
        };
    $prevResult = $result;
  }
  undef $searchResults;

  my @seqIDs     = $seqDB->getIDs();
  my $seq_cnt    = scalar( @seqIDs );
  my $sublength  = $seqDB->getSubtLength();
  my $totGCLevel = 100 * $seqDB->getGCLength() / $sublength;
  $totGCLevel = sprintf "%4.2f", $totGCLevel;
  my $totalSeqLen     = 0;
  my $totBPMasked     = 0;
  my $nonMaskedSeqLen = 0;
  my $workseq         = "";
  open OUTFILE, ">$outputFile";

  foreach my $seqID ( @seqIDs ) {
    my $seq = $seqDB->getSequence( $seqID );
    $totalSeqLen += length $seq;
    $workseq = $seq;
    $nonMaskedSeqLen += length $workseq;
    while ( $workseq =~ /([X,N]{20,})/ig ) {
      $nonMaskedSeqLen -= length( $1 );
    }
    foreach my $posRec ( @{ $annots{$seqID} } ) {
      my $beginPos = $posRec->{'begin'};
      my $endPos   = $posRec->{'end'};
      my $repLen   = $endPos - $beginPos + 1;
      substr( $workseq, $beginPos - 1, $repLen ) = "0" x ( $repLen );
      if ( $maskFormat eq 'xsmall' ) {
        substr( $seq, $beginPos - 1, $repLen ) =
            lc( substr( $seq, $beginPos - 1, $repLen ) );
      }
      elsif ( $maskFormat eq 'x' ) {
        substr( $seq, $beginPos - 1, $repLen ) = "X" x ( $repLen );
      }
      else {
        substr( $seq, $beginPos - 1, $repLen ) = "N" x ( $repLen );
      }
      $totBPMasked += $repLen;
    }
    print OUTFILE ">" . $seqID;
    my $desc = $seqDB->getDescription( $seqID );
    if ( $desc ne "" ) {
      print OUTFILE " " . $desc;
    }
    print OUTFILE "\n";
    $seq =~ s/(\S{50})/$1\n/g;
    $seq .= "\n"
        unless ( $seq =~ /.*\n+$/s );
    print OUTFILE $seq;
  }
  close OUTFILE;

  return ( $seq_cnt, $totalSeqLen, $nonMaskedSeqLen, $totGCLevel,
           $totBPMasked );
}
