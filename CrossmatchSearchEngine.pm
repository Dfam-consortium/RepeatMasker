#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) CrossmatchSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      the cross_match search engine.
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
#      'CrossmatchSearchEngine' );
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

CrossmatchSearchEngine

=head1 SYNOPSIS

use CrossmatchSearchEngine

  my $cEngine = CrossmatchSearchEngine->new( 
                       pathToEngine=>"/usr/local/bin/crossmatch" );

  $cEngine->setMatrix( "/users/bob/simple.matrix" );
  $cEngine->setQuery( "/users/bob/query.fasta" );
  $cEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $cEngine->search();

Usage: 

=head1 DESCRIPTION

An encapsulation class for the CrossMatch search engine.  It is capable
for running searches and returning the results as a SearchResultCollection.

The CrossMatch parser built into this object captures several types of i
data from the voluminous cross_match output stream.  The first is 
the high-scoring pair line or "hit" line:

Forward Strand:

 SW    perc perc perc qry   qry   qry  qry   subj           subj  subj subj 
 score div. del. ins. seq   begin end (left) seq            begin end (left) OV 
 ------------------------------------------------------------------------------
 2334  8.44 0.00 3.25 Human 127   737 (8222) AluSx#SINE/Alu  1    298 (14)   *  

Reverse Strand:

 SW    perc perc perc qry   qry   qry  qry     subj       subj  subj subj 
 score div. del. ins. seq   begin end (left) C seq       (left) end  begin OV 
 ------------------------------------------------------------------------------
 2334  8.44 0.00 3.25 Human 127   737 (8222) C AluSx#SINE (14)  298  1       


  SW score    = smith-waterman score of the match (complexity-adjusted, 
                by default).
  perc div.   = %substitutions in matching region.
  perc del.   = %deletions (in query seq rel to subject) in matching region.
  perc ins.   = %insertions (in query seq rel to subject) in matching region.
  qry seq     = id of query sequence.
  qry begin   = starting position of match in query sequence.
  qry end     = ending position of match in query sequence.
  qry (left)  = no. of bases in query sequence past the ending position of 
                match (so 0 means that the match extended all the way to 
                the end of the query sequence).
  C           = match is with the Complement of subject sequence.
  subj seq    = id of the subject sequence.
  subj (left) = The remaining bases in (complement of) subject sequence 
                prior to beginning of the match.
  subj end    = starting position of match in subject sequence (using 
                top-strand numbering).
  subj begin  = ending position of match in subject sequence.
  OV          = A "*" in this field indicates that there is a higher-scoring 
                match whose domain partly includes the domain of this match.

RepeatMasker has added two new fields to this format for it's standard
annotation output:

 SW    perc perc perc qry qry   qry  qry   subj subj  subj subj     Lineage
 score div. del. ins. seq begin end (left) seq  begin end (left) Id   Id    OV 
 ------------------------------------------------------------------------------

  Id         = A unique identifier for the repeat which corresponds
               to the repeat in the .align files.
  LineageId  = ?

The second type of data this object collects (optionally) is the
alignment data.  If the cross_match output contains alignment
data; typically in the form:

 NT_004321_1       1047 CACCCACATGCACACACACACGCGCGCACACACGCACACGCACACACATG 1096
                           v    ii           i i i       i     i        ii
 (CA)n#Simple_re      1 CACACACACACACACACACACACACACACACACACACACACACACACACA 50

 NT_004321_1       1097 CACACACGCGCAC--ACACGCACACATATGCACACACAAACGCACA 1142
                               i i         i      i ii        v  i
 (CA)n#Simple_re     51 CACACACACACACACACACACACACACACACACACACACACACACA 96

This object will parse these lines and include the
query and subject alignment sequences in the object.  NOTE: This can
greatly increase memory usage as the sequences are not compressed and
will include the gap characters ("-").

If the full cross_match output is available the matrix filename will be 
read from the cross_match invocation line.  If the input file is from 
RepeatMasker and contains the:

  Assumed background GC level in scoring matrices is 49 %

line, the GC attribute will be set.

Lastly, if cross_match was run with alignments turned on the 
transitions and transversion information will be read from
the following line:

  Transitions / transversions = 2.00 (42 / 21)


=head1 SEE ALSO

=over 4

SearchEngineI, SearchResultCollection

=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package CrossmatchSearchEngine;
use strict;
use POSIX qw(:sys_wait_h);
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use Carp;
use FileHandle;
use File::Basename;
use IPC::Open3;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "CrossmatchSearchEngine";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/cross_match\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak(
              "$CLASS" . "::set: Instance variable $name doesn't exist." . "" );
      }
      $this->$method( $value );
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setPathToEngine()

  Use: my $value    = getPathToEngine( );
  Use: my $oldValue = setPathToEngine( $value );

  Get/Set the fully qualified path to the search engine
  binary file.

=cut

##-------------------------------------------------------------------------##
sub getPathToEngine {
  my $this = shift;

  return $this->{'pathToEngine'};
}

sub setPathToEngine {
  my $this  = shift;
  my $value = shift;

  croak $CLASS. "::setPathToEngine( $value ): Program does not exist!"
      if ( not -x $value || `which $value` );

  my $result = `$value 2>&1`;
  if ( $result =~ /.*cross_match version ([\.\d]+).*/ ) {
    $this->{'version'} = $1;
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 getParameters()

  Use: my  $paramString  = getParameters( );

  Convert object parameters into cross_match command line parameters.

=cut

##-------------------------------------------------------------------------##
sub getParameters {
  my $this = shift;

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  # Generate parameter line
  my $parameters;
  my $value;
  if ( ( $value = $this->getScoreMode() ) ) {
    if ( $value == SearchEngineI::basicScoreMode ) {
      $parameters .= " -raw";
    }
  }
  if ( ( $value = $this->getWordRaw() ) ) {
    $parameters .= " -word_raw" if ( $value > 0 );
  }
  if ( ( $value = $this->getGenerateAlignments() ) ) {
    $parameters .= " -alignments" if ( $value > 0 );
  }
  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -gap_init $value";
  }
  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -ins_gap_ext $value";
  }
  if ( defined( $value = $this->getDelGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -del_gap_ext $value";
  }
  if ( ( $value = $this->getMinMatch() ) ) {
    $parameters .= " -minmatch $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMinScore() ) ) {
    $parameters .= " -minscore $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getBandwidth() ) ) {
    $parameters .= " -bandwidth $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMaskLevel() ) ) {
    $parameters .= " -masklevel $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMatrix() ) ) {

    # test if matrix exists
    if ( -f $value ) {
      $parameters .= " -matrix $value";
    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }
  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }
  if ( ( $value = $this->getSubject() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS. "::search: Error...subject ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  return ( "$engine $parameters" );
}

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

  Run the search and return a SearchResultsCollectionI.

=cut

##-------------------------------------------------------------------------##
sub search {
  my $this           = shift;
  my %nameValuePairs = @_;

  if ( %nameValuePairs ) {

    # TODO: Consider deprecating this way of setting up the object
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::search: Instance variable $name doesn't exist." );
      }
      $this->$method( $value );
    }
  }

  # Form the command line
  my $cmdLine = $this->getParameters();

  # Get matrix name for post-processing results
  my @matrixName = ();
  if ( ( my $value = $this->getMatrix() ) ) {

    # test if matrix exists
    if ( -f $value ) {
      my @path = split( /[\\\/]/, $value );
      my $matrix = $path[ $#path ];
      $matrix =~ s/[\n\r]//g;
      @matrixName = ( matrixName => $matrix );
    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  my $outputDirName;
  if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
    $outputDirName = $this->getTempDir();
  }
  else {
    $outputDirName = dirname( $this->getQuery() );
  }

  # Invoke engine and handle errors
  #print "Did run parameters=$parameters\n";
  my $PINPUT  = new FileHandle;
  my $POUTPUT = new FileHandle;
  my $PERROR  = new FileHandle;
  my $pid;

  #eval { $pid = open3( $PINPUT, $POUTPUT, $PERROR, "$engine $parameters" ) };
  #$pid = open3( $PINPUT, $POUTPUT, $PERROR, "$engine $parameters" );
  # NOTE: Here I wanted to use open3.  The problem with open3 is that
  #       programs writing output to both stderr and stdout will block
  #       if either buffer becomes full while the other is being read from.
  #       So....IO::Select needs to be used to alternate between reading
  #       from each one.

  print $CLASS
      . "::search(): Running crossmatch as:\n    "
      . "$cmdLine 2>/dev/null |\n"
      if ( $this->getDEBUG() );

  $pid = open( $POUTPUT, "$cmdLine 2>/dev/null |" );

  my $resultCode = 0;
  my $searchResultsCollection;

  # Create SearchResultCollection object from
  # the engine results.
  if ( $this->getDEBUG() ) {
    ## Create a debug file
    my $outFile;
    do {
      $outFile = $outputDirName . "/cmResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $resultCode = $? >> 8;
    $searchResultsCollection = parseOutput(
                                            searchOutput => $outFile,
                                            debug        => $this->getDEBUG(),
                                            @matrixName
    );
  }
  else {
    $searchResultsCollection = parseOutput(
                                            searchOutput => $POUTPUT,
                                            debug        => $this->getDEBUG(),
                                            @matrixName
    );
    close $POUTPUT;
    $resultCode = $? >> 8;
  }

  for ( my $i = 0 ; $i < $searchResultsCollection->size() ; $i++ ) {
    my $result = $searchResultsCollection->get( $i );

    #
    # Fix pctIns and pctDel values.  Traditionally crossmatch computes these
    # based on the query length.  We think it makes more sense to base pctDel
    # as a percentage of deletions in the aligned consensus and the pctIns
    # as the percentage of insertions in the aligned genomic sequence.
    # This fix only involves fixing the percDel
    #
    # NOTE: Rounding error can occur here
    #
    my $qryLen = $result->getQueryEnd() - $result->getQueryStart() + 1;
    my $sbjLen = $result->getSubjEnd() - $result->getSubjStart() + 1;
    if ( $qryLen != 0 && $sbjLen != 0 ) {
      my $deletions = ( $result->getPctDelete() / 100 ) * $qryLen;
      $result->setPctDelete(
                           sprintf( "%4.2f", ( $deletions * 100 ) / $sbjLen ) );
    }

    #
    # Calculate Kimura divergence using the CpG modification described
    # in SearchResult.pm
    #
    my ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
        $result->calcKimuraDivergence( divCpGMod => 1 );
    $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
  }

  if ( $resultCode == 141 ) {

    # It seems that in cross_match 1.08+ Phil changed the result code
    # behaviour.  Now it returns 141 in cases where no results are
    # found in the search, 0 if results are found, and some other
    # value if something goes wrong.
    $resultCode = 0;
  }

  return ( $resultCode, $searchResultsCollection );
}

##---------------------------------------------------------------------##

=head1 CLASS METHODS

=cut

##---------------------------------------------------------------------##

##---------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = parseOutput(
                                     searchOutput => $filename|$FH,
                                     [matrixName => $matrixName],
                                     [excludeAlignments => 1],
                                               );

  Parse the result of a search and return a SearchResultCollection.

=cut

##---------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $CMFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle|IO::File/ ) {
    open $CMFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $CMFILE = $nameValueParams{'searchOutput'};
  }

  my $callbackFunc = undef;
  if ( defined $nameValueParams{'callback'}
       && ref( $nameValueParams{'callback'} ) == /CODE/ )
  {
    $callbackFunc = $nameValueParams{'callback'};
  }

  #
  # Three versions of RM/Crossmatch output
  #
  # Out Files, and old Align format
  my @outfileFwdStrandKeys = qw( score pctDiverge pctDelete pctInsert queryName
      queryStart queryEnd queryRemaining orientation subjName subjType
      subjStart subjEnd subjRemaining id overlap );
  my @outfileRevStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining orientation
      subjName subjType subjRemaining subjEnd
      subjStart id overlap );

  # New align format, Cat File format ( id field not used )
  my @alignfileFwdStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining subjName
      subjStart subjEnd subjRemaining lineageId id );
  my @alignfileRevStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining orientation
      subjName subjRemaining subjEnd
      subjStart lineageId id  );

  # Standard crossmatch format
  my @crossmatchFwdStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining subjName
      subjStart subjEnd subjRemaining id overlap );
  my @crossmatchRevStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining orientation
      subjName subjRemaining subjEnd
      subjStart id overlap );

  my $seqName;
  my $querySeq;
  my $transI;
  my $transV;
  my $kimura;
  my $subjSeq;
  my $result;
  my $alignPos          = 0;
  my $sourceFileIndex   = 0;    # The index of the alignment
                                #  within the source file.
  my $queryBackgroundGC = 0;
  my $queryComplemented = 0;
  my $matrix            = "";
  my @hdrLineArray      = ();
  my $alignmts;

  my $resultColl = SearchResultCollection->new();

  while ( <$CMFILE> ) {

    # Crossmatch 1.080812 added some reporting features at the end
    # of the output that look to the parser like additional annotations.
    # This statements ends our search before we hit the end of the file
    # so we avoid having to diambiguate these lines.
    last if ( /^Score histogram:/ );

    #
    # RepeatMasker adds this to the crossmatch alignment output
    #  TODO: REMOVE: I don't believe this is the case anymore
    $queryBackgroundGC = $1 if ( /^Assumed background.*is (\d+) %/ );

    #
    # If we have complete output we can also grab the matrix
    # used.
    #
    # It used to be that I output "Score matrix" before all
    # alignments that shared the same matrix.  There is probably
    # some code out there that referenced this.  Just a note
    # in case we ever ressurect that code.  If so put this
    # back and don't clear the $matrix parameter after each
    # hit below.
    #if ( /Score matrix\s+(\S+)/ ) {
    #
    # Now I store the matrix field per alignment.
    if ( /^Matrix\s*=\s*(\S+)/ ) {
      my @path = split( /[\\\/]/, $1 );
      $matrix = $path[ $#path ];
      $matrix =~ s/[\n\r]//g;
    }

    # Kimura Divergence
    #  - Non crossmatch field
    if ( /^Kimura.*=\s*(\S+)/ ) {
      $kimura = $1;
    }

    #
    # Look for transition/transversion ratio
    #
    if ( /Transitions.*\((\d+)\s*\/\s*(\d+)\)/ ) {
      $transI = $1;
      $transV = $2;
    }

    # Look for a header line (e.g.):
    if ( /^\s*\d+\s+\d+(\.\d+)?/ ) {
      @hdrLineArray = split;    # Create an array out of this line
           # Check to see if we should filter out low scoring hits
           #unless (  ( exists $this->{scoreHighThresh} &&
           #            $hdrLineArray[0] >= $this->{scoreHighThresh}) ||
           #          ( exists $this->{scoreLowThresh} &&
           #            $hdrLineArray[0] <= $this->{scoreLowThresh}) ||
           #          ( exists $this->{subjPattern} &&
           #           (($hdrLineArray[8] =~ /$this->{subjPattern}/o)<1 &&
           #           ($hdrLineArray[9] =~ /$this->{subjPattern}/o)<1)) ||
           #          ( exists $this->{queryPattern} &&
           #            ($hdrLineArray[4] =~ /$this->{queryPattern}/o)<1 ) ||
           #          ( exists $this->{gcLowThresh} &&
           #           $queryBackgroundGC < $this->{gcLowThresh}) ||
           #          ( exists $this->{gcHighThresh} &&
           #           $queryBackgroundGC > $this->{gcHighThresh})) {

      if ( $#hdrLineArray > 10 ) {
        my @dataArray      = @hdrLineArray;
        my @fieldKeys      = ();
        my %nameValuePairs = ();

        # Is this a crossmatch line, an outfile line, or
        # and *.align line?
        #   crossmatch output:
        #     11/12 fields
        #   *.out files:
        #     Forward strand results have a "+" orientation field
        #     which is not the case in crossmatch files.
        #     Out files also have a non-numeric type field following
        #     the subject name.
        #         15/16 fields ( overlap field is 16th )
        #   *.align files:
        #     Old files ended in simple ids ( same id as in the *.out file )
        #        15 fields always
        #     New files end in batchID followed by *.out id:
        #         .... m_b2s452i0 1  ( 14/15 fields - 15 if rev strand )
        #   *.cat files:
        #      Old..
        #      New same as *.align except for no id column
        #              12/13 fields
        #
        if ( $hdrLineArray[ 8 ] eq "+"
             || !( $hdrLineArray[ 10 ] =~ /^[\(\)\d]+$/ ) )
        {

          # Definitely an *.out or old *.align line
          if ( $hdrLineArray[ 8 ] eq "+" ) {
            $nameValuePairs{'orientation'} = "";
            @fieldKeys = @outfileFwdStrandKeys;
          }
          else {
            $nameValuePairs{'orientation'} = "C";
            @fieldKeys = @outfileRevStrandKeys;
          }
        }
        elsif (
                (
                  (
                       $hdrLineArray[ 8 ] eq "C"
                    && defined $hdrLineArray[ 13 ]
                    && $hdrLineArray[ 13 ] =~ /\[?[mc]_b.*/
                  )
                  || ( defined $hdrLineArray[ 12 ]
                       && $hdrLineArray[ 12 ] =~ /\[?[mc]_b.*/ )
                )
            )
        {

          # A new *.align file or cat file ( id field unused )
          if ( $hdrLineArray[ 8 ] eq "C" ) {
            $nameValuePairs{'orientation'} = "C";
            @fieldKeys = @alignfileRevStrandKeys;
          }
          else {
            $nameValuePairs{'orientation'} = "";
            @fieldKeys = @alignfileFwdStrandKeys;
          }
        }
        else {

          # Probably a crossmatch file
          if ( $hdrLineArray[ 8 ] eq "C" ) {
            $nameValuePairs{'orientation'} = "C";
            @fieldKeys = @crossmatchRevStrandKeys;
          }
          else {
            $nameValuePairs{'orientation'} = "";
            @fieldKeys = @crossmatchFwdStrandKeys;
          }
        }

        map {
          my $datum = shift( @dataArray );
          my $field = $_;
          if ( $field ne "orientation" ) {
            $datum =~ s/[\(\)]//g
                if (    $field eq "subjRemaining"
                     || $field eq "queryRemaining" );
            if ( $#dataArray == 0 && $datum eq "*" ) {
              $field = 'id';
            }
            $nameValuePairs{$field} = $datum;
          }
        } @fieldKeys;
        if ( defined $result && defined $callbackFunc ) {
          $callbackFunc->( $result );
        }

        $result = SearchResult->new( %nameValuePairs );
        if ( !defined $callbackFunc ) {
          $resultColl->add( $result );
        }
        $querySeq = "";
        $subjSeq  = "";
      }

      @hdrLineArray = ();

    }

    #
    # Concatenate the alignment sequences
    #
    if ( /^(C?)\s+(\S+)\s+\d+\s+(\S+)\s+\d+\s*$/
         && !exists $nameValueParams{'excludeAlignments'} )
    {
      if ( $alignPos == 0 ) {

        # This is the query sequence
        $querySeq .= $3;
        $queryComplemented = 1 if ( $1 eq "C" );
      }
      else {

        # This is the subj sequence
        $subjSeq .= $3;
      }
      $alignPos ^= 1;
    }

    #
    # Look for a signal for the end of an alignment
    #
    if ( /Gap_init rate/ && defined $result ) {

      # Store this alignment in our data structure.
      if ( $queryComplemented ) {

        # Crossmatch does not complement the subject
        # sequence in it's alignments ( only the query
        # sequence).  However our SearchResult object
        # requires that the query sequence be in
        # the forward direction.  So...do a simple
        # reversal of the sequence so that we can
        # store it in the object.
        $querySeq = reverse $querySeq;
        $querySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
        $subjSeq = reverse $subjSeq;
        $subjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;     # complement
      }
      $result->setQueryString( $querySeq );
      $result->setSubjString( $subjSeq );

      if ( $matrix ne "" ) {
        $result->setMatrixName( $matrix );
      }
      elsif ( $nameValueParams{'matrixName'} ) {
        $result->setMatrixName( $nameValueParams{'matrixName'} );
      }

      if ( defined $kimura ) {
        $result->setPctKimuraDiverge( $kimura );
      }

      if ( defined $callbackFunc ) {
        $callbackFunc->( $result );
      }

      $result            = undef;
      $querySeq          = "";
      $subjSeq           = "";
      $matrix            = "";
      $transV            = 0;
      $transI            = 0;
      $kimura            = undef;
      $alignPos          = 0;
      $queryComplemented = 0;
      @hdrLineArray      = ();
    }
  }
  close $CMFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Use:  my _systemint( $cmd );
##
##     Interruptible system call routine.
##
##  Returns
##
##  Globals Used: None
##-------------------------------------------------------------------------##?
sub _systemint {
  my ( $cmd ) = @_;
  my ( $pid );
  my ( $flag ) = 0;

  local $SIG{INT}  = sub { &_handler( @_ ) if ( $flag ) };    #^C
  local $SIG{QUIT} = sub { &_handler( @_ ) if ( $flag ) };    #^\
  local $SIG{TERM} =
      sub { &_handler( @_ ) if ( $flag ) };    #kill command or system crash
  local $SIG{HUP} = sub { &_handler( @_ ) if ( $flag ) };

FORK:
  {
    if ( $pid = fork ) {
      $flag = 1;
      waitpid( $pid, 0 );                      #Waits for child to finish...
      my ( $status ) = $?;
      if ( WIFSTOPPED( $status ) ) {
        my ( $signal ) = WSTOPSIG( $status );
        print "\nforksys:  Program terminated by a signal $signal.\n";
        print "The executing command was:  $cmd\n";
        return 1;
      }
      if ( WIFEXITED( $status ) ) {
        my ( $temp ) = WEXITSTATUS( $status );
        return $temp;
      }
      if ( WIFSIGNALED( $status ) ) {
        my ( $signal ) = WTERMSIG( $status );
        return $signal;
      }
      if ( WIFSIGNALED( $status ) ) {
        my ( $signal ) = WTERMSIG( $status );
        print "\nforksys:  Program terminated by a signal $signal.\n";
        print "The executing command was:  $cmd\n";
        return 1;
      }
    }
    elsif ( defined $pid ) {
      exec( "$cmd" ) or die "Exec $cmd failed\n";
    }
    elsif ( $! =~ /No more process/o ) {
      print "$!\n";
      sleep 5;
      redo FORK;
    }
    else {
      die "Can't fork...\n";
    }
  }
}
##-------------------------------------------------------------------------##
## Use: my _handler( $sig );
##
##  Interrupt handler used by systemint() ###
##
##  Returns
##
##  Globals Used: None
##-------------------------------------------------------------------------##?
sub _handler {
  my ( $sig ) = @_;

  print $CLASS. "_handler(): Aborting with a SIG$sig\n";
  exit( -1 );
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
