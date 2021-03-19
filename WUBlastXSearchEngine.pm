#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) WUBlastXSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      the WU-BlastX search engine.
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
#      'WUBlastXSearchEngine' );
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

WUBlastXSearchEngine

=head1 SYNOPSIS

use WUBlastXSearchEngine

Usage: 

  use SearchEngineI;
  use WUBlastXSearchEngine;
  use SearchResultCollection;

  my $wuEngine = WUBlastXSearchEngine->new( 
                    pathToEngine=>"/usr/local/wublast/blastx" );

  $wuEngine->setMatrix( "/users/bob/simple.matrix" );
  $wuEngine->setQuery( "/users/bob/query.fasta" );
  $wuEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $wuEngine->search();

=head1 DESCRIPTION

  A concrete implementation of the abstract class / interface SearchEngineI
  which use the WUBlastX sequence search engine.

=head1 INSTANCE METHODS

=cut 

package WUBlastXSearchEngine;
use strict;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use FileHandle;
use File::Basename;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Constants
#
use constant Query   => 1;
use constant Subject => 2;

#
# Version
#
my $VERSION = 0.1;
my $CLASS   = "WUBlastXSearchEngine";
my $DEBUG   = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/blastp\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setOverrideParameters()

  Use: my $value    = getOverrideParameters( );
  Use: my $oldValue = setOverrideParameters( $value );

  Get/Set the the override paramters.  These are used instead
  of all the SearchEngineI default parameters if set.

=cut

##-------------------------------------------------------------------------##
sub getOverrideParameters {
  my $this = shift;

  return $this->{'overrideParameters'};
}

sub setOverrideParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'overrideParameters'};
  $this->{'overrideParameters'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 MaskLevelSequence()

  Use: my $value    = getMaskLevelSequence( );
  Use: my $oldValue = setMaskLevelSequence( $value );

  Get/Set the MaskLevelSequence paramter.  This is the
  sequence ( SearchResult::Query or SearchResult::Subject )
  which will be considered when applying the mask level.  The
  default is SearchResult::Query.

=cut

##-------------------------------------------------------------------------##
sub getMaskLevelSequence {
  my $this = shift;

  return $this->{'maskLevelSequence'};
}

sub setMaskLevelSequence {
  my $this  = shift;
  my $value = shift;

  croak $CLASS
      . "::setMaskLevelSequence: Invalid value ( $value ). "
      . "Should be either SearchResult::Query or "
      . "SearchResult::Subject\n"
      if (    $value != SearchResult::Query
           && $value != SearchResult::Subject );

  my $oldValue = $this->{'maskLevelSequence'};
  $this->{'maskLevelSequence'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 PValueCutoff()

  Use: my $value    = getPValueCutoff( );
  Use: my $oldValue = setPValueCutoff( $value );

  Get/Set the PValueCutoff.

=cut

##-------------------------------------------------------------------------##
sub getPValueCutoff {
  my $this = shift;

  return $this->{'PValueCutoff'};
}

sub setPValueCutoff {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'PValueCutoff'};
  $this->{'PValueCutoff'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 FilterWords()

  Use: my $value    = getFilterWords( );
  Use: my $oldValue = setFilterWords( $value );

  Turn on / off the Dust/Seg/Xnu screening of words.

=cut

##-------------------------------------------------------------------------##
sub getFilterWords {
  my $this = shift;

  return $this->{'filterWords'};
}

sub setFilterWords {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'filterWords'};
  $this->{'filterWords'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 TempDir()

  Use: my $value    = getTempDir( );
  Use: my $oldValue = setTempDir( $value );

  Set the directory to use as a temp directory for a search.  
  The default is to use the directory which contains the query sequence.

=cut

##-------------------------------------------------------------------------##
sub getTempDir {
  my $this = shift;

  return $this->{'tempDir'};
}

sub setTempDir {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'tempDir'};
  $this->{'tempDir'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAdditionalParameters()

  Use: my $value    = getAdditionalParameters( );
  Use: my $oldValue = setAdditionalParameters( $value );

  Get/Set the additional paramters.  These are used in addition
  to the existing parameter set. 

=cut

##-------------------------------------------------------------------------##
sub getAdditionalParameters {
  my $this = shift;

  return $this->{'additionalParameters'};
}

sub setAdditionalParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'additionalParameters'};
  $this->{'additionalParameters'} = $value;

  return $oldValue;
}

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
  if ( $result =~ /^(BLASTX) (\S+ \[.*\])/m ) {
    $this->{'version'}    = $2;
    $this->{'engineName'} = $1;
  }
  else {
    croak $CLASS
        . "::setPathToEngine( $value ): Cannot determine "
        . "engine variant and version!\n";
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## Instance Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or                                                                            
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

  Run the search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub search {
  my $this           = shift;
  my %nameValuePairs = @_;

  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::search: Instance variable $name doesn't exist." );
      }
      $this->$method( $value );
    }
  }

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  # Generate parameter line
  my $parameters    = "";
  my $spanParameter = "";
  my $value;
  if ( ( $value = $this->getSubject() ) ) {

    # Make sure we have the compressed form of the database handy
    if ( -f "$value.ahd" || -f "$value.xns" || -f "$value.xps" ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS
          . "::search: Error...compressed subject "
          . "database ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  my $outputDirName;
  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
    if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
      $outputDirName = $this->getTempDir();
    }
    else {
      $outputDirName = dirname( $value );
    }
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }

  # Constant parameters
  $parameters .= " -warnings -T=1000000 hspmax=2000000000 B=100000000";
  if ( $this->{'filterWords'} >= 1 ) {
    $parameters .= " wordmask=xnu";
  }

  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $parameters .= " Q=" . abs( $value );
  }

  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " R=" . abs( $value );
  }

  if ( defined( $value = $this->getMinMatch() )
       && $value =~ /\d+/ )
  {
    $parameters .= " W=$value";
  }

  #     S:  Overall score threshold...not sure how this relates to the others
  #    S2:  Score threshold for ungapped HSPs
  # gapS2:  Score threshold for gapped HSPs
  #     X:  Dropoff score for ungapped HSPs
  #  gapX:  Dropoff score for gapped HSPs
  #  if ( defined( $value = $this->getMinScore() )
  #       && $value =~ /\d+/ )
  #  {
  #    $parameters .= " S=$value";
  #    $parameters .= " gapS2=$value";
  #    $parameters .= " S2=" . int( $value / 2 );
  #    $parameters .= " X=$value";
  #    $parameters .= " gapX=" . ( $value * 2 );
  #  }

  if ( defined( $value = $this->getBandwidth() )
       && $value =~ /\d+/ )
  {
    $parameters .= " gapW=" . ( ( $value * 3 ) + 1 );
  }

  if ( defined( $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # WUBLAST requires that the matrix filename parameter
      # be relative to a directory path specified in
      # environment variables.
      my @path   = split( /[\\\/]/, $value );
      my $matrix = pop @path;

      # WUBLAST expects the path to be above "aa" or "nt" sub
      # directories
      if ( $path[ $#path ] eq "aa" || $path[ $#path ] eq "nt" ) {
        pop @path;
      }

      # Set the environment
      $ENV{BLASTMAT}   = join( "/", @path );
      $ENV{WUBLASTMAT} = $ENV{BLASTMAT};

      $parameters .= " -matrix=$matrix";

    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  # Invoke engine and handle errors
  my $POUTPUT = new FileHandle;
  my $errFile;
  do {
    $errFile = $outputDirName . "/wuResults-" . time() . ".err";
  } while ( -f $errFile );
  my $pid;
  my $runParameters;
  if ( defined $this->{'overrideParameters'}
       && $this->{'overrideParameters'} ne "" )
  {
    $runParameters = $this->{'overrideParameters'};
  }
  else {
    $runParameters = $parameters . " ";
  }

  if ( defined $this->{'additionalParameters'}
       && $this->{'additionalParameters'} ne "" )
  {
    $runParameters .= " " . $this->{'additionalParameters'};
  }

  print $CLASS
      . "::search() Invoking search engine as: $engine "
      . "$runParameters 2>$errFile |\n"
      if ( $DEBUG );

  $pid = open( $POUTPUT, "$engine $runParameters 2>$errFile |" );

  # Create SearchResultCollection object from
  # the engine results.
  my $searchResultsCollection;
  if ( $DEBUG ) {
    ## Create a debug file
    my $outFile;
    do {
      $outFile = $outputDirName . "/wuResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $searchResultsCollection = parseOutput(
                                       searchOutput => $outFile,
                                       pValueCutoff => $this->getPValueCutoff(),
                                       scoreMode    => $this->getScoreMode(),
                                       scoreCutoff  => $this->getMinScore()
    );

  }
  else {

    # Just pipe to the parser
    $searchResultsCollection = parseOutput(
                                       searchOutput => $POUTPUT,
                                       pValueCutoff => $this->getPValueCutoff(),
                                       scoreMode    => $this->getScoreMode(),
                                       scoreCutoff  => $this->getMinScore()
    );
  }
  close $POUTPUT;

  my $resultCode = ( $? >> 8 );
  print "WUBlast returned a the following result code >$resultCode<\n"
      if ( $DEBUG );

  # This result code can occur if wublast encounters a sequence
  # (at the end of a query collection) which is short enough
  # that it would never score higher than the threshold. An
  # exit code of 16 | 17 occurs with a printed error of
  # "FATAL:   There are no valid contexts in the requested search."
  # Look for this and ignore.  NOTE: In WUBLast 2.0 ( non free version )
  # you can use the -novalidctxok option to warn but not return a bad exit
  # code.
  if ( $resultCode != 0 ) {
    open ERR, "<$errFile";
    my $errOk = 1;
    while ( <ERR> ) {
      if ( /^FATAL:/ ) {
        if ( !/context|cpus|P option|uses P|shorter/i ) {
          $errOk = 0;
          print STDERR $CLASS . "::search: $_\n";
        }
      }
    }
    close ERR;
    $resultCode = 0 if ( $errOk == 1 );
  }
  unlink $errFile unless ( $DEBUG );

  print $CLASS
      . "::search: "
      . $searchResultsCollection->size()
      . " hits " . "\n"
      if ( $DEBUG );

  #
  # Postprocess the results
  #
  if ( defined $searchResultsCollection
       && $searchResultsCollection->size() > 0 )
  {

    #
    # Mask level filtering
    #
    my $maskLevel;
    my $strand;
    if ( defined( $maskLevel = $this->getMaskLevel() ) ) {
      if ( defined( $strand = $this->getMaskLevelSequence() ) ) {
        $searchResultsCollection->maskLevelFilter( value  => $maskLevel,
                                                   strand => $strand );
      }
      else {
        $searchResultsCollection->maskLevelFilter( value => $maskLevel );
      }
      print $CLASS
          . "::search: "
          . $searchResultsCollection->size()
          . " hits "
          . "after masklevel filtering\n"
          if ( $DEBUG );
    }

    # Secondary sort by query position
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
            || $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart();
      }
    );

  }

  return ( $resultCode, $searchResultsCollection );
}

##-------------------------------------------------------------------------##

=head1 Class Methods

=cut

##-------------------------------------------------------------------------##

##---------------------------------------------------------------------##

=head2 parseOutput()

  Use: my \@results = &parseWublastXOutput( searchOutput => "",
                                            pValueCutoff => #,
                                            scoreMode => "",
                                            scoreCutoff => #
                                          );

        searchOutput   : WublastX results file or FH
        pValueCutoff   : Lowest pValue to include in results

    Returns
       A SearchResultCollection with all the results.

=cut

##---------------------------------------------------------------------##
sub parseOutput {
  my %parameters = @_;

  croak $CLASS. "::parseWublastXOutput() missing searchOutput parameter!\n"
      if ( !exists $parameters{'searchOutput'} );

  my $WUFILE;
  if ( ref( $parameters{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseWublastXOutput() Opening file "
        . $parameters{'searchOutput'} . "\n"
        if ( $DEBUG );
    open $WUFILE, $parameters{'searchOutput'}
        or die $CLASS
        . "::parseWublastXOutput: Unable to open "
        . "results file: $parameters{'searchOutput'} : $!";
  }
  else {
    $WUFILE = $parameters{'searchOutput'};
  }

  my $cutoffP = -1;
  if ( defined $parameters{'pValueCutoff'} ) {
    $cutoffP = $parameters{'pValueCutoff'};
  }

  my @results = ();
  my $queryID;
  my $qryStart;
  my $qrySeq;
  my $qryEnd;
  my $pValue;
  my $score;
  my $orientation;
  my $sbjStart;
  my $sbjEnd;
  my $sbjSeq;
  my $subjectID;
  my $resultColl = SearchResultCollection->new();

  my ( $matrixRef, $matrixAlphabet, $matrixFreqRef ) = _readMatrix();
  my $lambda = _calculateLambda( $matrixRef, $matrixFreqRef );

  #
  while ( <$WUFILE> ) {
    chomp;

    if ( defined $qryStart
         && ( /^Query=/ || /^>/ || /^ Score =/ || /Strand HSPs/ ) )
    {
      if ( $cutoffP < 0 || $pValue <= $cutoffP ) {

        # Complexity adjust the score
        if ( defined $parameters{'scoreMode'}
             && $parameters{'scoreMode'} ==
             SearchEngineI::complexityAdjustedScoreMode )
        {
          $score =
              _complexityAdjust( $score, $qrySeq, $sbjSeq, $lambda,
                                 $matrixAlphabet, $matrixRef, $matrixFreqRef );
        }

        if ( !defined $parameters{'scoreCutoff'}
             || $score >= $parameters{'scoreCutoff'} )
        {

          # Save the result
          my $result = SearchResult->new(
            queryName  => $queryID,
            queryStart => $qryStart,
            queryEnd   => $qryEnd,

            #queryRemaining => ( $qryLength - $qryEnd ),
            #queryString    => $qrySeq,
            #subjString     => $sbjSeq,
            orientation => $orientation,
            subjName    => $subjectID,
            subjStart   => $sbjStart,
            subjEnd     => $sbjEnd,

            #subjRemaining  => ( $sbjLength - $sbjEnd ),
            #queryString    => $qrySeq,
            #pctDiverge     => $percDiv,
            #pctInsert      => $percIns,
            #pctDelete      => $percDel,
            pValue => $pValue,
            score  => $score
          );
          $resultColl->add( $result );

        }
      }
      undef $qryStart;
      undef $qryEnd;
      $sbjStart = "";
      $sbjEnd   = "";
      $qrySeq   = "";
      $sbjSeq   = "";
    }

    if ( /^Query=\s+(\S+)/ ) {
      $queryID = $1;
    }
    elsif ( /^>(\S+)/ ) {
      $subjectID = $1;
    }
    elsif ( /^\s+(\w+) Strand HSPs:\s*$/ ) {

      # TODO: Change this to +/-
      if ( $1 eq "Minus" ) {
        $orientation = "C";
      }
      else {
        $orientation = "";
      }
    }
    elsif ( /^ Score = (\d+) .+ (\S+)\s*$/ ) {
      $score  = $1;
      $pValue = $2;
    }
    elsif ( /^Query:\s+(\d+)\s+(\S+)\s+(\d+)\s*$/ ) {

      $qrySeq .= $2;
      if ( $qryStart < 1 ) {
        $qryStart = _min( $1, $3 );
        $qryEnd   = _max( $1, $3 );
      }
      else {
        $qryStart = _min( _min( $1, $3 ), $qryStart );
        $qryEnd   = _max( _max( $1, $3 ), $qryEnd );
      }

      #$orientation = "C" if ( $1 > $3 );
    }
    elsif ( /Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/ ) {

      $sbjSeq .= $2;
      my $leftNum  = $1;
      my $rightNum = $3;
      if ( $sbjStart < 1 ) {
        $sbjStart = _min( $leftNum, $rightNum );
        $sbjEnd   = _max( $leftNum, $rightNum );
      }
      else {
        $sbjStart = _min( _min( $leftNum, $rightNum ), $sbjStart );
        $sbjEnd   = _max( _max( $leftNum, $rightNum ), $sbjEnd );
      }
    }
  }
  if ( defined $qryStart && ( $cutoffP < 0 || $pValue <= $cutoffP ) ) {

    # Complexity adjust the score
    if ( defined $parameters{'scoreMode'}
         && $parameters{'scoreMode'} ==
         SearchEngineI::complexityAdjustedScoreMode )
    {
      $score =
          _complexityAdjust( $score, $qrySeq, $sbjSeq, $lambda, $matrixAlphabet,
                             $matrixRef, $matrixFreqRef );
    }

    if ( !defined $parameters{'scoreCutoff'}
         || $score >= $parameters{'scoreCutoff'} )
    {

      # Save the result
      my $result = SearchResult->new(
        queryName  => $queryID,
        queryStart => $qryStart,
        queryEnd   => $qryEnd,

        #queryRemaining => ( $qryLength - $qryEnd ),
        #queryString    => $qrySeq,
        #subjString     => $sbjSeq,
        orientation => $orientation,
        subjName    => $subjectID,

        #subjStart      => $sbjStart,
        #subjEnd        => $sbjEnd,
        #subjRemaining  => ( $sbjLength - $sbjEnd ),
        #queryString    => $qrySeq,
        #pctDiverge     => $percDiv,
        #pctInsert      => $percIns,
        #pctDelete      => $percDel,
        pValue => $pValue,
        score  => $score
      );
      $resultColl->add( $result );
    }
  }

  close $WUFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _complexityAdjust( $swnScore, $qrySeq, $sbjSeq, $matLambda,
##                            $matAlphabet, $matScoresRef, $matFreqsRef );
##
##      $swnScore    : Smith-Waterman Raw Alignment Score
##      $qrySeq      : Query string from the alignment
##      $sbjSeq      : Subject string from the alignment
##      $matLambda   : Matrix Lambda parameter
##      $matAlphabet : Matrix alphabet string (in column order)
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $adj_score : The complexity adjusted score according to
##                   to Phil Green's swat/cross_match program.
##-------------------------------------------------------------------------##
sub _complexityAdjust {
  my ( $swnScore, $qrySeq, $sbjSeq, $matLambda, $matAlphabet, $matScoresRef,
       $matFreqsRef )
      = @_;
  my $mCol      = 0;
  my $t_factor  = 0;
  my $t_sum     = 0;
  my $t_counts  = 0;
  my $n_letters = 0;
  my $baseIndex = 0;
  my $qBase     = "";
  my $sBase     = "";
  my $pos_score = 0;
  my @matCounts = ();
  my $adj_score = 0;

  #
  # Recalculates the raw score ( $pos_score ) based on the matrix we
  # have been handed.
  #
  # Creates a vector of base counts from the query sequence.  It
  # ignores base instances when they are part of a deletion or
  # are insertion characters "-".
  #
  for ( $baseIndex = 0 ; $baseIndex < length( $qrySeq ) ; $baseIndex++ ) {
    $qBase = substr( $qrySeq, $baseIndex, 1 );
    $sBase = substr( $sbjSeq, $baseIndex, 1 );
    if ( $qBase ne "-" && $sBase ne "-" ) {
      $pos_score +=
          $$matScoresRef[ index( $matAlphabet, $qBase ) ]
          [ index( $matAlphabet, $sBase ) ];
      $matCounts[ index( $matAlphabet, $qBase ) ]++
          if ( ( index $matAlphabet, $qBase ) >= 0 );
    }
  }

  #print "Recalculated raw score = $pos_score\n";

  #
  #
  #
  for ( $mCol = 0 ; $mCol < length( $matAlphabet ) ; $mCol++ ) {
    if ( defined $matCounts[ $mCol ] && $matCounts[ $mCol ] > 0 ) {
      if ( $$matFreqsRef[ $mCol ] > 0 && log( $$matFreqsRef[ $mCol ] ) != 0 ) {
        $t_factor += $matCounts[ $mCol ] * log( $matCounts[ $mCol ] );
        $t_sum    += $matCounts[ $mCol ] * log( $$matFreqsRef[ $mCol ] );
        $t_counts += $matCounts[ $mCol ];
        $n_letters++;
      }
    }
  }

  #
  #
  #
  $t_factor -= $t_counts * log( $t_counts );
  $t_sum    -= $t_factor;

  #
  # Looks like Phil changed his mind here
  #
  #my $complexity_factor = 0.25;
  #if ( $n_letters > ( 1 / $complexity_factor )) {
  #  $t_factor /= $t_counts * log(1 / $n_letters);
  #}else {
  #  $t_factor /= $t_counts * log($complexity_factor);
  #}
  #
  #$old_adj_score = $swnScore + ($t_factor * $pos_score) - $pos_score + .5;

  #
  #
  #
  $adj_score = sprintf( "%0.0d", $swnScore + $t_sum / $matLambda + .999 );

  $adj_score = 0 if ( !( $adj_score =~ /\d+/ ) || $adj_score < 0 );
  return ( $adj_score );
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
  my ( $matScoresRef, $matFreqsRef ) = @_;
  my $lambda_upper = 0;
  my $lambda_lower = 0;
  my $lambda       = 0.5;
  my $sum          = 0;

  do {
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum < 1.0 ) {
      $lambda_lower = $lambda;
      $lambda *= 2.0;
    }
  } while ( $sum < 1.0 );

  $lambda_upper = $lambda;

  while ( $lambda_upper - $lambda_lower > .00001 ) {
    $lambda = ( $lambda_lower + $lambda_upper ) / 2.0;
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum >= 1.0 ) {
      $lambda_upper = $lambda;
    }
    else {
      $lambda_lower = $lambda;
    }
  }
  return ( $lambda );
}

##-------------------------------------------------------------------------##
## Use: my _getSum( $lambda, $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet frequency vector
##
## Returns
##      $sum : Good question....???  This is used by the
##             lambda estimation procedure.  _getSum is
##             derived from Phil Green's swat/cross_match
##             programs.
##-------------------------------------------------------------------------##
sub _getSum {
  my ( $lambda, $matScoresRef, $matFreqsRef ) = @_;
  my $check = 0;
  my $total = 0;
  my $i     = 0;
  my $j     = 0;

  for ( $i = 0 ; $i <= $#$matFreqsRef ; $i++ ) {
    for ( $j = 0 ; $j <= $#$matFreqsRef ; $j++ ) {
      if ( $matFreqsRef->[ $i ] && $matFreqsRef->[ $j ] ) {
        $total +=
            $matFreqsRef->[ $i ] * $matFreqsRef->[ $j ] *
            exp( $lambda * $matScoresRef->[ $i ]->[ $j ] );
        $check += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ];
      }
    }
  }
  die "error in _getSum!! check=$check\n"
      if (    $check > 1.001
           || $check < .999 );

  return ( $total );
}

## TODO: Use SequenceSimilarityMatrix.pm here
##-------------------------------------------------------------------------##
## Use: my ( \@matrix, $alphabet, \@freqArray ) =
##                                        _readMatrix( $matrixFileName );
##
##      $matFileName : Name of the file containing the WUBlast Matrix
##
## Returns
##      This is a very specific matrix reader for WUBlast matrices
##      which have been created for RepeatMasker.  This routine
##      assumes that the standard crossmatch matrix frequency line
##      is included in the comment section. Ie.:
##
##        # FREQS A 0.325 C 0.175 G 0.175 T 0.325
##
##      The routine returns the matrix stored as a 2-d array.  It also
##      returns the alphabet as a compressed string and the frequencies
##      as an array ( sizeof alphabet ).
##
##-------------------------------------------------------------------------##
sub _readMatrix {
  my @matrix    = ();
  my $alphabet  = "";
  my @freqArray = ();

  $alphabet = "ARNDCQEGHILKMFPSTWYVBZX";
  @matrix = (
    [ qw(4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0) ],
    [
      qw(-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1)
    ],
    [
      qw(-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1)
    ],
    [
      qw(-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1)
    ],
    [
      qw( 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2)
    ],
    [
      qw(-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1)
    ],
    [
      qw(-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1)
    ],
    [
      qw( 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1)
    ],
    [
      qw(-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1)
    ],
    [
      qw(-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1)
    ],
    [
      qw(-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1)
    ],
    [
      qw(-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1)
    ],
    [
      qw(-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1)
    ],
    [
      qw(-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1)
    ],
    [
      qw(-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2)
    ],
    [
      qw( 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0)
    ],
    [
      qw( 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0)
    ],
    [
      qw(-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2)
    ],
    [
      qw(-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1)
    ],
    [
      qw( 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1)
    ],
    [
      qw(-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1)
    ],
    [
      qw(-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1)
    ],
    [
      qw( 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1)
    ]
  );

  @freqArray =
      qw( 0.0777 0.0627 0.0336 0.0542 0.0078 0.0315 0.0859 0.0730 0.0192 0.0666 0.0891 0.0776 0.0241 0.0361 0.0435 0.0466 0.0487 0.0102 0.0300 0.0817 0.0002 );

  return ( \@matrix, $alphabet, \@freqArray );
}

##-------------------------------------------------------------------------##
## Use: my _min( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The minimum of the two numbers
##
##-------------------------------------------------------------------------##
sub _min {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num1 );
  }
  else {
    return ( $num2 );
  }
}

##-------------------------------------------------------------------------##
## Use: my _max( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The maximum of the two numbers
##
##-------------------------------------------------------------------------##
sub _max {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num2 );
  }
  else {
    return ( $num1 );
  }
}

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
