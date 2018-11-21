#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) HMMERSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      HMMER3.1 nucleotide search engines nhmmer/nhmmscan
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2012 Developed by
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
#      'HMMERSearchEngine' );
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

HMMERSearchEngine

=head1 SYNOPSIS

use HMMERSearchEngine

Usage: 

  use SearchEngineI;
  use HMMERSearchEngine;
  use SearchResultCollection;

  my $hmmerEngine = HMMERSearchEngine->new( 
                    pathToEngine=>"/usr/local/hmmer3.1/nhmmer" );

or 
  my $hmmerEngine = HMMERSearchEngine->new( 
                    pathToEngine=>"/usr/local/hmmer3.1/nhmmscan" );

  $hmmerEngine->setQuery( "/users/bob/query.fasta" );
  $hmmerEngine->setSubject( "/users/bob/subject.hmm" );
  my $searchResults = $hmmerEngine->search();

=head1 DESCRIPTION

  A concrete implementation of the abstract class / interface SearchEngineI
  which use the HMMER sequence search engine.

=head1 INSTANCE METHODS

=cut 

package HMMERSearchEngine;
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
# Version
#
my $VERSION = 1.0;
my $CLASS   = "HMMERSearchEngine";

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

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::set: Instance variable $name doesn't exist." . "" );
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

  my $result = `$value -h 2>&1`;

  # HMMER 3.0 (March 2010); http://hmmer.org/
  # HMMER 3.1dev (March 2010); http://hmmer.org/
  # HMMER 3.1dev_0.30_size (March 2010); http://hmmer.org/
  # HMMER hmmer3.1-snap20120830 (August 2012); http://hmmer.org/
  while ( $result =~ /([^\n\r]*)[\n\r]/ig ) {
    my $line = $1;
    if (    $line =~ /^#\s+HMMER\s+(\d+\S+\s+\(.*\)).*/
         || $line =~ /^#\s+HMMER\s+(hmmer\S+\s+\(.*\)).*/ )
    {
      $this->{'version'}    = $1;
      $this->{'engineName'} = "nhmmer";
      last;
    }
    last if ( $line !~ /^#/ );
  }

  if ( $this->{'version'} eq "" ) {
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

=head2 getParameters()

  Use: my  $wuBlastParamString  = getParameters( );

  Convert object parameters into WUBlast command line parameters.

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

  # TESTING: Setting BG
  #
  my $bgFile = undef;
  if ( 0 ) {
    my $value;
    if ( $value = $this->getMatrix() ) {
      if ( $value =~ /(\d\d)g/ ) {
        my $GCFrac = $1 / 200;
        my $ATFrac = ( 100 - $1 ) / 200;
        my $outputDirName;
        my $currentTime;
        if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
          $outputDirName = $this->getTempDir();
        }
        else {
          $outputDirName = dirname( $this->getQuery() );
        }
        do {
          $currentTime = time();
          $bgFile      = $outputDirName . "/hmmerBGFile-$currentTime-$$.dat";
        } while ( -f $bgFile );
        open BGFILE, ">$bgFile"
            or die "Could not open hmmer bgfile!: $bgFile\n";
        print BGFILE "DNA\n";
        print BGFILE "A $ATFrac\n";
        print BGFILE "C $GCFrac\n";
        print BGFILE "G $GCFrac\n";
        print BGFILE "T $ATFrac\n";
        close BGFILE;
      }
    }

    # END Testing
  }

  # Generate parameter line
  #   Each model should have three cutoffs.  Currently we are
  #   using the Gathering Threshold.
  #
  # GA = Gathering Treshold
  # TC = Trusted Cutoff
  # NC = Noise Cutoff
  my $parameters = " --cut_ga ";

  #
  # Enable restrictions on thread usage. If SearchEngine::setCores
  # is not initialized then nhmmer will run without a "-cpu #"
  # setting.  Without "-cpu #" nhmmer will attempt to use as many
  # cores as are found on the system.  So two simultaneous nhmmer
  # runs on a 4 core machine will attempt to startup 8 threads ( 4
  # per process ).  Note that in it's current state a nhmmer + Dfam
  # thread could easily use up to 1gb of RAM ( 6k position model
  # aligning to a 3k base long window produces a DP of ~18 million
  # cells x 2 tables x 24 bytes per cell ).
  #
  #
  if ( $this->getCores() =~ /\d+/ ) {
    $parameters .= " --cpu " . $this->getCores() . " ";
  }

  # TESTING:
  if ( $bgFile ne "" ) {
    $parameters .= " --bgfile $bgFile ";
  }

  # End TESTING
  my $spanParameter = "";
  my $value;
  if ( ( $value = $this->getSubject() ) ) {
    if ( !-f "$value" ) {
      if ( -s "$value.hmm" ) {
        $parameters .= " $value.hmm";
      }
      else {
        croak $CLASS
            . "::search: Error...hmm subject "
            . "database ($value) does not exist!\n";
      }
    }
    else {
      $parameters .= " $value";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
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

  return ( "$engine $runParameters" );
}

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or                                                                            
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( query=>"file.fa",
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

  # Form the command line
  my $cmdLine = $this->getParameters();

  my $outputDirName;
  if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
    $outputDirName = $this->getTempDir();
  }
  else {
    $outputDirName = dirname( $this->getQuery() );
  }

  # Invoke engine and handle errors
  my $POUTPUT = new FileHandle;
  my $currentTime;
  my $errFile;
  do {
    $currentTime = time();
    $errFile     = $outputDirName . "/hmmerResults-$currentTime-$$.err";
  } while ( -f $errFile );
  my $pid;

  print $CLASS
      . "::search() Invoking search engine as: $cmdLine "
      . " 2>$errFile |\n"
      if ( $this->getDEBUG() );

  $pid = open( $POUTPUT, "$cmdLine 2>$errFile |" );

  my %parseParams = ();
  $parseParams{'debug'} = $this->getDEBUG() if ( $this->getDEBUG() );

  # Create SearchResultCollection object from
  # the engine results.
  my $searchResultsCollection;
  ## Create a debug file
  my $outFile = $outputDirName . "/hmmerResults-$currentTime-$$.out";
  open OUT, ">$outFile";
  while ( <$POUTPUT> ) {
    print OUT $_;
  }
  close OUT;
  close $POUTPUT;
  $parseParams{'searchOutput'} = $outFile;
  $searchResultsCollection = parseOutput( %parseParams );

  my $resultCode = ( $? >> 8 );
  print "HMMER returned a the following result code >$resultCode<\n"
      if ( $this->getDEBUG() );

  ##
  ## DEBUGING for Kaitlin
  ##
#if ( $resultCode )
#{
#  my $lowerCode = ( $? & 255 );
#  print "ERROR returned from nhmmer invocation: 16bit code = $?, high 8bits = $resultCode, lower 8bits = $lowerCode,  error file = $errFile, output file = $outFile\n";
#}

  print $CLASS
      . "::search: "
      . $searchResultsCollection->size()
      . " hits before postprocessing\n"
      if ( $this->getDEBUG() );

  #
  # Postprocess the results
  #
  if ( defined $searchResultsCollection
       && $searchResultsCollection->size() > 0 )
  {

    # The final result collection should be sorted by
    #  queryname and secondarily by query start position.
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
            || $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart();
      }
    );

    #
    # Mask level filtering ( using bitScores )
    #
    my $maskLevel;
    if ( defined( $maskLevel = $this->getMaskLevel() ) && $maskLevel < 101 ) {
      my $strand;
      print "Using: masklevel = $maskLevel\n" if ( $this->getDEBUG() );
      print ""
          . $searchResultsCollection->toString( SearchResult::OutFileFormat )
          . "\n";
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
          if ( $this->getDEBUG() );
      print ""
          . $searchResultsCollection->toString( SearchResult::OutFileFormat )
          . "\n";
    }

    for ( my $i = 0 ; $i < $searchResultsCollection->size() ; $i++ ) {
      my $result = $searchResultsCollection->get( $i );

      #
      # Calculate Kimura divergence using the CpG modification described
      # in SearchResult.pm
      #
      my ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
          $result->calcKimuraDivergence( divCpGMod => 1 );
      $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    }

  }

  unless ( $this->getDEBUG() || $resultCode ) {
    unlink $errFile;
    unlink $outFile;
    return ( $resultCode, $searchResultsCollection );
  }
  else {
    return ( $resultCode, $searchResultsCollection, $outFile, $errFile );
  }

}

##-------------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = HMMERSearchEngine::parseOutput(
                                     searchOutput => $filename|$FH,
                                     [excludeAlignments => 1]
                                                                    );

  Parse the result of a search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $HMMERFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $nameValueParams{'debug'} );
    open $HMMERFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $HMMERFILE = $nameValueParams{'searchOutput'};
  }

  my $callbackFunc = undef;
  if ( defined $nameValueParams{'callback'}
       && ref( $nameValueParams{'callback'} ) == /CODE/ )
  {
    $callbackFunc = $nameValueParams{'callback'};
  }

  my $isHMMSCAN;
  my $inAlignState = 0;
  my $sbjID        = "";
  my $qryID        = "";
  my $score        = 0;
  my $evalue       = 0;
  my $adjScore     = 0;
  my $sbjSeq       = "";
  my $qrySeq       = "";
  my $ppSeq        = "";
  my $sbjOrient    = "";
  my $qryOrient    = "";
  my $sbjStart     = 0;
  my $sbjEnd       = 0;
  my $qryStart     = 0;
  my $qryEnd       = 0;
  my $qryLength    = 0;
  my $sbjLength    = 0;

  my $resultColl = SearchResultCollection->new();

  while ( <$HMMERFILE> ) {

    if ( $inAlignState
         && /^\s*([0-9\*\.]+)\s+PP\s*$/ )
    {
      $ppSeq .= $1;
    }

    #print "processing: $_";

    #
    # Conditions for the end of a hit record:
    #   o Must have seen a score
    #   o Must be in the alignment state (ie. have seen Query: and
    #     Subj: recently)
    #   o Must see something which isn't either "Query:" "Subj:" or " "
    #     or must see the end of the file
    #
    if ( $inAlignState == 2 ) {

      #
      # Reorient if this is a reverse strand
      # hit.
      #
      if ( $qryStart > $qryEnd ) {
        $qryOrient = "C" if ( $qryStart > $qryEnd );
        my $tmp = $qryStart;
        $qryStart = $qryEnd;
        $qryEnd   = $tmp;
      }

      if ( 1 ) {
        my $orientation = "";
        if ( $qryOrient eq "C" ) {

          # Fix the sequence orientation so that it
          # matches the SearchResult.pm convention of
          # the query being in the forward direction.
          $qrySeq = reverse $qrySeq;
          $qrySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $sbjSeq = reverse $sbjSeq;
          $sbjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $ppSeq       = reverse $ppSeq;
          $orientation = "C";
        }

        # HMMER uses "." as the gap character in places
        $qrySeq =~ s/\./-/g;
        $sbjSeq =~ s/\./-/g;

        #
        # Calculate percent divergence
        #           percent insertions
        #           percent deletions
        #
        my %baseFreq = ();
        my $mismatch = 0;
        for ( my $i = 0 ; $i < length( $qrySeq ) ; $i++ ) {
          my $qryBase = substr( $qrySeq, $i, 1 );
          my $sbjBase = substr( $sbjSeq, $i, 1 );
          next if ( $qryBase =~ /-|x/i || $sbjBase =~ /-|x/i );
          $baseFreq{$qryBase}++;
          $mismatch++ if ( $qryBase ne $sbjBase );
        }
        my $percDiv =
            sprintf( "%4.2f", $mismatch * 100 / ( $qryEnd + 1 - $qryStart ) );
        my $qgap = $qrySeq =~ tr/-/-/;
        my $sgap = $sbjSeq =~ tr/-/-/;
        my $percIns =
            sprintf( "%4.2f", $sgap * 100 / ( ( $sbjEnd + 1 ) - $sbjStart ) );
        my $percDel =
            sprintf( "%4.2f", $qgap * 100 / ( ( $qryEnd + 1 ) - $qryStart ) );

        my $result = SearchResult->new(
                                     queryName      => $qryID,
                                     queryStart     => $qryStart,
                                     queryEnd       => $qryEnd,
                                     queryRemaining => ( $qryLength - $qryEnd ),
                                     queryString    => $qrySeq,
                                     subjString     => $sbjSeq,
                                     orientation    => $orientation,
                                     subjName       => $sbjID,
                                     subjStart      => $sbjStart,
                                     subjEnd        => $sbjEnd,
                                     subjRemaining  => ( $sbjLength - $sbjEnd ),
                                     pctDiverge     => $percDiv,
                                     pctInsert      => $percIns,
                                     pctDelete      => $percDel,
                                     matrixName     => "HMM",
                                     score          => int( $score ),
                                     EValue         => $evalue,
                                     PPString       => $ppSeq
        );

        if ( defined $callbackFunc ) {
          $callbackFunc->( $result );
        }
        else {
          $resultColl->add( $result );
        }

        #print "ALIGN: " .
        #$result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";

      }

      if ( $isHMMSCAN ) {
        $sbjID = "";
      }
      else {
        $qryID     = "";
        $qryLength = 0;
      }
      $score        = "";
      $evalue       = "";
      $sbjSeq       = "";
      $qrySeq       = "";
      $qryOrient    = "";
      $ppSeq        = "";
      $sbjStart     = 0;
      $sbjEnd       = 0;
      $qryStart     = 0;
      $qryEnd       = 0;
      $inAlignState = 0;
    }

    if ( /^Query:\s*(\S+)\s+\[([LM])=(\d+)\].*/ ) {

      # The two programs report the Query differently.  The designation
      # for the query length is used to distinguish between the two:
      #  nhmmscan
      #     Query:       Human:44186-44485-c_b1s352i14  [L=300]
      #  nhmmer
      #     Query:       7SLRNA#SINE/Alu  [M=320]
      #
      # Remember that RepeatMasker uses the "query" label to designate
      # the sequence being searched.  HMMER designates the query as the
      # HMM, HMMSCAN designates the subj as the HMM.  So depending on the
      # program we need to reverse the labels here.
      if ( $2 eq "L" ) {
        $qryID     = $1;
        $qryLength = $3;
        $isHMMSCAN = 1;
      }
      else {
        $sbjID     = $1;
        $sbjLength = $3;
        $isHMMSCAN = 0;
      }
    }

    # Start of a search result
    if ( /^>>\s+(\S+)/ ) {

      # Remember that RepeatMasker uses the "query" label to designate
      # the sequence being searched.  HMMER designates the query as the
      # HMM, HMMSCAN designates the subj as the HMM.  So depending on the
      # program we need to reverse the labels here.
      if ( $isHMMSCAN ) {
        $sbjID = $1;
      }
      else {
        $qryID = $1;
      }
    }

    # Typical header fields:
    #   ?/!  score  bias Evalue  hmm_from  hmm_to
    #   ali_from  ali_to env_from env_to sq_len acc
    if ( /^\s+([\?\!]\s+[-\d\.]+.*)/ ) {
      my $null;
      my @hdrFields = split( /\s+/, $1 );
      if ( @hdrFields == 15 ) {
        if ( $isHMMSCAN ) {
          (
            $null,   $score, $null,     $evalue, $sbjStart,
            $sbjEnd, $null,  $qryStart, $qryEnd, $null,
            $null,   $null,  $null,     $sbjLength
              )
              = @hdrFields;
        }
        else {
          (
            $null,   $score, $null,     $evalue, $sbjStart,
            $sbjEnd, $null,  $qryStart, $qryEnd, $null,
            $null,   $null,  $null,     $qryLength
              )
              = @hdrFields;
        }
      }
      else {
        warn "Can't parse: $_\n";
      }
    }

    # Start of an alignment
    if ( /^\s+Alignment:/ ) {
      die "$CLASS: Error - Alignment found but queryEnd wasn't recorded yet!"
          if ( $qryEnd == 0 );
      $inAlignState = 1;
    }

    #
    # Alignment lines
    #
    if ( $inAlignState
         && /^\s+(\S+)\s+(\d+|-)\s+([acgtnbdhvrykmsw\.\-]+)\s+(\d+|-)/i )
    {

      # Query or Subject
      if ( $1 eq $qryID ) {
        $qrySeq .= uc( $3 );
        if ( $4 == $qryEnd ) {
          $inAlignState = 2;
        }
      }
      elsif ( $1 eq $sbjID ) {
        $sbjSeq .= uc( $3 );
      }
      else {
        croak "ERROR: Cannot determine alignment sequence:\n"
            . "       qryID = $qryID, sbjID = $sbjID\n"
            . "       Line ID = $1\n";
      }
    }

  }
  close $HMMERFILE;

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
