#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) rmToTrackHub.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
# ChangeLog
#
#     $Log$
#
###############################################################################
#

=head1 NAME

rmToTrackHub.pl - generate track hub files for a UCSC bigRmsk track 

=head1 SYNOPSIS

  rmToTrackHub.pl [-version]
                  [-only_beds][-force]
                  -genome <ucsc genome> 
                  -out <repeatmasker *.out[.gz]> 
                  -hubname <track_hub_name>
                  [-align <*.align[.gz]> ]

  e.g:

  ./rmToTrackHub.pl -hubname myHub -genome hg38 -out hg38.fa.out -align hg38.fa.align

=head1 DESCRIPTION

Generate the necessary files for visualizing the output of RepeatMasker on
a UCSC hosted genome using the new BigRmsk TE visualization track.  This
tool will query the UCSC site for details on the sizes of the hosted sequences,
build the bigBed files for the visualization and generate track hub template
files (trackDB.txt, hub.txt and genomes.txt).

The options are:

=over 4

=item -version

Displays the version of the program

=item -force

Force overwriting files in the hub directory if they already exist.

=item -only_beds

Only generate the bigBed files for the RepeatMasker output files.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2022-2024 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use IO::File;
use FindBin;
use File::Basename;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use CrossmatchSearchEngine;
use SearchResultCollection;
use SearchResult;

my $Version = "0.3";
my $oldFormat = 0;  # Must be global for callback
my @alignMap  = (); # Must be global for callback
my $alignPos  = 0;  # dito
my $DEBUG     = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-genome=s',
                    '-hubname=s',
                    '-out=s',
                    '-force|f',
                    '-only_beds',
                    '-align=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

# Dependencies need to be in path (for now)
my $retVal = `whereis bedToBigBed`;
my $BEDTOBIGBED_PRGM = "";
if ( $retVal =~ /^bedToBigBed:\s+(\S+)/ ){
  $BEDTOBIGBED_PRGM = $1;
}else {
  print "\n\nMissing bedToBigBed UCSC tool.  Make sure the tool is in your path.\n\n";
  usage();
}

if ( ! exists $options{'out'} || ! -s $options{'out'} ) {
  print "\n\nMissing or empty RepeatMasker *.out file\n\n";
  usage();
}
my $outFile = $options{'out'};

if ( ! exists $options{'hubname'} || $options{'hubname'} =~ /\s+/ ) {
  print "\n\nMissing hubname parameter or it contains spaces!\n\n";
  usage();
}
my $hubname = $options{'hubname'};

if ( ! exists $options{'genome'} || $options{'genome'} =~ /\s+/ ) {
  print "\n\nMissing genome parameter or it contains spaces!\n\n";
  usage();
}
my $genome = $options{'genome'};

my $alignFile;
$alignFile = $options{'align'} if ( $options{'align'} );

print "\n";
print "#\n";
print "# rmToTrackHub.pl\n";
print "#\n";

print "# Obtaining sequences sizes for $genome from UCSC....\n";
&getChromSizes($genome);
if ( ! -s "$genome.chrom.sizes" ) {
  print "Failed to download $genome.chrom.sizes file for $genome from UCSC!\n";
  exit(1);
}

if ( ! -d $hubname ) {
  mkdir($hubname);
}elsif ( ! $options{'f'} ) {
  print "The directory $hubname already exists.  If you would like to overwrite it use the -f automatically use the -f flag.\n";
}

print "# Building TSV files...\n";
my ($joinTSVFile, $alignTSVFile) =  &buildTSVFiles($outFile,$alignFile);

$outFile = basename($outFile);
if ( $alignFile ) {
  $alignFile = basename($alignFile);
}
print "# Building bigBed files...\n";
# Make out file bigBed
my $cmd = "$BEDTOBIGBED_PRGM -tab -as=$FindBin::RealBin/bigRmskBed.as -type=bed9+5 $joinTSVFile $genome.chrom.sizes $hubname/$outFile.bb";
system($cmd);

# Make align file bigBed
if ( $alignFile ) {
$cmd = "$BEDTOBIGBED_PRGM -tab -as=$FindBin::RealBin/bigRmskAlignBed.as -type=bed3+14 $alignTSVFile $genome.chrom.sizes $hubname/$alignFile.bb";
system($cmd);
}

unless ( $options{'only_beds'} ) {
  print "# Making trackDB.txt\n";
  open OUT,">$hubname/trackDB.txt" or die "Could not create trackDB.txt file!\n";
  print OUT "track myTrackName\n";
  print OUT "superTrack on show\n";
  print OUT "shortLabel myShortLabel\n";
  print OUT "longLabel My Long Label\n\n";
  print OUT "\ttrack subTrackLabel\n";
  print OUT "\tparent myTrackName\n";
  print OUT "\ttype bigRmsk 9 +\n";
  print OUT "\tvisibility hide\n";
  print OUT "\tbigDataUrl $outFile.bb\n";
  if ( $alignFile ) {
    print OUT "\txrefDataUrl $alignFile.bb\n";
  }
  print OUT "\tshortLabel Subb Track Short Label\n";
  print OUT "\tlongLabel  Sub Track Long Label with up to 80 characters...\n";
  print OUT "\thtml subTrackDescription.html\n";
  print OUT "\tpriority 30\n";
  close OUT;

  print "# Making subTrackDescription.html\n";
  open OUT,">$hubname/subTrackDescription.html" or die "Could not create index.html file!\n";
  print OUT "<html><h1>My Sub Track Description</h1></html>\n";
  close OUT;

  print "# Making hub.txt\n";
  open OUT,">$hubname/hub.txt" or die "Could not create hub.txt file!\n";
  print OUT "hub myHubWithoutSpaces\n";
  print OUT "shortLabel My Hubs Name\n";
  print OUT "longLabel My Hub Long Label up to 80 characters vs shortLabel limited to 17\n";
  print OUT "genomesFile genomes.txt\n";
  print OUT "email myEmail\@myDomain\n";
  print OUT "descriptionUrl index.html\n";
  close OUT;

  print "# Making genomes.txt\n";
  open OUT,">$hubname/genomes.txt" or die "Could not create genomes.txt file!\n";
  print OUT "genome $genome\n";
  print OUT "trackDb trackDB.txt\n";
  close OUT;

  print "# Making index.html\n";
  open OUT,">$hubname/index.html" or die "Could not create index.html file!\n";
  print OUT "<html><h1>My Track Hub</h1></html>\n";
  close OUT;


  print "\nDone building trackHub!  Please move/copy the directory '$hubname/'\n";
  print "to a webserver and point the UCSC genome browser at your trackhub using\n";
  print "this page: https://genome.ucsc.edu/cgi-bin/hgHubConnect.  Your trackhub\n";
  print "URL will look something like this: 'http:://yourdomain/$hubname/hub.txt'\n";
  print "For more help with UCSC trackhubs and settings please see the help page\n";
  print "at: https://genome.ucsc.edu/goldenPath/help/hubQuickStart.html\n\n";
}else {
  print "\nDone building bigBED file(s).\n";
}

exit;

######################## S U B R O U T I N E S ############################

sub buildTSVFiles  {
  my $outFile = shift;
  my $alignFile = shift;

  # need to set locale to C for output sort to work correctly.
  $ENV{'LC_COLLATE'} = "C";
  $ENV{'LC_ALL'} = "C";

  if ( defined $alignFile ) {
    my $alignTSVFile = basename($alignFile);
    $alignTSVFile =~ s/\.align.*/.align.tsv/;
    open ALIGNTSV, "| sort -k1,1 -k2,2n >$alignTSVFile" or die "Could not open $alignTSVFile\n";
  
    ## Determine if the alignment file was lifted correctly.
    ##   - Old alignment files didn't have the ID field lifted
    if ( $alignFile =~ /(\S+)\.gz/ ) {
      open IN, "gunzip -c $alignFile |"
          or die "Could not open compressed alignment file: $alignFile!\n";
    }
    else {
      open IN, "<$alignFile"
          or die "Could not open alignment file: $alignFile!\n";
    }
    my $lineIdx    = 1;
    my $prevIDLine = 0;
  
    print "Checking $alignFile alignment file...\n";
    while ( <IN> ) {
      if ( /^\s*\d+\s+\d+\.\d+\s+/ ) {
        my ( $ID ) = ( /(\d+)\s*$/ );
        if ( $ID == 1 ) {
          if ( $prevIDLine && ( $lineIdx - $prevIDLine ) > 500 ) {
            print "Hmmm...lineIdx = $lineIdx ID = $ID and prevIDLine = $prevIDLine\n";
            print $_;
            # Not likely
            $oldFormat = 1;
            last;
          }
          else {
            $prevIDLine = $lineIdx;
          }
        }
        $lineIdx++;
      }
    }
    close IN;
  
    if ( $oldFormat ) {
      warn "\n    The "
          . $options{'align'}
          . " file has repeating identifiers ( ID field ).\n"
          . "    Either this is an older UCSC run or lifting of the ID field\n"
          . "    failed.  Making a last ditch effort to re-link *.align\n"
          . "    annotations with *.out annotations using bedTools.\n"
          . "    This is an imperfect process and some data loss or\n"
          . "    mismatches may occur.  You may want to skip building/loading\n"
          . "    the *.align data for this run.\n\n";
  
      print "Creating temporary bed files...\n";
      my $outBedFile = basename($outFile);
      $outBedFile =~ s/\.(rm)?out.*/\.out.bed/;
      _rmToBedFile( rmFile => $outFile, bedFile => $outBedFile );
      my $alignBedFile = basename($alignFile);
      $alignBedFile =~ s/\.(rm)?align.*/\.align.bed/;
      _rmToBedFile( rmFile => $alignFile, bedFile => $alignBedFile );
  
      print "Comparing ranges using bedTools...\n";
      open DATA, "intersectBed -f 0.9 -loj -a $alignBedFile -b $outBedFile |"
          or die "Could not open external program 'intersectBed' "
          . "from the BedTools package!";
  
      ## Inexplicably there are exact duplicate lines in entries in
      ## the UCSC *.align files: ie mm10: chr1:26251593-26252732
      ## Class check won't work as there are too many good but not
      ## exact class joins.
      ##
      ## TODO: Another way to do this would be to consider the job
      ## fragmentation size and how RMPart batches were built.  It
      ## could be this will be far more accurate way of detecting ID resets.
      my $unmappable      = 0;
      my $totalAlignments = 0;
      while ( <DATA> ) {
        my @fields = split( /\s+/, $_ );
  
        my $alignPos = $fields[ 18 ];
        $alignPos =~ s/[\[\]]//g;
  
        if ( $fields[ 19 ] eq "." ) {
          $alignMap[ $alignPos ] = -1;
          $unmappable++;
          next;
        }
        my $outID = $fields[ 36 ];
  
        if ( $alignPos <= $#alignMap ) {
          if ( $alignMap[ $alignPos ] == -1 ) {
          }
          elsif ( $alignMap[ $alignPos ] != $outID ) {
            $alignMap[ $alignPos ] = -1;
            $unmappable++;
          }
        }
        else {
          $alignMap[ $alignPos ] = $outID;
        }
      }
      close DATA;
      unlink( $alignBedFile );
      unlink( $outBedFile );
      print "Total unmappable alignments = $unmappable\n";
  
    } # if oldformat
    print "Writing out $alignTSVFile...\n";
    $alignPos = 0;
    my $ALIGN = IO::File->new();
    if ( $options{'align'} =~ /(\S+)\.gz/ ) {
      open $ALIGN, "gunzip -c $options{'align'} |"
          or die "Could not open compressed alignment file: $options{'align'}!\n";
    }
    else {
      open $ALIGN, "<$options{'align'}"
          or die "Could not open alignment file: $options{'align'}!\n";
    }
  
    #open $ALIGN, "gunzip -c $options{'align'} |"
    #    or die "Could not open $options{'align'} file!\n";
    CrossmatchSearchEngine::parseOutput( searchOutput => $ALIGN,
                                         callback     => \&procAlignResult );
  
    close ALIGNTSV;
  }
  
  # *.out files are typically lifted and IDs are fixed by ucsc
  #
  
  #track db='hg19' name='rmskFamilies' priority='70' description='Repeat Masker family graph' visibility=pack itemRgb='On' exonArrows='on'
  #chr21	33031600	33032600	AluXyz	0	+	33031600	33032400	225,124,122	2	800,1	0,999
  #chr21	33034600	33040000	AluSz,_Family_Alu,_20%_Kimura	0	+	33034800	33039000	50,93,129	4	1,400,300,1	0,200,3099,5399
  #chr21	33035000	33037500	AluAbc	0	-	33035200	33037500	190,228,116	3	1,400,500	0,300,2000
  #chr21	33035400	33037400	AluDef	0	+	33035400	33037400	225,124,122	3	1,1450,1	0,280,1999
  #chr21	33038100	33039100	AluGhi	0	+	33038100	33039100	252,152,72	2	1,500	0,500
  #chr21	33039300	33041000	AluJkl	0	-	33039300	33041000	100,153,200	1	1700	0
  
  #print "track db='hg19' name='rmskFamilies' priority='70' "
  #   .  "description='Repeat Masker family graph' visibility=pack "
  #   .  "itemRgb='On' exonArrows='on'\n";
  
  my $joinTSVFile = basename($outFile);
  $joinTSVFile =~ s/\.(rm)?out.*/\.join.tsv/;

  open JOINTSV, "| sort -k1,1 -k2,2n >$joinTSVFile"
      or die "Could not open $joinTSVFile for writing!\n";
  
  my @bedKeys = (
                  "chrom",   "chromStart", "chromEnd",   "name",
                  "score",   "strand",     "thickStart", "thickEnd",
                  "itemRgb", "blockCount", "blockSizes", "blockStarts",
                  "id", "description"
  );
  
  my %dataHash   = ();
  my $maxIDDelta =
      1000;    ## absurdly high because of a bug in versions prior to 4.0
  my %ids         = ();
  my $idx         = 1;
  my $highestID   = 0;
  my $prevSeqName = "";
  
  if ( $outFile =~ /(\S+)\.gz/ ) {
    open IN, "gunzip -c $outFile |"
        or die "Could not open compressed alignment file: $outFile!\n";
  }
  else {
    open IN, "<$outFile"
        or die "Could not open alignment file: $outFile!\n";
  }
  
  print "Writing out $joinTSVFile...\n";
  while ( <IN> ) {
    if ( /^\s*(\d+\s+\d+\.\d+.*)/ ) {
      my $outLine = $1;
      my $leftUnalignedSize;
      my $rightUnalignedSize;
      my $rec    = $1;
      my @fields = split( /\s+/, $rec );
      my $outStr = join( "\t", @fields );
  
      #$outStr =~ s/[\(\)]//g;
      #print OUTTSV "$outStr\n";
      my $score = $fields[ 0 ];
      my $div   = $fields[ 1 ];
      my $pDel  = $fields[ 2 ];
      my $pIns  = $fields[ 3 ];
      my $seq   = $fields[ 4 ];
  
      # Output files use 0-based, half open coordinates
      my $qStart     = $fields[ 5 ] - 1;
      my $qEnd       = $fields[ 6 ];
      my $qRemaining = $fields[ 7 ];
      $qRemaining =~ s/[\(\)]//g;
      my $orient  = $fields[ 8 ];
      my $sbjName = $fields[ 9 ];
      my $class   = $fields[ 10 ];
      my ( $type, $noop, $subtype ) = ( $class =~ /([^\/]+)(\/(\S+))?/ );
      my $sbjStart     = $fields[ 11 ];
      my $sbjEnd       = $fields[ 12 ];
      my $sbjRemaining = $fields[ 13 ];
      my $id           = $fields[ 14 ];
      if ( $id !~ /\d+/ || $id < 0 ) {
        warn "Identifier in out file is invalid for: $_   - Setting value to 999999999\n";
        $id = 999999999;
      }
  
      if ( $orient eq "C" ) {
        $sbjRemaining = $fields[ 11 ];
        $sbjEnd       = $fields[ 12 ];
        $sbjStart     = $fields[ 13 ];
      }
      $sbjRemaining =~ s/[\(\)]//g;
      # Work around some issues with ProcessRepeats
      $sbjRemaining = 0 if ( $sbjRemaining < 0 );
      $sbjStart = 1 if ( $sbjStart < 1 );
  
      if ( $seq ne $prevSeqName ) {
        %ids = ();
        foreach my $idKey ( sort { $a <=> $b } keys( %dataHash ) ) {
          if ( $dataHash{$idKey}->{'chromStart'} > $dataHash{$idKey}->{'thickStart'} ||
               $dataHash{$idKey}->{'chromEnd'} < $dataHash{$idKey}->{'thickEnd'} )
          {
            warn "WARNING chromStart/End does not bound thickStart/End: " . Dumper($dataHash{$idKey}) . "\n";
          }
          foreach my $bedKey ( @bedKeys ) {
            if ( $bedKey eq "blockSizes" || $bedKey eq "blockStarts" ) {
              print JOINTSV ""
                  . join( ",", @{ $dataHash{$idKey}->{$bedKey} } ) . "\t";
            }
            elsif ( $bedKey eq "score" ) {
            my $avgDiv = sprintf(
                                      "%0.1f",
                                  (
                                    $dataHash{$idKey}->{$bedKey} /
                                        $dataHash{$idKey}->{'divCount'}
                                  )
              ) * 10;
              if ( $avgDiv > 1000 ) {
                warn "WARNING:  $dataHash{$idKey}->{'id'} produced an div score of $avgDiv setting to 100%.\n";
                $avgDiv = 1000;
              }
              print JOINTSV "$avgDiv\t";
            }
            elsif ( $bedKey eq "description" ) {
              print JOINTSV "$dataHash{$idKey}->{$bedKey}\n";
            }
            else {
              print JOINTSV "$dataHash{$idKey}->{$bedKey}\t";
            }
          }
        }
        %dataHash = ();
      } # if ( $seq ne $prevSeqName )
  
      ## ID Checks
      if ( exists $ids{$id} && $ids{$id} < ( $idx - $maxIDDelta ) ) {
        warn "Found an ID renumbering candidate: $seq: id=$id \@ $idx and "
            . $ids{$id} . " ("
            . ( $idx - $ids{$id} )
            . ") $_\n";
      }
      if ( $id > $highestID ) {
        if ( $id - $highestID > 1 ) {
          warn "ID $id was skipped in this input file!\n";
        }
        $highestID = $id;
      }
      $ids{$id} = $idx;
      $idx++;
      $prevSeqName = $seq;
  
      if ( exists $dataHash{$id} ) {
  
        #  Hash key exists: append to data
        my $entry = $dataHash{$id};
        if ( $orient eq "+" ) {
          $entry->{'chromEnd'}   = $qEnd + $sbjRemaining;
          $leftUnalignedSize     = $sbjStart - $entry->{'prevSbjPos'} - 1;
          $rightUnalignedSize    = $sbjRemaining;
          $entry->{'prevSbjPos'} = $sbjEnd;
        }
        else {
          $entry->{'chromEnd'}   = $qEnd + $sbjStart - 1;
          $leftUnalignedSize     = $entry->{'prevSbjPos'} - $sbjEnd - 1;
          $rightUnalignedSize    = $sbjStart - 1;
          $entry->{'prevSbjPos'} = $sbjStart;
        }
  
        if ( $entry->{'chromEnd'} > ( $qEnd + $qRemaining ) ) {
          $entry->{'chromEnd'} = $qEnd + $qRemaining;
        }
        $entry->{'thickEnd'} = $qEnd;
        $entry->{'score'} += $div;
        $entry->{'divCount'}++;
  
        $entry->{'blockSizes'}->[ $entry->{'blockCount'} - 1 ] =
            $leftUnalignedSize;
        push @{ $entry->{'blockSizes'} },  ( $qEnd - $qStart );
        push @{ $entry->{'blockStarts'} }, ( $qStart - $entry->{'chromStart'} );
        $entry->{'blockCount'}++;
        push @{ $entry->{'blockSizes'} },  $rightUnalignedSize;
        push @{ $entry->{'blockStarts'} }, -1;
        $entry->{'blockCount'}++;
        $outLine =~ s/[\n\r]+//g;
        $outLine =~ s/\s+/ /g;
        $entry->{'description'} .= ",".$outLine;
      }
      else {
  
        #  Create a new hash entry
        my $entry = {};
        $entry->{'chrom'} = $seq;
        if ( $orient eq "+" ) {
          $entry->{'chromStart'} = $qStart - $sbjStart;
          $entry->{'chromEnd'}   = $qEnd + $sbjRemaining;
          $leftUnalignedSize     = $sbjStart - 1;
          $rightUnalignedSize    = $sbjRemaining;
          $entry->{'prevSbjPos'} = $sbjEnd;
        }
        else {
          $entry->{'chromStart'} = $qStart - $sbjRemaining;
          $entry->{'chromEnd'}   = $qEnd + $sbjStart - 1;
          $leftUnalignedSize     = $sbjRemaining;
          $rightUnalignedSize    = $sbjStart - 1;
          $entry->{'prevSbjPos'} = $sbjStart;
        }
        if ( $entry->{'chromStart'} < 0 ) {
          $entry->{'chromStart'} = 0;
        }
        if ( $entry->{'chromEnd'} > ( $qEnd + $qRemaining ) ) {
          $entry->{'chromEnd'} = $qEnd + $qRemaining;
        }
        $entry->{'name'}  = "$sbjName#$class";
        $entry->{'score'} = $div;
        $entry->{'divCount'}++;
        $entry->{'strand'}     = "+";
        $entry->{'strand'}     = "-" if ( $orient eq "C" );
        $entry->{'thickStart'} = $qStart;
        $entry->{'thickEnd'}   = $qEnd;
        $entry->{'itemRgb'}    = "0";
        $entry->{'blockSizes'} = [];
        $entry->{'blockStart'} = [];
        push @{ $entry->{'blockSizes'} },  $leftUnalignedSize;
        push @{ $entry->{'blockStarts'} }, -1;
        $entry->{'blockCount'}++;
        push @{ $entry->{'blockSizes'} },  ( $qEnd - $qStart );
        push @{ $entry->{'blockStarts'} }, ( $qStart - $entry->{'chromStart'} );
        $entry->{'blockCount'}++;
        push @{ $entry->{'blockSizes'} },  $rightUnalignedSize;
        push @{ $entry->{'blockStarts'} }, -1;
        $entry->{'blockCount'}++;
        $entry->{'id'} = $id;
        $outLine =~ s/[\n\r]+//g;
        $outLine =~ s/\s+/ /g;
        $entry->{'description'} = $outLine;
        $dataHash{$id} = $entry;
      }
    }
  }
  ## Trailing case...TODO: make this a subroutine to reduce the code duplication
      %ids = ();
      foreach my $idKey ( sort { $a <=> $b } keys( %dataHash ) ) {
        if ( $dataHash{$idKey}->{'chromStart'} > $dataHash{$idKey}->{'thickStart'} ||
             $dataHash{$idKey}->{'chromEnd'} < $dataHash{$idKey}->{'thickEnd'} )
        {
          warn "WARNING chromStart/End does not bound thickStart/End: " . Dumper($dataHash{$idKey}) . "\n";
        }
        foreach my $bedKey ( @bedKeys ) {
          if ( $bedKey eq "blockSizes" || $bedKey eq "blockStarts" ) {
            print JOINTSV ""
                . join( ",", @{ $dataHash{$idKey}->{$bedKey} } ) . "\t";
          }
          elsif ( $bedKey eq "score" ) {
            my $avgDiv = sprintf(
                                  "%0.1f",
                                  (
                                    $dataHash{$idKey}->{$bedKey} /
                                        $dataHash{$idKey}->{'divCount'}
                                  )
            ) * 10;
            if ( $avgDiv > 1000 ) {
              warn "WARNING:  $dataHash{$idKey}->{'id'} produced an div score of $avgDiv setting to 100%.\n";
              $avgDiv = 1000;
            }
            print JOINTSV "$avgDiv\t";
          }
          elsif ( $bedKey eq "description" ) {
            print JOINTSV "$dataHash{$idKey}->{$bedKey}\n";
          }
          else {
            print JOINTSV "$dataHash{$idKey}->{$bedKey}\t";
          }
        }
      }
      %dataHash = ();
  ## end trailing case
  close JOINTSV;
  return ( $joinTSVFile, $alignTSVFile);
}


##
## Process an alignment record
##
sub procAlignResult {
  my $result = shift;
  my $id     = $result->getId();

  if ( $oldFormat ) {
    if ( defined $alignMap[ $alignPos ] && $alignMap[ $alignPos ] >= 1 ) {
      $result->setId( $alignMap[ $alignPos ] );
    }
    else {
      $alignPos++;
      return;
    }
  }

  #print "". $result->toStringFormatted( SearchResult::AlignWithQuerySeq );
  if (    $result->getSubjType() eq ""
       && $result->getSubjName() =~ /^(\S+)\#(\S+)/ )
  {
    $result->setSubjName( $1 );
    $result->setSubjType( $2 );
  }

  my $rec = $result->toStringFormatted( SearchResult::CompressedAlignCSV );
  my @recFields = split( ',', $rec );
  my $sequence  = $recFields[ 16 ];

  my $cRec =
      $result->getQueryName() . "\t"
      . ( $result->getQueryStart() - 1 ) . "\t" # zero-based, half open
      . $result->getQueryEnd() . "\t"           # zero-based, half open
      . $result->getQueryRemaining() . "\t"
      . $result->getScore() . "\t"
      . $result->getPctDiverge() . "\t"
      . $result->getPctDelete() . "\t"
      . $result->getPctInsert() . "\t";

  if ( $result->getOrientation =~ /C|c/ ) {
    $cRec .= "-\t";
  }
  else {
    $cRec .= "+\t";
  }

  $cRec .= $result->getSubjName() . "\t";

  if ( $result->getSubjType() =~ /(\S+)\/(\S+)/ ) {
    $cRec .= $1 . "\t" . $2 . "\t";
  }
  else {
    $cRec .= $result->getSubjType() . "\t\t";
  }

  $cRec .=
        $result->getSubjStart() . "\t"
      . $result->getSubjEnd() . "\t"
      . $result->getSubjRemaining() . "\t";

  $cRec .= $result->getId() . "\t" . $sequence;

  print ALIGNTSV "$cRec\n";

  $alignPos++;
}


##-------------------------------------------------------------------------##
## Use: _rmToBedFile( rmFile => value, bedFile => value );
##
##      $rmFile       : A parameter to the method
##      $bedFile       : A parameter to the method
##
##    Convert an existing RepeatMasker *.out|*.align file to a bed file.
##
##  Returns
##      Nothing
##
##-------------------------------------------------------------------------##
sub _rmToBedFile {
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " . @{ [ %parameters ] } . "): Called\n" if ( $DEBUG );

  if ( defined $parameters{'rmFile'} && -s $parameters{'rmFile'} ) {
    if ( $parameters{'rmFile'} =~ /(\S+)\.gz/ ) {
      open RM, "gunzip -c $parameters{'rmFile'} |"
          or die "$subroutine: Could not open compressed "
          . "alignment file: $parameters{'rmFile'}!\n";
    }
    else {
      open RM, "<$parameters{'rmFile'}"
          or die "$subroutine: Could not open alignment "
          . "file: $parameters{'rmFile'}!\n";
    }

    my $recIndex = 0;
    if ( defined $parameters{'bedFile'} ) {
      open BED, "| sort -k1,1 -k2,2n >$parameters{'bedFile'}"
          or die
          "$subroutine: Could not open output file $parameters{'bedFile'}!\n";
      while ( <RM> ) {
        if ( /^\s*(\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+.*)/ ) {
          my $rec = $1;
          my @fields = split( /\s+/, $rec );
          $rec =~ s/[\n\r]+//g;
          $rec =~ s/\s+/ /g;

          # Zero based half open coords
          print BED "$fields[4]\t"
              . ( $fields[ 5 ] - 1 )
              . "\t$fields[6]\t$rec [$recIndex]\n";
          $recIndex++;
        }
      }
      close RM;
      close BED;
    }
    else {
      die "$subroutine: Missing -bedFile parameter!\n";
    }
  }
  else {
    die "$subroutine: Missing -rmFile parameter or file is empty!\n";
  }
}

sub getChromSizes {
  my $DB = shift;

  my $url="http://hgdownload.soe.ucsc.edu/goldenPath/$DB/bigZips/$DB.chrom.sizes";
  if ( $DB =~ /hg19|hg38/i ) {
    $url="http://hgdownload.soe.ucsc.edu/goldenPath/$DB/bigZips/latest/$DB.chrom.sizes";
  }

  my $retVal = `whereis curl`;
  if ( $retVal =~ /^curl:\s(\S+)/ ) {
    system("curl $url -s > $DB.chrom.sizes 2> /dev/null"); 
  }else {
    $retVal = `whereis wget`;
    if ( $retVal =~ /^wget:\s(\S+)/ ) { 
      system("wget $url -O $DB.chrom.sizes 2> /dev/null");
    }else {
      die "Could not find curl or wget in the users path!\n";
    }
  }
}
