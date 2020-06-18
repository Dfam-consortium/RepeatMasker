#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) EMBL.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A simple EMBL object.  This is a simple datastructure
##      for encapsulating the EMBL formatted sequence data and
##      annotations.  **THIS IS NOT A COMPLETE IMPLEMENTATION OF
##      THE EMBL FORMAT** but rather restricted to the features
##      currently used by RepeatMasker and Repbase. In addition
##      this is an implementation of the ArrayList class which enables
##      EMBL datasets to be concatenated etc..
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2017 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
#    The internal data storage is defined below.  Data accessor
#    methods have not been built ( but should be before external
#    access is allowed ).  Currently the file's records ( defined
#    below ) are stored in an ArrayList:
#
#    bless ArrayList( { 'EMBL_elements'} => [ {#EMBL Record}, {#EMBL Record}, .. ]},
#          'EMBL' );
#
#    RepeatMasker/EMBL Record                            EMBL Tag
#    ------------------------------------------------    --------
#    { 'id' => blah,                                     ID
#      'dataclass' => blah,                              ID
#      'molecule' => blah,                               ID
#      'divisions' => [ blah, blah, blah ],              ID
#      'accessions' => [ blah, blah, blah ],             AC
#      'dateCreated' => blah,                            DT
#      'dateUpdated' => blah,                            DT
#      'comments' => [ blah, blah, blah, blah ],         CC
#      'keywords' => [ blah, blah, blah, blah ],         KW
#      'description' => blah,                            DE
#      'speciesName' => blah,                            OS
#      'commonSpeciesName' => blah,                      OS
#      'classification' => [ blah, blah, blah, blah ],   OC
#      'references' => [ { 'number' = > blah,            RN
#                          'comment' => blah,            RC
#                          'position' => blah,           RP
#                          'xref' => blah,               RX
#                          'authors' => blah,            RA
#                          'title' => blah,              RT
#                          'location' => blah } ],       RL
#       'xref' => blah,                                  DR
#       'length' => blah,                                ID
#       'composition' => { 'A' => blah,                  SQ
#                          'C' => blah,                  SQ
#                          'G' => blah,                  SQ
#                          'T' => blah,                  SQ
#                          'other' => 'blah' },          SQ
#       'sequence' => blah };                            SQ
#
#******************************************************************************
#

=head1 NAME

EMBL


=head1 SYNOPSIS

use EMBL

Usage:

    $DB = EMBL->new( fileName => "humrep.ref" );

=head1 DESCRIPTION

A simple EMBL object.  This is a simple datastructure
for encapsulating the EMBL formatted sequence data and
annotations.  **THIS IS NOT A COMPLETE IMPLEMENTATION OF
THE EMBL FORMAT** but rather restricted to the features
currently used by RepeatMasker and Repbase.  In addition 
this is an implementation of the ArrayList class which enables 
EMBL datasets to be concatenated etc..

Caveats:
   - Currently no attempt is made to store the ordering
     of tag groups.  The output ordering is fixed and
     may not reflect the original input ordering.  Ie.
     it appears there is a convention of using a "CC"
     tag following an "OS" tag to store the common name
     of a speciess:

         OS   Anopheles gambiae
         CC   African malaria mosquito.

     The EMBL Standard requires that this information be
     placed on the "OS" line and surrounded by parantheses.
     The reader does not currently differentiate between
     "CC" locations in a record.  Therefore this will be
     reformatted upon output as:

         OS   Anopheles gambiae
         OC   ..
         ..
         XX
         CC   African malaria mosquito.

   - The line breaks/spacing in comment blocks are potentially
     important.  This object attempts to keep save the line
     break and spacing information for several types of
     data including comment tags.  It will be the responsibility
     of the data user to verify line width and reformat paragraphs
     as the data is manipulated.

=head1 INSTANCE METHODS

=cut

package EMBL;
use strict;
use PubRef;
use RepeatRecord;
use Carp;
use Data::Dumper;
use Text::Wrap;
use ArrayList;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(ArrayList Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = "0.1";
my $CLASS   = "EMBL";

#
# EMBL Tags used by Repbase
#
my %validTagBits = (
                     'ID' => 1,
                     'CC' => 2,
                     'AC' => 4,
                     'DT' => 8,
                     'DE' => 16,
                     'KW' => 32,
                     'OS' => 64,
                     'OC' => 128,
                     'RN' => 256,
                     'RP' => 512,
                     'RC' => 1024,
                     'RA' => 2048,
                     'RT' => 4096,
                     'RL' => 8192,
                     'DR' => 16384,
                     'SQ' => 32768,
                     'RX' => 65536,
                     'FT' => 131072,
                     'RG' => 262144,
                     'FH' => 524288,
                     'NM' => 1048576
);

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  my $this = $class->SUPER::new( @_ );

  my %nameValuePairs = @_;

  # Bless this hash in the name of the father, the son...
  #bless $this, $class;

  if ( %nameValuePairs ) {
    my %validKeys = ( fileName => 1 );
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      unless ( exists $validKeys{$name} ) {
        croak "$CLASS"
            . "::new: Error $name is not a valid "
            . "parameter to the constructor!\n\n"
            . "Usage: \$DB = $CLASS->new();\n"
            . "       \$DB = $CLASS->new( fileName=>\"humrep.ref\" );\n";
      }
      $this->{$name} = $value;
    }
  }

  if ( exists $this->{'fileName'} ) {
    _parseFromFile( $this, $this->{'fileName'} );
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 Use: $DB->writeEMBLFile( $fileName, $headerStr );

  Write out this object as an EMBL file;

=cut

##-------------------------------------------------------------------------##
sub writeEMBLFile {
  my $this      = shift;
  my $fileName  = shift;
  my $headerStr = shift;

  open OUT, ">$fileName"
      or die "$CLASS"
      . "::writeEMBLFile() Unable to open "
      . "file $fileName for writing: $!\n";

  if ( $headerStr && $headerStr ne "" ) {
    print OUT "$headerStr";
    if ( $headerStr !~ /.*\n$/ ) {
      print OUT "\n";
    }
    print OUT "XX\n";
  }

  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    print OUT "" . &toEMBLString( $this->get( $i ) );
  }
  close OUT;

}

##-------------------------------------------------------------------------##

=head2 Use: my $string = toEMBLString( $recordRef );

  Convert a single EMBL record back to EMBL format.
  NOTE: This is a static method.

=cut

##-------------------------------------------------------------------------##
sub toEMBLString {
  my $recRef = shift;    # EMBL Data Structure

  my $outStr   = toEMBLHeader( $recRef );
  my $seq      = $recRef->getSequence();
  my $seqIndex = 0;
  while ( ( $seq =~ s/^(.{60})// ) ) {
    $seqIndex += 60;
    my $lineSeq = $1;
    $lineSeq =~ s/(\w{10})/$1 /g;
    $outStr .= "     " . $lineSeq . "  $seqIndex\n";
  }
  $seqIndex += length( $seq );
  $seq =~ s/(\w{10})/$1 /g;
  my $numPadding = 68 - length( $seq );
  $outStr .= "     " . $seq . ( " " x $numPadding ) . $seqIndex . "\n";
  $outStr .= "//\n";

  return $outStr;

}

sub toEMBLHeader {
  my $recRef = shift;    # EMBL Data Structure

  my $tmpStr = "";
  my $outStr = "";

  #   ID   DF0000001; SV 4; linear; DNA; STD; UNC; 262 BP.
  $outStr = "ID   "
      . $recRef->getId() . "; " . "SV "
      . $recRef->getSeqVersion() . "; "
      . $recRef->getTopology() . "; "
      . $recRef->getMolecule() . "; "
      . $recRef->getDataclass() . "; "
      . $recRef->getDivisionsElement( 0 ) . "; "
      . $recRef->getLength()
      . " BP.\n";
  $outStr .= "XX\n";

  my $tmpVal = $recRef->getName();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $outStr .= "NM   $tmpVal\n";
    $outStr .= "XX\n";
  }

  if ( $recRef->getAccessionsCount() > 0 ) {
    my $accessionStr = join( "; ", $recRef->getAccessionsArray() );
    ## TODO Need a routine to break lines here
    $accessionStr =~ s/^/AC   /g;
    $outStr .= $accessionStr . "\n";
    $outStr .= "XX\n";
  }

  #
  my $dateOut = 0;
  $tmpVal = $recRef->getDateCreated();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $outStr .= "DT   $tmpVal\n";
    $dateOut = 1;
  }
  $tmpVal = $recRef->getDateUpdated();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $outStr .= "DT   $tmpVal\n";
    $dateOut = 1;
  }
  $outStr .= "XX\n" if ( $dateOut );

  #
  $tmpVal = $recRef->getDescription();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $tmpVal =~ s/^/DE   /gm;
    $outStr .= $tmpVal;
    $outStr .= "XX\n";
  }

  #
  # Keywords
  #
  if ( defined $recRef->getKeywordsArray() ) {
    $tmpStr = join( "; ", $recRef->getKeywordsArray() );
    $outStr .= "KW   $tmpStr.\n";
    $outStr .= "XX\n";
  }

  #
  # Organism Name
  #
  $tmpVal = $recRef->getSpeciesName();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $outStr .= "OS   $tmpVal";
    $tmpVal = $recRef->getCommonSpeciesName();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $outStr .= "($tmpVal)";
    }
    $outStr .= "\n";
  }

  #
  # Organism Classification
  #
  if ( $recRef->getClassificationCount() > 0 ) {
    $outStr .= "OC   ";
    $tmpStr = join( "; ", $recRef->getClassificationArray() ) . ".";

    # TODO Break up lines!

    $outStr .= $tmpStr . "\n";
    $outStr .= "XX\n";
  }

  #
  # References
  #
  foreach my $reference ( $recRef->getReferencesArray() ) {
    $outStr .= "RN   " . $reference->getNumber();
    $tmpVal = $reference->getPosition();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RP   /gm;
      $outStr .= $tmpVal;
    }
    $tmpVal = $reference->getComment();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RC   /gm;
      $outStr .= $tmpVal;
    }
    $tmpVal = $reference->getAuthors();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RA   /gm;
      $outStr .= $tmpVal;
    }
    $tmpVal = $reference->getTitle();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RT   /gm;
      $outStr .= $tmpVal;
    }
    $tmpVal = $reference->getLocation();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RL   /gm;
      $outStr .= $tmpVal;
    }
    $tmpVal = $reference->getXref();
    if ( defined $tmpVal && $tmpVal ne "" ) {
      $tmpVal =~ s/^/RX   /gm;
      $outStr .= $tmpVal;
    }
    $outStr .= "XX\n";
  }

  #
  # Database references
  #
  $tmpVal = $recRef->getXref();
  if ( defined $tmpVal && $tmpVal ne "" ) {
    $tmpVal =~ s/^/DR   /gm;
    $outStr .= $tmpVal;
    $outStr .= "XX\n";
  }

  #
  # Comments
  #
  foreach my $comment ( $recRef->getCommentsArray() ) {
    $comment =~ s/^/CC   /gm;
    $outStr .= $comment;
    $outStr .= "XX\n";
  }

  #
  # FT Lines
  #
  #  TODO: Some work needs to be done here to format this
  #        info with FT/FH lines.
  #my $ftLines = $recRef->getFTLines();
  #if ( $ftLines ) {
  #  $outStr .= $ftLines;
  #}

  #
  # Sequence
  #
  $outStr .= "SQ   Sequence " . $recRef->getLength() . " BP;";
  foreach my $base ( "A", "C", "G", "T", "other" ) {
    if ( defined $recRef->getCompositionElement( $base ) ) {
      $outStr .= " " . $recRef->getCompositionElement( $base ) . " $base;";
    }
  }
  $outStr .= "\n";

  return ( $outStr );

}

sub getRecordCount {
  my $this = shift;

  return ( $this->size() );

}

sub getRecord {
  my $this  = shift;
  my $index = shift;

  return $this->get( $index );
}

sub getRelease {
  my $this = shift;

  return $this->{'release'};
}

##-------------------------------------------------------------------------##
## Private Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _parseFromFile( $obj, $filename );
##
##   Parse a FASTA output file into this datastructure.
##
##-------------------------------------------------------------------------##
sub _parseFromFile {
  my $this     = shift;
  my $filename = shift;

  open IN, "<$filename"
      or die "$CLASS"
      . "::_parseFromFile() Unable to open "
      . "RepeatMasker/EMBL file $filename: $!\n";

  my $prevTag          = "";
  my $tag              = "";
  my $data             = "";
  my $recordRef        = RepeatRecord->new();
  my $firstRecordFound = 0;
  while ( <IN> ) {

    ## Ignore spacers ( and RepBase bug where they placed ;'s before tags )
    next if ( /^;?XX/ );

    ## Get file release if it exists
    if ( /^CC\s+RELEASE\s+(\d+)/ ) {
      $this->{'release'} = $1;
    }

    ## Process the end record tag
    if ( /^\/\/\s*$/ || eof( IN ) ) {
      if ( $data ne "" ) {

        # Parse last tag
        _parseTagData( $tag, $data, $recordRef );

        # Save record!
        $this->add( $recordRef );
        $recordRef = RepeatRecord->new();
      }
      $data = "";
      next;
    }

    ## Capture the tag!
    if ( /^(;?)(\S\S)(?:\s\s\s(.*))?/ ) {
      if ( $1 eq ";" ) {

        # Another Repbase error.  Special case to trap it.
        warn
"WARNING: Found error in EMBL file.  Invalid tag in line $_.  Recovering...\n";
      }
      $prevTag = $tag;
      $tag     = $2;
      my $tmpData = $3;
      ## Check validity of the tag
      if ( defined $validTagBits{$tag} ) {

        # If this tag group has ended check to see
        # if we have data to parse
        if ( $tag ne $prevTag ) {
          if ( $data ne "" ) {

            # Parse last tag
            _parseTagData( $prevTag, $data, $recordRef );
          }

          # New tag so overwrite old data
          if ( defined $tmpData ) {
            $data = $tmpData . "\n";
          }
          else {
            $data = "\n";
          }
        }
        else {

          # This is a tag continuation so
          # append the data
          if ( defined $tmpData ) {
            $data .= $tmpData . "\n";
          }
          else {
            $data .= "\n";
          }
        }
      }
      else {
        croak "Error - invalid tag $tag on line:\n$_\n   In file: $filename";
      }
    }
    else {
      if ( /^\s+[\w\s]+\d+$/ ) {
        $data .= $_;
      }
      else {
        print STDERR "Warning...unknown stuff <$_>\n";
      }
    }
  }
  close IN;
}

##-------------------------------------------------------------------------##
## Use: my _parseTagData( $tag, $data, $recRef );
##
##  Parse a RepeatMasker/Repbase EMBL data element.  This is a helper function
##  for the general parser.  This function expects data for
##  multiline tags to already be compressed into a single
##  blob.  Ie.
##
##     CC  line1
##     CC  line2
##
##  Should be passed to this function as:
##
##     $tag = "CC"
##     $data = "line1\nline2"
##
##-------------------------------------------------------------------------##
sub _parseTagData {
  my $tag    = shift;
  my $data   = shift;
  my $recRef = shift;
  return if ( $tag ne "ID" && !defined $recRef->getId() );
  if ( $tag eq "ID" ) {

# IDentification line
#
# From ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
#
#   3.4.1  The ID Line
#   The ID (IDentification) line is always the first line of an entry. The
#   format of the ID line is:
#   ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
#   The tokens represent:
#      1. Primary accession number
#      2. Sequence version number
#      3. Topology: 'circular' or 'linear'
#      4. Molecule type (see note 1 below)
#      5. Data class (see section 3.1)
#      6. Taxonomic division (see section 3.2)
#      7. Sequence length (see note 2 below)
#
#   Note 1 - Molecule type: this represents the type of molecule as stored and can
#   be any value from the list of current values for the mandatory mol_type source
#   qualifier. This item should be the same as the value in the mol_type
#   qualifier(s) in a given entry.
#   Note 2 - Sequence length: The last item on the ID line is the length of the
#   sequence (the total number of bases in the sequence). This number includes
#   base positions reported as present but undetermined (coded as "N").
#   An example of a complete identification line is shown below:
#   ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.
#
# Caveats:
#   - RepBase breaks the EMBL standard by adding "repbase" after the ID.
#   - RepBase does not provide the Sequence Version (#2), the
#     Topology (#3), or the data class (#5) tokens.
#   - Dfam/RepBase use "DNA" as the molecule type.
#   - Dfam uses "STD" it's data class.
#   - Dfam uses "UNC" as it's taxonomic division and relies on the OS/OC
#     tags to give the correct level of detail.
#   - Dfam includes a header in the file before the first record.
#     This is non-spec but we do use the valid "CC" tag to prefix each of these lines.

    # Dfam.embl Example:
    #   ID   DF0000001; SV 4; linear; DNA; STD; UNC; 262 BP.
    if ( $data =~
/(D[FR]\d+)\s*\;\s+SV\s+(\d+)\s*\;\s+linear\s*\;\s+DNA\s*\;\s+STD\s*\;\s+UNC\s*\;\s+(\d+)\s+BP.*/
        )
    {
      $recRef->setId( $1 );
      $recRef->setSeqVersion( $2 );
      $recRef->setTopology( "linear" );
      $recRef->setMolecule( "DNA" );
      $recRef->setDataclass( "STD" );
      $recRef->pushDivisions( "UNC" );
    }

    # Generic Reformatted example
    #   ID   ACROBAT1; SV 0; linear; DNA; STD; ???; 262 BP.
    elsif ( $data =~
/(\S+)\s*\;\s+SV\s+(\d+)\s*\;\s+(\S+)\s*\;\s+(\S+)\s*\;\s+(\S+)\s*\;\s+(\S+)\s*\;\s+(\d+)\s+BP.*/
        )
    {
      $recRef->setId( $1 );
      $recRef->setSeqVersion( $2 );
      $recRef->setTopology( $3 );
      $recRef->setMolecule( $4 );
      $recRef->setDataclass( $5 );
      $recRef->pushDivisions( $6 );
    }

    #
    # Legacy DfamConsensus.embl example:
    #   ID   IS1     repeatmasker; DNA;  ???;  768 BP.
    # RepBase Example
    #   ID   MER57F      repbase;    DNA;    HUM; 435 BP.
    elsif ( $data =~ /(\S+)\s+\S+\;\s+(\S+)\s*\;\s+(\S+)\s*\;\s+(\d+)\s+BP.*/ )
    {
      $recRef->setId( $1 );
      $recRef->setSeqVersion( 0 );
      $recRef->setTopology( "linear" );
      $recRef->setMolecule( $2 );
      $recRef->setDataclass( "STD" );
      $recRef->pushDivisions( $3 );
    }

    # RepBase Error:
    #   foo bar; DNA; ???; HUM BP.
    # Deprecated - they fixed the issue
    #elsif ( $data =~ /(\S+)\s+(\S+)\;\s+(\S+)\s*\;\s+(\S+)\s*\;\s+HUM\s+BP.*/ )
    #{
    #  warn "parseTagData: ID field not to EMBL"
    #      . " spec \"$_\"...recoverable\".\n";
    #  $recRef->setId( $1 );
    #  $recRef->setDataclass( $2 );
    #  $recRef->setMolecule( $3 );
    #  $recRef->pushDivisions( $4 );
    #}
    else {
      croak "parseTagData: ID field not to EMBL" . " spec \"$data\" from $_\n";
    }
  }
  elsif ( $tag eq "CC" ) {
    if ( $data =~ /RepeatMasker\s+Annotations/i ) {
      if ( $data =~ /Type:[ \t]*(\S+)[ \t]*[\n\r]/i ) {
        $recRef->setRMType( $1 );
      }
      if ( $data =~ /SubType:[ \t]*(\S+)[ \t]*[\n\r]/i ) {
        $recRef->setRMSubType( $1 );
      }

      # TODO: Decide if we want to be consistent with
      #       other lists...ie. place a coma between entities
      if ( $data =~ /Species:[ \t]*([\S \t]+)[ \t]*[\n\r]/i ) {
        my $speciesStr = $1;

        $speciesStr =~ s/\s//g;
        $speciesStr =~ s/(.*),$/$1/g;
        my @species = split( /\,/, $speciesStr );
        foreach my $species ( @species ) {
          if ( $species =~ /Search/ ) { die "Hmmm $speciesStr\n"; }
          $species =~ s/\s//g;
          $recRef->pushRMSpecies( $species );
        }
      }
      if ( $data =~ /Refineable/ ) {
        $recRef->setRMRefineable( 1 );
      }
      if ( $data =~ /SearchStages:\s*([\d\s*\,\s*]+)[\n\r]/i ) {
        my $stageStr = $1;
        $stageStr =~ s/\s//g;
        $stageStr =~ s/(.*),$/$1/g;
        my @stages = split( /\,/, $stageStr );
        foreach my $stage ( @stages ) {
          $recRef->pushRMSearchStages( $stage );
        }
      }
      if ( $data =~ /BufferStages:([^\n\r]*)[\n\r]+/i ) {
        my $bufData = $1;
        while ( $bufData =~ s/\s*(\d+\s*(?:\[\s*\d+\s*-\s*\d+\s*\]\s*)*),?//m )
        {
          my $bufStr = $1;
          $bufStr =~ s/\s//g;
          if ( $bufStr ne "" ) {
            $recRef->pushRMBufferStages( $bufStr );
          }
        }
      }
      if ( $data =~ /Description:\s+(\S.*)/i ) {
        $recRef->setRMDescription( $1 );
      }
    }

    # Coment line
    # Multi-record/multi-line comment allowed
    $recRef->pushComments( $data );
  }
  elsif ( $tag eq "AC" ) {

    # ACession line
    my @accessions = split( /\s*[;\n\r]\s*/, $data );
    foreach my $accession ( @accessions ) {
      $recRef->pushAccessions( $accession );
    }
  }
  elsif ( $tag eq "DT" ) {

    # DaTe line
    # Parse lines
    while ( $data =~ /^(.*)$/gm ) {
      $_ = $1;
      if ( /(\d+)[-\.](\S+)(-\d+)?\s+\(Rel.\s+(\d+(\.\d+)*), Created\)/ ) {
        $recRef->setDateCreated( $_ );
      }
      elsif (
/(\d+)[-\.](\w+)(-\d+)?\s+\(Rel. (\d+(\.\d+)*), Last updated, Version (\-?\d+)\)/
          )
      {
        if ( /Version\s+(\-?\d+)/ ) {
          $recRef->setSeqVersion( abs( $1 ) );
        }
        if ( /Version \-1/ ) {
          s/Version \-1/Version 1/g;
        }
        $recRef->setDateUpdated( $_ );
      }
      elsif ( /(\d+)[-\.](\w+)-(\d+)\s+\(Created\)/ ) {

        # Bug in RepBase "05-MAY-2015 (Created)"
        # Try to guess.
        #   Version major: 2015 = 20 ( so year - 1995 )
        #   Version minor: month
        my %months = (
                       "JAN" => 1,
                       "FEB" => 2,
                       "MAR" => 3,
                       "APR" => 4,
                       "MAY" => 5,
                       "JUN" => 6,
                       "JUL" => 7,
                       "AUG" => 8,
                       "SEP" => 9,
                       "OCT" => 10,
                       "NOV" => 11,
                       "DEC" => 12
        );
        my $majorVersion = $3 - 1995;
        my $minorVersion = $months{$2};
        warn "parseTagData: DT field not to EMBL"
            . " spec \"$_\"...fixing by adding \"Rel. $majorVersion.$minorVersion\".\n";
        s/Created/Rel\. $majorVersion\.$minorVersion, Created/;

        $recRef->setDateUpdated( $_ );
      }
      elsif ( /^DT$/ ) {

        # Empty DT tags...another RepBase bug
        warn "WARNING: Missing DT tag data.  Known RepBase bug...continuing\n";
      }
      else {
        croak "parseTagData: DT field not to EMBL" . " spec \"$data\"\n";
      }
    }
  }
  elsif ( $tag eq "DE" ) {

    # Description line
    $data =~ s/[\n\r]//g;
    $recRef->setDescription( $data );
  }
  elsif ( $tag eq "KW" ) {

    # KeyWord line
    my @keywords = split( /\s*[;\n\r]\s*/, $data );

    # Strip last dot ( used as punctuation in RB EMBL format )
    if ( $keywords[ $#keywords ] =~ /(\S+)\./ ) {
      $keywords[ $#keywords ] = $1;
    }
    foreach my $keyword ( @keywords ) {
      $recRef->pushKeywords( $keyword );
    }
  }
  elsif ( $tag eq "OS" ) {
    ## It looks like GIRI has been using a "CC" line following
    ## the "OS" line to record additional information about the
    ## organism name.  Almost looks like they are using it to store
    ## the common name.  There is an EMBL standard for using the "OC"
    ## line and a pair of parentheses to store the common name.
    $data =~ s/[\n\r]//g;
    if ( ( $data =~ s/\(\s*(.*)\s*\)// ) ) {
      $recRef->setCommonSpeciesName( $1 );
    }
    $recRef->setSpeciesName( $data );
  }
  elsif ( $tag eq "OC" ) {
    my @classifications = split( /\s*[;\n\r]\s*/, $data );

    # Strip last dot ( used as punctuation in RB EMBL format )
    if ( $classifications[ $#classifications ] =~ /(\S+)\./ ) {
      $classifications[ $#classifications ] = $1;
    }
    foreach my $classification ( @classifications ) {
      $recRef->pushClassification( $classification );
    }
  }
  elsif ( $tag eq "RN" ) {
    my $newRef = PubRef->new();
    $newRef->setNumber( $data );
    $recRef->pushReferences( $newRef );
  }
  elsif ( $tag eq "RP" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setPosition( $data );
  }
  elsif ( $tag eq "RC" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setComment( $data );
  }
  elsif ( $tag eq "RA" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setAuthors( $data );
  }
  elsif ( $tag eq "RT" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setTitle( $data );
  }
  elsif ( $tag eq "RL" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setLocation( $data );
  }
  elsif ( $tag eq "RX" ) {
    my $lastRef =
        $recRef->getReferencesElement( $recRef->getReferencesCount() - 1 );
    $lastRef->setXref( $data );
  }
  elsif ( $tag eq "DR" ) {
    $recRef->setXref( $data );
  }

  # SQ   Sequence 2781 BP; 892 A; 479 C; 522 G; 885 T; 3 other;
  # SQ   Sequence 262 BP; 70 A; 54 C; 62 G; 76 T; 0 other;
  elsif ( $tag eq "SQ" ) {
    my $seq = "";
    while ( $data =~ /^(.*)$/gm ) {
      $_ = $1;
      if (
/Sequence\s+(\d+)\s+BP;(?:\s+(\d+)\s+A;\s+(\d+)\s+C;\s+(\d+)\s+G;\s+(\d+)\s+T;\s+(\d+)\s+other;?)?/
          )
      {
        $recRef->setLength( $1 );
        $recRef->setComposition( 'A',     $2 );
        $recRef->setComposition( 'C',     $3 );
        $recRef->setComposition( 'G',     $4 );
        $recRef->setComposition( 'T',     $5 );
        $recRef->setComposition( 'other', $6 );
      }
      elsif ( /^\s+([acgtbdhvrykmswnxACGTBDHVRYKMSWNX\s]+)\d+$/ ) {
        my $seqLine .= $1;
        $seqLine =~ s/\s//g;
        $seq .= $seqLine;
      }
      elsif ( /^\s+con([acgtbdhvrykmswnxACGTBDHVRYKMSWNX\s]+)\d+$/ ) {

        # Special case: Fix a bug in GIRIs pipeline.  They
        # sometimes leave the prefix "con" at the start of their
        # consensus sequences.  Remove, warn, and move on.
        warn "SQ record for "
            . $recRef->getId()
            . " contains 'con' non-sequence prefix!\n";
        my $seqLine .= $1;
        $seqLine =~ s/\s//g;
        $seq .= $seqLine;
      }
      elsif ( /^\s+([acgtbdhvrykmswnxACGTBDHVRYKMSWNX\s\*lqu]+)\d+$/ ) {

        # Special case: Fix a bug in GIRIs pipeline.  They
        # sometimes have placed "*","l","q" or "u" in the sequence.  Remove,
        # warn, and move on.
        warn "SQ record for "
            . $recRef->getId()
            . " contains *|l|q|u. All known cases of extra characters in RepBase records. Removing offending characters and moving on.  Beware.\n";
        my $seqLine .= $1;
        $seqLine =~ s/[\s\*lqu]//g;
        $seq .= $seqLine;

       # TODO: Probably not sane to just assume these were insertions.  It could
       # be that the overwrote some of the consensus.  Is it still the right
       # length?
      }
      else {
        croak "SQ record for " . $recRef->getId() . " is not valid! $_\n";
      }
    }
    $recRef->setSequence( $seq );
  }
  elsif ( $tag eq "FT" ) {
    my $FTLines = $recRef->getFTLines();
    $FTLines = "" if ( !defined $FTLines );
    $recRef->setFTLines( $FTLines . $data );
  }
  elsif ( $tag eq "RG" ) {

    # Do nothing for now
  }
  elsif ( $tag eq "FH" ) {

    # Do nothing for now
  }
  elsif ( $tag eq "NM" ) {
    if ( $data =~ /(\S+)/ ) {
      $recRef->setName( $1 );
    }
  }
  else {
    croak "Not sure how to deal with $tag\n";
  }

}

1;
