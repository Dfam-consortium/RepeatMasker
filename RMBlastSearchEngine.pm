#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RMBlastSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the 3.x series of rmblastn.
##
##      The 3.x series is a reimplementation of rmblastn.  Its alignment
##      semantics match the 2.x series (same scoring, same x-drop
##      parameters, same tab delimited output columns), but its command
##      line syntax does not, and it searches a FASTA or 2bit file
##      directly rather than a database built by makeblastdb.
##
##      This class therefore inherits everything from
##      NCBIBlastSearchEngine except the two things that differ: how
##      parameters are spelled on the command line, and what it means to
##      prepare a subject.  NCBIBlastSearchEngine::new() substitutes this
##      class when it probes a 3.x binary, so callers do not choose.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2026 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
#
# ChangeLog
#
#     $Log$
#
###############################################################################

=head1 NAME

RMBlastSearchEngine

=head1 SYNOPSIS

use RMBlastSearchEngine;

  my $engine = RMBlastSearchEngine->new(
                                 pathToEngine => "/usr/local/rmblast/bin/rmblastn" );

  my $subject = $engine->prepareSubject( "mylib.fa" );
  $engine->setSubject( $subject );
  $engine->setQuery( "myseq.fa" );
  my ( $status, $results ) = $engine->search();

Usually obtained indirectly: NCBIBlastSearchEngine->new() returns one of
these when the binary it probes reports major version 3 or above.

=head1 DESCRIPTION

A SearchEngineI implementation for rmblastn 3.x.

=head1 SEE ALSO

=over 4

SearchEngineI, NCBIBlastSearchEngine

=back

=head1 COPYRIGHT

Copyright 2026 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut

package RMBlastSearchEngine;
use strict;
use SearchEngineI;
use NCBIBlastSearchEngine;
use SearchResultCollection;
use Data::Dumper;
use FileHandle;
use File::Basename;
use File::Spec;
use Carp;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

# Not really a NCBI produced engine, but rather a port of the NCBI codebase,
# but this simplifies code that expected a NCBIBlastSearchEngine.
@ISA = qw(Exporter NCBIBlastSearchEngine);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $VERSION = 0.1;
my $CLASS   = "RMBlastSearchEngine";

##
## Option spellings for the 3.x command line.
##
## The 3.x series uses GNU style long options where 2.x used single-dash
## NCBI toolkit options.  The keys below are the canonical parameter names
## produced by NCBIBlastSearchEngine::_computeSearchParameters(); the
## values are how this engine spells them.
##
## Checked against rmblastn 3.0.6 --help on 2026-09-03.  Note that 3.x
## also accepts the 2.x spellings ( -word_size, -word-size, --word_size
## and --word-size are all equivalent ), so the long forms below are a
## preference rather than a requirement.  The one option that genuinely
## does not exist is num_alignments; see _renderParameters().
##
my %OPT = (
            'db'                   => '--db',
            'query'                => '--query',
            'gapopen'              => '--gapopen',
            'gapextend'            => '--gapextend',
            'mask_level'           => '--mask-level',
            'complexity_adjust'    => '--complexity-adjust',
            'word_size'            => '--word-size',
            'xdrop_ungap'          => '--xdrop-ungap',
            'xdrop_gap'            => '--xdrop-gap',
            'xdrop_gap_final'      => '--xdrop-gap-final',
            'min_raw_gapped_score' => '--min-raw-gapped-score',
            'dust'                 => '--dust',
            'outfmt'               => '--outfmt',
            'num_threads'          => '--num-threads',
            'mt_mode'              => '--mt-mode',
            'matrix'               => '--matrix',
            'gilist'               => '--gilist',
);

##
## 3.x takes the matrix as a path and has no BLASTMAT lookup, unlike the
## 2.x convention of a bare filename resolved through the environment.
## Confirmed against 3.0.6.
##
my $MATRIX_AS_PATH = 1;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # NCBIBlastSearchEngine::new() will not re-bless when it is called on
  # this class, so this simply inherits.
  return $class->SUPER::new( %nameValuePairs );
}

##-------------------------------------------------------------------------##
## Subject preparation
##-------------------------------------------------------------------------##

=head2 prepareSubject()

  Use: my $subjectPath = prepareSubject( $seqFile, ... );

  The 3.x series reads a FASTA or 2bit file directly, so there is nothing
  to build.  The outputDir and dbName parameters are accepted and ignored
  so that callers need not branch on engine version.

  The sequence file itself is the database.  A caller that previously
  wrote index files into a scratch directory and searched against that
  directory must use the path returned here, which points at the original
  sequence file.

=cut

##-------------------------------------------------------------------------##
sub prepareSubject {
  my $this    = shift;
  my $seqFile = shift;

  croak $CLASS
      . "::prepareSubject(): Sequence file ($seqFile) does not "
      . "exist or is empty!\n"
      if ( !-s $seqFile );

  return $seqFile;
}

##-------------------------------------------------------------------------##

=head2 isSubjectPrepared()

  Use: my $bool = isSubjectPrepared( $path );

  True if $path is a readable, non-empty sequence file.

=cut

##-------------------------------------------------------------------------##
sub isSubjectPrepared {
  my $this = shift;
  my $path = shift;

  return 0 if ( !defined $path );

  return ( -s $path );
}

##-------------------------------------------------------------------------##

=head2 getSubjectArtifacts()

  Preparation produces no files, so this is always empty.

=cut

##-------------------------------------------------------------------------##
sub getSubjectArtifacts {
  return ();
}

##-------------------------------------------------------------------------##

=head2 getPathToDBFormatter()

  There is no database formatting program in the 3.x distribution.

=cut

##-------------------------------------------------------------------------##
sub getPathToDBFormatter {
  my $this = shift;

  croak $CLASS
      . "::getPathToDBFormatter(): The 3.x series does not use a "
      . "separate database formatting program.\n";
}

##-------------------------------------------------------------------------##

=head2 setSubjectIDList()

  Use: my $oldValue = setSubjectIDList( $file );

  The 3.x series reads the text list directly, so the file is used as
  given.  A line holding a bare number matches a subject named "gi|N".

=cut

##-------------------------------------------------------------------------##
sub setSubjectIDList {
  my $this  = shift;
  my $value = shift;

  croak $CLASS
      . "::setSubjectIDList(): List file ($value) does not exist or "
      . "is empty!\n"
      if ( defined $value && !-s $value );

  return $this->SearchEngineI::setSubjectIDList( $value );
}

sub _subjectIDListForEngine {
  my $this = shift;

  return $this->getSubjectIDList();
}

##-------------------------------------------------------------------------##
## Command line rendering
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
##  Use: my $string = $this->_renderParameters( $paramsRef );
##
##  Render the computed parameter values using the 3.x option syntax.  The
##  values themselves come from NCBIBlastSearchEngine, so the scoring,
##  x-drop and threading decisions are shared with the 2.x engine and are
##  not restated here.
##-------------------------------------------------------------------------##
sub _renderParameters {
  my $this = shift;
  my $p    = shift;

  my $parameters = "";

  # Emitted in the same order as the 2.x renderer, so the two command
  # lines line up when comparing them during validation.
  #
  # num_alignments is deliberately absent.  2.x needed -num_alignments
  # 9999999 to lift the default reporting cap; 3.x has no such cap and no
  # such option, and passing it is a hard argument error.
  $parameters .= " $OPT{'db'} " . $p->{'db'};
  $parameters .= " $OPT{'query'} " . $p->{'query'};
  $parameters .= " $OPT{'gapopen'} " . $p->{'gapopen'};
  $parameters .= " $OPT{'gapextend'} " . $p->{'gapextend'};
  $parameters .= " $OPT{'mask_level'} " . $p->{'mask_level'}
      if ( defined $p->{'mask_level'} );
  $parameters .= " $OPT{'complexity_adjust'}" if ( $p->{'complexity_adjust'} );
  $parameters .= " $OPT{'word_size'} " . $p->{'word_size'};

  if ( defined $p->{'xdrop_ungap'} ) {
    $parameters .= " $OPT{'xdrop_ungap'} " . $p->{'xdrop_ungap'};
    $parameters .= " $OPT{'xdrop_gap_final'} " . $p->{'xdrop_gap_final'};
    $parameters .= " $OPT{'xdrop_gap'} " . $p->{'xdrop_gap'};
  }

  if ( defined $p->{'min_raw_gapped_score'} ) {
    $parameters .=
        " $OPT{'min_raw_gapped_score'} " . $p->{'min_raw_gapped_score'};
    $parameters .= " $OPT{'dust'} " . $p->{'dust'};
  }

  $parameters .=
        " $OPT{'outfmt'}=\"6 "
      . join( " ", @{ $p->{'outfmt_fields'} } ) . "\""
      if ( defined $p->{'outfmt_fields'} );

  $parameters .= " $OPT{'num_threads'} " . $p->{'num_threads'};
  $parameters .= " $OPT{'mt_mode'} " . $p->{'mt_mode'}
      if ( defined $p->{'mt_mode'} );

  if ( defined $p->{'matrix'} ) {
    my $matrix = $p->{'matrix'};
    $matrix = $p->{'matrix_dir'} . "/" . $matrix
        if ( $MATRIX_AS_PATH && defined $p->{'matrix_dir'} );
    $parameters .= " $OPT{'matrix'} " . $matrix;
  }
  $parameters .= " $OPT{'gilist'} " . $p->{'gilist'}
      if ( defined $p->{'gilist'} );

  return $parameters;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

1;
