#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RepeatMaskerConfig.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##      Arian Smit <asmit@systemsbiology.org>
##  Description:
##      This is the main configuration file for the RepeatMasker
##      program suite.  Before you can run the programs included
##      in this package you will need to edit this file either
##      through the use of the configure script ( preferred ) or
##      by hand-editing this file.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2005-2019 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
package RepeatMaskerConfig;
use FindBin;
use File::Basename;
use Data::Dumper;
require Exporter;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
@ISA         = qw(Exporter);
my $CLASS = "RepeatMaskerConfig";

BEGIN {
##----------------------------------------------------------------------##
##     CONFIGURE THE FOLLOWING PARAMETERS FOR YOUR INSTALLATION         ##
##                                                                      ##
##  This file may be hand edited or preferably, modified with the       ##
##  the use of the "configure" script.                                  ##
##                                                                      ##
##  In the following section default values for paths/programs          ##
##  may be hardcoded. Each parameter appears in a block of text         ##
##  below as:                                                           ##
##                                                                      ##
##    PARAMETER_NAME => {                                               ##
##                 ...                                                  ##
##              value => "/usr/local/bin/foo"                           ##
##                      }                                               ##
##                                                                      ##
##  only change the "value" field for each parameter. If you are        ##
##  unsure of how to edit this file, simply use the the ./configure     ##
##  script provided with the package to perform the same task.  For     ##
##  more details on configuring this package please see the README      ##
##  file.                                                               ##
##                                                                      ##
  ## STCFG --do-not-remove--
  $configuration = {
          'ABBLAST_DIR' => {
                             'command_line_override' => 'abblast_dir',
                             'description' => 'The path to the installation of the ABBLAST sequence alignment program.',
                             'environment_override' => 'ABBLAST_DIR',
                             'expected_binaries' => [
                                                      'xdformat',
                                                      'blastp'
                                                    ],
                             'expected_files' => [],
                             'param_type' => 'directory',
                             'required' => 0,
                             'value' => ''
                           },
          'CROSSMATCH_DIR' => {
                                'command_line_override' => 'crossmatch_dir',
                                'description' => 'The path Phil Green\'s cross_match program ( phrap program suite ).',
                                'environment_override' => 'CROSSMATCH_DIR',
                                'expected_binaries' => [
                                                         'cross_match'
                                                       ],
                                'expected_files' => [],
                                'param_type' => 'directory',
                                'required' => 0,
                                'value' => ''
                              },
          'DEFAULT_SEARCH_ENGINE' => {
                                       'command_line_override' => 'default_search_engine',
                                       'description' => 'The default search engine to use',
                                       'param_type' => 'value',
                                       'required' => 1,
                                       'value' => ''
                                     },
          'HMMER_DIR' => {
                           'command_line_override' => 'hmmer_dir',
                           'description' => 'The path to the HMMER profile HMM search software.',
                           'environment_override' => 'HMMER_DIR',
                           'expected_binaries' => [
                                                    'nhmmer'
                                                  ],
                           'expected_files' => [],
                           'param_type' => 'directory',
                           'required' => 0,
                           'value' => ''
                         },
          'LIBDIR' => {
                        'command_line_override' => 'libdir',
                        'description' => 'Path to the RepeatMasker libraries directory.',
                        'environment_override' => 'LIBDIR',
                        'expected_binaries' => [],
                        'expected_files' => [
                                              'RepeatAnnotationData.pm',
                                              'RepeatPeps.lib',
                                              'RMRBMeta.embl'
                                            ],
                        'param_type' => 'directory',
                        'required' => 0,
                        'value' => ''
                      },
          'RMBLAST_DIR' => {
                             'command_line_override' => 'rmblast_dir',
                             'description' => 'The path to the installation of the RMBLAST sequence alignment program.',
                             'environment_override' => 'RMBLAST_DIR',
                             'expected_binaries' => [
                                                      'rmblastn',
                                                      'dustmasker',
                                                      'makeblastdb',
                                                      'blastdbcmd',
                                                      'blastdb_aliastool',
                                                      'blastn'
                                                    ],
                             'expected_files' => [],
                             'param_type' => 'directory',
                             'required' => 0,
                             'value' => ''
                           },
          'TRF_PRGM' => {
                          'command_line_override' => 'trf_prgm',
                          'description' => 'The full path including the name for the TRF program.',
                          'environment_override' => 'TRF_PRGM',
                          'param_type' => 'program',
                          'required' => 1,
                          'value' => ''
                        }
        };

  ## EDCFG --do-not-remove--
  

##                                                                      ##
##                      END CONFIGURATION AREA                          ##
##----------------------------------------------------------------------##
##----------------------------------------------------------------------##
##  Do not edit below this line                                         ##
##----------------------------------------------------------------------##

#
# Current version of the software
#
$VERSION = "4.1.5";

#
# Set this flag to default to debug mode for the entire package
#
$DEBUGALL = 0;

# 
# Prompt for a specific parameter and update the object
#
sub promptForParam{
  my $param = shift;
  my $screenHdr = shift;

  if ( ! exists $configuration->{$param} ) {
    return;
  }

  # Grab defaults
  my $defaultValue;
  if ( $configuration->{$param}->{'param_type'} eq "directory" )
  {
    if ( exists $configuration->{$param}->{'expected_binaries'} &&
         @{$configuration->{$param}->{'expected_binaries'}} ) 
    {
      my $binary = $configuration->{$param}->{'expected_binaries'}->[0];
      $defaultValue = `/usr/bin/which $binary`;
      if ( $defaultValue !~ /\/usr\/bin\/which:/ )
      {
        $defaultValue =~ s/[\n\r\s]+//g;
        $defaultValue =~ s/^(.*)\/$binary/$1/;
      }else {
        $defaultValue = "";
      }
    }
    if ( $defaultValue eq "" && exists $configuration->{$param}->{'value'} &&
         -d $configuration->{$param}->{'value'} ) {
      $defaultValue = $configuration->{$param}->{'value'};
    }
  }elsif ( $configuration->{$param}->{'param_type'} eq "program" )
  {
    # The program type is used in cases where a single
    # script/binary is referenced and may not have the
    # exact name we expect.  TRF is a good example of this
    # as the binary is often distributed with names like:
    # trf409.linux64 etc..
    if ( exists $configuration->{$param}->{'value'} ) {
      my($binary, $dirs, $suffix) = fileparse($configuration->{$param}->{'value'});
      $defaultValue = `/usr/bin/which $binary`;
      if ( $defaultValue !~ /\/usr\/bin\/which:/ )
      {
        $defaultValue =~ s/[\n\r\s]+//g;
      }else {
        $defaultValue = $configuration->{$param}->{'value'};
      }
    }
  }

  my $value = "";
  my $validParam;
  do { 
    $validParam = 1;
    system("clear");
    if ( $screenHdr ) {
      print "$screenHdr\n";
    }else {
      print "\n\n\n\n";
    }
    print "" . $configuration->{$param}->{'description'} . "\n";

    # Prompt and get the value
    if ( $defaultValue ) {
      print "$param [$defaultValue]: ";
    }else {
      print "$param: ";
    }
    $value = <STDIN>;
    $value =~ s/[\n\r]+//g;
    if ( $value eq "" && $defaultValue )
    {
      $value = $defaultValue;
    }

    if ( $configuration->{$param}->{'param_type'} eq "directory" )
    {
      if ( -d $value ) {
        foreach my $file ( @{$configuration->{$param}->{'expected_files'}} )
        {
          if ( ! -s "$value/$file" ) 
          {
            print "\nCould not find the required file \"$file\" inside\n"
               .  "the directory \"$value\"!\n\n";
            $validParam = 0;
            last;
          }
        }
        foreach my $binary ( @{$configuration->{$param}->{'expected_binaries'}} )
        {
          if ( ! -x "$value/$binary" )
          {
            print "\nCould not find the required program \"$binary\" inside\n" 
                . "the directory \"$value\"!\n\n";
            $validParam = 0;
            last;
          }elsif ( -d "$value/$binary" )
          {
            print "\nCould not find the required program \"$binary\" inside\n"
                . "the directory \"$value\"!  It appears to be the name of a\n"
                . "subdirectory.\n\n";
            $validParam = 0;
            last;
          }
        }
      }else { 
          print "\nCould not find the \"$value\" directory.\n\n";
          $validParam = 0;
      }   
    }elsif ( $configuration->{$param}->{'param_type'} eq "program" )
    {
      if ( ! -x $value ) {
        print "\nThe program \"$value\" doesn't appear to exist\n"
            . "or it's not executable!\n\n";
        $validParam = 0;
      }elsif ( -d $value ){
        print "\nThe value \"$value\" appears to be a directory rather\n" 
            . "than an executable binary or script!\n\n";
        $validParam = 0;
      }
    }

    if ( $validParam == 0 )
    {
      print "<PRESS ENTER TO CONTINUE, CTRL-C TO BREAK>\n";
      <STDIN>;
    }
  }while ( $validParam == 0 );
  $configuration->{$param}->{'value'} = $value;
}

#
# Validate parameter
#
sub validateParam{
  my $param = shift;
  my $new_setting = shift;

  if ( ! exists $configuration->{$param} ) {
    return 0;
  }

  my $value = $configuration->{$param}->{'value'};
  $value = $new_setting if ( defined $new_setting );

  # Always assume the "good" in parameters...
  my $validParam = 1;
  if ( $configuration->{$param}->{'param_type'} eq "directory" )
  {
    if ( -d $value ) {
      foreach my $file ( @{$configuration->{$param}->{'expected_files'}} )
      {
        if ( ! -s "$value/$file" ) 
        {
          $validParam = 0;
          last;
        }
      }
      foreach my $binary ( @{$configuration->{$param}->{'expected_binaries'}} )
      {
        if ( ! -x "$value/$binary" )
        {
          $validParam = 0;
          last;
        }elsif ( -d "$value/$binary" )
        {
          $validParam = 0;
          last;
        }
      }
    }else { 
        $validParam = 0;
    }   
  }elsif ( $configuration->{$param}->{'param_type'} eq "program" )
  {
    if ( ! -x $value ) {
        $validParam = 0;
    }elsif ( -d $value ){
      $validParam = 0;
    }
  }

  return $validParam;
}
 
#
# Clear all settings
#
sub clearValues{
  foreach my $param ( keys %$configuration ) {
    if ( $configuration->{$param}->{'value'} )
    {
      $configuration->{$param}->{'value'} = '';
    }
  }
  &updateConfigFile();
}

#
# Update this file ( beware: self modifying code ) new
# paramter settings.
#
sub updateConfigFile{
  open IN,"<$CLASS.pm" or die;
  open OUT,">new-$CLASS.pm" or die;
  my $inCfg;
  $Data::Dumper::Sortkeys = 1;
  while (<IN>){
    if ( /##\s+STCFG/ ) {
      $inCfg = 1;
      print OUT;
      my $cStr = Dumper($configuration);
      $cStr =~ s/\$VAR1/  \$configuration/;
      print OUT "$cStr\n";
    }elsif ( /##\s+EDCFG/ ) {
      $inCfg = 0;
      print OUT;
    }elsif ( ! $inCfg ) {
      print OUT;
    }
  }
  close IN;
  close OUT;
  rename("$CLASS.pm", "$CLASS.pm.bak");
  rename("new-$CLASS.pm", "$CLASS.pm");
}

# 
# Create a GetOpt list for the command-line parameters defined
# in this configuration file.  These may be appended to a program's
# existing GetOpt parameter as:
#
#     push @getopt_args, RepeatMaskerConfig::getCommandLineOptions();
#
sub getCommandLineOptions {
  my @options = ();
  foreach my $param ( keys %$configuration ) {
    if ( exists $configuration->{$param}->{'command_line_override'} ) {
      push @options, "-" . $configuration->{$param}->{'command_line_override'} . "=s";
    }
  }
  return @options;
}

#
# Get POD documentation to add to the existing program POD stored in 
# the main script. 
#
sub getPOD {
  my $pod_str;
  foreach my $param ( keys %$configuration ) {
    if ( exists $configuration->{$param}->{'command_line_override'} ) {
      $pod_str .= "=item -" . $configuration->{$param}->{'command_line_override'} . " <string>\n\n";
      $pod_str .= $configuration->{$param}->{'description'} . "\n\n";
    }
  }
  if ( $pod_str ) {
    return( "\n=over 4\n\n" . $pod_str . "=back\n\n" );
  }
  return;
}

#
# After GetOpt has filled in the options hash simply pass it to
# this function to perform resolution.  The following precedence
# is used:
# 
#    1. Command Line Parameter
#    2. Environment Variable
#    3. Configuration File
#
# This will update the configuration{param}->{'value'} for use 
# in the main program.
#
sub resolveConfiguration {
  my $options = shift;
  
  foreach my $param ( keys %$configuration ) {
    if ( exists $options->{$configuration->{$param}->{'command_line_override'}} )
    {
      $configuration->{$param}->{'value'} = $options->{$configuration->{$param}->{'command_line_override'}};
    } elsif ( exists $ENV{$configuration->{$param}->{'environment_override'}} )
    {
      $configuration->{$param}->{'value'} = $ENV{$configuration->{$param}->{'environment_override'}};
    }
  }
}


}

1;
