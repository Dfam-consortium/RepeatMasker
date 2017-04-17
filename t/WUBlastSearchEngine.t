use strict;
use WUBlastSearchEngine;
use SearchResult;
use SearchResultCollection;
use SearchEngineI;
use RepeatMaskerConfig;
use Data::Dumper;
use Test::More;

my $wrkDir = "wrkdir";
if (! -d $wrkDir ) { 
  system( "mkdir $wrkDir" );
}


#------------------------------------------------------------------------##
# Tests
#------------------------------------------------------------------------##


##
## Object Creation
##
my $wuse = WUBlastSearchEngine->new(
                             pathToEngine => $RepeatMaskerConfig::WUBLASTP_PRGM,
                             DEBUG        => 0 );

ok( $wuse, "Testing object creation ->new( pathToEngine => $RepeatMaskerConfig::WUBLASTP_PRGM )" );
  

##
## Build a database
##
unlink( "libraries/alu.fa.ahd" ) if ( -e "libraries/alu.fa.ahd" );
unlink( "libraries/alu.fa.atb" ) if ( -e "libraries/alu.fa.atb" );
unlink( "libraries/alu.fa.bsq" ) if ( -e "libraries/alu.fa.bsq" );
system(   "$RepeatMaskerConfig::SETDB_PRGM libraries/alu.fa > "
        . "$wrkDir/setdb.log 2>&1" );

ok( ( -s "libraries/alu.fa.ahd" &&
      -s "libraries/alu.fa.atb" &&
      -s "libraries/alu.fa.bsq" ), 
    "Testing wublast setdb command on libraries/alu.fa");

##
## Search  
##
my ( $resultCode, $rc ) = $wuse->search( matrix=>"matrices/aa/14p35g.matrix",
                                 scoreMode=>SearchEngineI::complexityAdjustedScoreMode,
                                 gapInit => -30,
                                 insGapExt => -6,
                                 delGapExt => -6,
                                 minScore => 200,
                                 maskLevel=>101,
                                 query=>"seqs/general/alu-ex.fa",
                                 subject=>"libraries/alu.fa" );

ok( $rc->size() > 1, "Testing basic run of WUBlastSearchEngine.pm" );


##
## Masklevel
##
##./MaskerAid -masklevel 10 -minscore 200 -gap_init -30 -ins_gap_ext -6 -del_gap_ext -6 -matrix 14p35g.matrix ../seqs/general/NT_001493.seq ../libraries/alu.fa  > ../seqs/general/NT_001493-wucmp-1/run1.out
##
##  1985 14.04 0.34 0 ref|NT_001493|HsX_1591 13509 13800 (18871) AluSz#SINE/Alu 1 293 (19)
##  2051 13 2.33 0 ref|NT_001493|HsX_1591 14385 14684 (17987) AluSx1#SINE/Alu 1 307 (5)
##  574 18.44 0 1.98 ref|NT_001493|HsX_1591 27345 27447 (5224) AluJr#SINE/Alu 1 101 (211)
##  472 29.29 5.73 0.6 ref|NT_001493|HsX_1591 30054 30210 (2461) AluYGibG6#SINE/Alu 117 281 (30)
##
( $resultCode, $rc ) = $wuse->search( matrix=>"/home/rhubley/projects/RepeatMasker/Matrices/wublast/aa/14p35g.matrix",
                                 scoreMode=>SearchEngineI::complexityAdjustedScoreMode,
                                 gapInit => -30,
                                 insGapExt => -6,
                                 delGapExt => -6,
                                 minScore => 200,
                                 maskLevel=> 10,
                                 query=>"seqs/general/NT_001493.seq",
                                 subject=>"libraries/alu.fa" );

ok( $rc->size() == 4, "Object: masklevel output" ) or
  diag(" search: " . $wuse->getParameters() . " size = " . $rc->size() );
ok( $rc->get(0)->getScore() == 1985 &&
    $rc->get(1)->getScore() == 2051 &&
    $rc->get(2)->getScore() == 574 &&
    $rc->get(3)->getScore() == 472, "Object: Complexity adjusted scores" ) or
  diag(" search: " . $wuse->getParameters() . " size = " . $rc->size() );

#for ( my $i = 0; $i < $rc->size(); $i++ )
#{
#  my $r = $rc->get( $i );
#  print "" . $r->toStringFormatted( SearchResult::AlignWithSubjSeq ) . "";
#  print "" . $r->toStringFormatted( ) . "";
#}

done_testing();

1;
