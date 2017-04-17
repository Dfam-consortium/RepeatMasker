use strict;
use WUBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchResult;
use SearchResultCollection;
use SearchEngineI;
use Data::Dumper;
use Test::More;


my $wrkDir = "wrkdir";
if (! -d $wrkDir ) { 
  system( "mkdir $wrkDir" );
}


#------------------------------------------------------------------------##
# Tests
#------------------------------------------------------------------------##


my $resultCol = CrossmatchSearchEngine::parseOutput( searchOutput => "seqs/masklevel/seq1-2-ml101.out" );
ok( defined $resultCol && $resultCol->size() == 12, "Building a SearchResultColleciton from Crossmatch output" );

$resultCol->maskLevelFilter( value => 30 );
# was 97% = 2
ok( $resultCol->size() == 2 , "Testing seq1-2-ml101.out with masklevel 30" ) or
  diag("Actual collection size is " . $resultCol->size() );

$resultCol = CrossmatchSearchEngine::parseOutput( searchOutput => "seqs/masklevel/seq1-2-ml101.out" );
$resultCol->maskLevelFilter( value => 98 );
ok( $resultCol->size() == 3 , "Testing seq1-2-ml101.out with masklevel 98" ) or
  diag("Actual collection size is " . $resultCol->size() );

$resultCol = CrossmatchSearchEngine::parseOutput( searchOutput => "seqs/masklevel/hum-1-sine-ml101.out" );
ok( defined $resultCol && $resultCol->size() == 3882, "Testing hum-1-sine-ml101.out with masklevel 101" ) or
  diag("Actual collection size is " . $resultCol->size() );

$resultCol->maskLevelFilter( value => 90 );
ok( defined $resultCol && $resultCol->size() == 176, "Testing hum-1-sine-ml101.out with masklevel 90" ) or
  diag("Actual collection size is " . $resultCol->size() );

# TODO: A more comprehensive comparison with the actual output results.
#$resultCol->write( "tmp.out", SearchResult::NoAlign );
#$resultCol = CrossmatchSearchEngine::parseOutput( searchOutput => "seqs/masklevel/hum-1-sine-ml90.out" );
#$resultCol->write( "tmp2.out", SearchResult::NoAlign );
#
done_testing();

1;
