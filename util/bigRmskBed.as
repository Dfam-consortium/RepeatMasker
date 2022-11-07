table bigRmskBed
"Repetitive Element Annotation"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of visualization on chromosome"
uint    chromEnd;	"End position of visualation on chromosome"
string  name;		"Name repeat, including the type/subtype suffix"
uint    score;		"Divergence score"
char[1] strand;		"+ or - for strand"
uint    thickStart;	"Start position of aligned sequence on chromosome"
uint    thickEnd;	"End position of aligned sequence on chromosome"
uint  	reserved;	"Reserved"
uint    blockCount;	"Count of sequence blocks"
lstring blockSizes;     "A comma-separated list of the block sizes(+/-)"
lstring blockStarts;    "A comma-separated list of the block starts(+/-)"
uint    id;             "A unique identifier for the joined annotations in this record"
lstring description;    "A comma separated list of technical annotation descriptions"
)
