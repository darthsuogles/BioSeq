// Ultraconservative Interspecies DNA Sequence Detection
#include <time.h>

//#include "pre_parser.h"
#include "parser.h"
#include "eval.h"

using namespace std;

#define debug 0

// change these variables as you need
const int SPECIE_THRESH = 40;
const int LEN_THRESH = 100;
const double PERCENT_THRESH = 0.8;


int main()
{
	time_t timer_start=0, timer_end=0;
	timer_start = time(NULL);
	//pre_parser(SPECIE_THRESH, LEN_THRESH, PERCENT_THRESH);
	int block_count = parser(SPECIE_THRESH, LEN_THRESH, PERCENT_THRESH);
	timer_end = time(NULL);
	printf("\n\n%d blocks are valid", block_count);
	printf("\ntime spent in total %f secs\n\n", difftime(timer_end, timer_start));
	return 0;
}

