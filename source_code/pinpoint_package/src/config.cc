#include "config.h"
#include <string>
#include <cstdlib>

int block_size = 10;
double min_prob = 0.6;
int block_size1 = 10;
double min_prob1 = 0.90;
int max_correct_distance = 100;
double min_accept_expression = 5.0;
int min_block_distance = 10000;
bool extended = false;

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-w")
		{
			block_size = 2 * atoi(argv[i + 1]);
			i++;
		}
		if(string(argv[i]) == "-p")
		{
			min_prob = atof(argv[i + 1]);
			i++;
		}
	}
	return 0;
}
