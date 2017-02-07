#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <utility>
using namespace std;

extern int block_size;
extern int min_block_distance;
extern int max_correct_distance;
extern double min_prob;
extern double min_accept_expression;

extern int block_size1;
extern double min_prob1;

extern bool extended;

typedef pair<int, int> PI;

int parse_arguments(int argc, const char ** argv);

#endif
