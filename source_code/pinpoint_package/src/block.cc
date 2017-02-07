#include "block.h"
#include <cstdio>
#include <sstream>

block::block(int p, double pr)
{
	pos = p;
	prob = pr;
	match = -1;
}

bool block::operator<(const block &b) const
{
	if(prob > b.prob) return true;
	else return false;
}

int block::print() const
{
	printf("position = %d, prob = %.6lf, match = %d\n", pos, prob, match);
	return 0;
}

int block::distance(const block &b) const
{
	if(pos <= b.pos) return b.pos - pos;
	else return pos - b.pos;
}
