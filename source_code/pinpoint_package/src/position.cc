#include "position.h"
#include <cstdio>
#include <sstream>
#include <cassert>
#include <cmath>

position::position(const string &s)
{
	char line[1024000];
	stringstream sstr(s);
	sstr >> tlab >> line;
	sstr >> xyz[0] >> xyz[1] >> xyz[2] >> line;
	sstr >> pp >> pred >> pabd;
	pabd = exp(pabd);
}

int position::print()
{
	printf("%d -> %.2lf %.2lf %.2lf -> %.2lf %d, tabd = %.3lf rabd = %.3lf [pabd = %.3lf]\n", 
			tlab, xyz[0], xyz[1], xyz[2], pp, pred, tabd, rabd, pabd);
	return 0;
}
