#include "pinpoint.h"
#include "config.h"

#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc == 1) 
	{
		printf("usage: %s: <sample-file> <prediction-file> <output-file> [-p probability-threshold] [-w window-size]\n", argv[0]);
		return 0;
	}

	parse_arguments(argc, argv);

	pinpoint pp(argv[1], argv[2], argv[3]);
	pp.solve();

	return 0;
}
