#ifndef __POSITION_H__
#define __POSITION_H__

#include <vector>
#include <string>

using namespace std;

class position
{
public:
	position(const string &s);

public:
	// load from prediction file
	int tlab;		// true label
	double xyz[3];	// prob for 0, 1, and 2
	double pp;		// largest prob
	int pred;		// predicted label
	double pabd;	// predicted abd

	// load from sample file
	double tabd;	// true abundance
	double rabd;	// read abundance

public:
	int print();
};

#endif
