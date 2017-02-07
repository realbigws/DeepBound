#ifndef __BLOCK_H__
#define __BLOCK_H__

using namespace std;

class block
{
public:
	block(int p, double pr);
	bool operator<(const block &b) const;

public:
	int pos;
	int match;
	double prob;

public:
	int print() const;
	int distance(const block &b) const;
};

#endif
