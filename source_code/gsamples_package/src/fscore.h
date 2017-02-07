#ifndef __FSCORE_H__
#define __FSCORE_H__

#include <vector>
#include <fstream>

using namespace std;

class fscore
{
public:
	int w;
	vector<int> sa;
	vector<int> sb;
	vector<double> va;
	vector<double> vb;

public:
	int init(int k);
	int write(ofstream &fout);
};

#endif
