#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include "position.h"
#include "block.h"
#include <fstream>
#include <vector>

using namespace std;

class sample
{
public:
	string header;
	vector<position> positions;

	vector<block> blocks0;
	vector<block> blocks2;
	vector<int> blocks1;

	bool accept;
	int correct0;
	int correct1;
	int correct2;
	int label0;
	int label1;
	int label2;
	int predict0;
	int predict1;
	int predict2;

	vector<int> splist;
	vector<double> vabd;
	double abdratio;

public:
	int add_position(const string &s);
	int clear();
	int process();
	int print(int index);
	int print_result(int index);
	int write(ofstream &fout);

private:
	int qualify();
	int build_blocks(int ff, vector<block> &blocks);
	int align_blocks(int ff, vector<block> &blocks, int &ncorrect, int &nlabel);
	int build_blocks1();
	int build_splice_positions();
	int build_abundance();
	int assess_abundance();
};

#endif
