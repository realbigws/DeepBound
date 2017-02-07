#ifndef __BLOCK_H__
#define __BLOCK_H__

#include <string>
#include <vector>
#include <set>

#include "fscore.h"
#include "interval_map.h"

using namespace std;

class block
{
public:
	string chrm;
	int32_t pos;
	static int index;
	int ltype;
	int rtype;

	vector<int> labels;	// labels
	vector<int> abd;	// real abundance
	vector<int> abl;	// abundance label

	// features
	string seq;			// sequence
	vector<int> s;		// abundance
	vector<double> q;	// ave-quality
	fscore fs20;		// window 20 scores
	fscore fs50;		// window 50 scores
	fscore fs100;		// window 100 scores

	vector<int32_t> ss;	// start-boundaries
	vector<int32_t> tt;	// end-boundaries

public:
	int clear();
	int evaluate(int a, int b, double &ave, double &dev);
	int build_labels();
	int build_labels(const set<int32_t> &ss, const set<int32_t> &tt);
	int build_abundance(const join_interval_map &jmap);
	int build_feature_score(fscore &fs);
	int build_features();

	bool qualify_boundary_training();
	bool qualify_abundance_training();
	int write_samples(ofstream &fout);
	int write_index();

	int print_labels();
};

#endif
