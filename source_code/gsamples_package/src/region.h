#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "block.h"
#include "binomial.h"
#include "fasta.h"

using namespace std;

class region
{
public:
	region(const string &_chrm, int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap, const split_interval_map *_qmap, fasta &_fa);
	~region();

private:
	string chrm;					// chrm
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *imap;	// pointer to a interval map
	const split_interval_map *qmap;	// pointer to a quality map
	fasta &fa;						// fasta

	int32_t lcore;					// left core position
	int32_t rcore;					// right core position
	bool empty;						// whether this region is completely spliced

public:
	int build_block(block &b);
	int print(int index) const;

private:
	int init();
	int evaluate_rectangle(int ll, int rr, double &ave, double &dev);
	int evaluate_triangle(int ll, int rr, double &ave, double &dev);
};

#endif
