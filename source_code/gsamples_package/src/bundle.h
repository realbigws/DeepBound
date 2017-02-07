#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "interval_map.h"
#include "bundle_base.h"
#include "junction.h"
#include "region.h"
#include "block.h"
#include "fasta.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb, fasta &_fa);
	virtual ~bundle();

public:
	fasta &fa;						// fasta
	vector<junction> junctions;		// splice junctions
	vector<region> regions;			// regions
	split_interval_map imap;		// interval map
	split_interval_map qmap;		// quality map
	vector<block> blocks;			// blocks

	split_interval_map bmap;		// block map
	vector<int32_t> mss;			// ground-truth start-positions
	vector<int32_t> mtt;			// ground-truth end-positions

public:
	int build();
	int print(int index) const;
	int build_boundaries(const set<int32_t> &ss, const set<int32_t> &tt);
	int assign_boundaries();

protected:
	int check_left_ascending();
	int check_right_ascending();
	int infer_junctions();
	int process_single_hits();
	int build_regions();
	int build_blocks();
	int build_block_map();
	int locate_block(int32_t x);
};

#endif
