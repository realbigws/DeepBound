#ifndef __GTF_H__
#define __GTF_H__

#include <fstream>
#include <string>
#include <set>

#include "genome.h"
#include "interval_map.h"

using namespace std;

typedef pair< string, set<int32_t> > PSSI;
typedef map< string, set<int32_t> > MSSI;
typedef pair< string, join_interval_map> PSJIM;
typedef map< string, join_interval_map> MSJIM;

class gtf
{
public:
	gtf(const string &file);
	~gtf();

public:
	genome gm;
	MSSI mss1;		// positive strand
	MSSI mss2;		// negative strand
	MSSI mtt1;
	MSSI mtt2;
	MSJIM jmap;

protected:
	int build_boundary_positions();
	int build_interval_map();
	int print();
};

#endif
