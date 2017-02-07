#include <cstdio>
#include <cassert>
#include <sstream>

#include "gtf.h"
#include "config.h"

gtf::gtf(const string &file)
	: gm(file)
{
	build_boundary_positions();
	build_interval_map();
	print();
}

gtf::~gtf()
{
}

int gtf::build_boundary_positions()
{
	mss1.clear();
	mss2.clear();
	mtt1.clear();
	mtt2.clear();

	for(int i = 0; i < gm.genes.size(); i++)
	{
		gene &g = gm.genes[i];
		string chrm = g.get_seqname();
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript &t = g.transcripts[k];
			PI32 p = t.get_bounds();
			if(p.first < 0 || p.second < 0) continue;

			//printf("transcript %d [%d-%d), expression = %d\n", k, p.first, p.second, t.expression);

			if(t.coverage < min_transcript_expression) continue;

			if(t.strand == '+')
			{
				if(mss1.find(chrm) == mss1.end())
				{
					set<int32_t> s;
					s.insert(p.first);
					mss1.insert(PSSI(chrm, s));
				}
				else
				{
					mss1[chrm].insert(p.first);
				}

				if(mtt1.find(chrm) == mtt1.end())
				{
					set<int32_t> s;
					s.insert(p.second);
					mtt1.insert(PSSI(chrm, s));
				}
				else
				{
					mtt1[chrm].insert(p.second);
				}
			}
			else if(t.strand == '-')
			{
				if(mss2.find(chrm) == mss2.end())
				{
					set<int32_t> s;
					s.insert(p.first);
					mss2.insert(PSSI(chrm, s));
				}
				else
				{
					mss2[chrm].insert(p.first);
				}

				if(mtt2.find(chrm) == mtt2.end())
				{
					set<int32_t> s;
					s.insert(p.second);
					mtt2.insert(PSSI(chrm, s));
				}
				else
				{
					mtt2[chrm].insert(p.second);
				}
			}
		}
	}

	return 0;
}

int gtf::build_interval_map()
{
	jmap.clear();
	for(int i = 0; i < gm.genes.size(); i++)
	{
		gene &g = gm.genes[i];
		string chrm = g.get_seqname();
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript &t = g.transcripts[k];
			int abd = t.coverage;

			for(int e = 0; e < t.exons.size(); e++)
			{
				int p1 = t.exons[e].first;
				int p2 = t.exons[e].second;

				if(jmap.find(chrm) == jmap.end())
				{
					join_interval_map m;
					m += make_pair(ROI(p1, p2), abd);
					jmap.insert(PSJIM(chrm, m));
				}
				else
				{
					jmap[chrm] += make_pair(ROI(p1, p2), abd);
				}
			}
		}
	}

	return 0;
}

int gtf::print()
{
	int cnt1 = 0;
	int cnt2 = 0;
	for(MSSI::iterator it = mss1.begin(); it != mss1.end(); it++)
	{
		printf("#map of chrm %s has %lu start positions on positive strands\n", it->first.c_str(), it->second.size());
		cnt1 += it->second.size();
	}
	for(MSSI::iterator it = mss2.begin(); it != mss2.end(); it++)
	{
		printf("#map of chrm %s has %lu start positions on negative strands\n", it->first.c_str(), it->second.size());
		cnt1 += it->second.size();
	}
	for(MSSI::iterator it = mtt1.begin(); it != mtt1.end(); it++)
	{
		printf("#map of chrm %s has %lu end positions on positive strands\n", it->first.c_str(), it->second.size());
		cnt2 += it->second.size();
	}
	for(MSSI::iterator it = mtt2.begin(); it != mtt2.end(); it++)
	{
		printf("#map of chrm %s has %lu end positions on positive strands\n", it->first.c_str(), it->second.size());
		cnt2 += it->second.size();
	}
	printf("#total %d start boundaries, %d end boundaries\n", cnt1, cnt2);
	return 0;
}
