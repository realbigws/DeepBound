#include "region.h"
#include "config.h"
#include "util.h"
#include <algorithm>

using namespace std;

region::region(const string &_chrm, int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap, const split_interval_map *_qmap, fasta &_fa)
	:chrm(_chrm), lpos(_lpos), rpos(_rpos), imap(_imap), qmap(_qmap), ltype(_ltype), rtype(_rtype), fa(_fa)
{
	lcore = lpos;
	rcore = rpos;
}

region::~region()
{}

int region::build_block(block &b)
{
	b.clear();

	init();
	if(empty == true) return 0;

	b.chrm = chrm;
	b.pos = lcore;
	b.s.clear();
	b.q.clear();

	for(int i = lcore; i < rcore; i++)
	{
		int32_t w = compute_overlap(*imap, i);
		double q = compute_overlap(*qmap, i);
		if(w <= 0) q = -1;
		else q = q * 1.0 / w; 
		b.s.push_back(w);
		b.q.push_back(q);
	}

	/*
	if(ltype == START_BOUNDARY) b.ltype = true;
	else b.ltype = false;
	if(rtype == END_BOUNDARY) b.rtype = true;
	else b.rtype = false;
	*/


	b.ltype = ltype;
	b.rtype = rtype;
	b.seq = fa.get_seq(chrm, lcore, rcore - lcore);

	return 0;
}

int region::init()
{
	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);

	if(lit == imap->end() || rit == imap->end())
	{
		empty = true;
		return 0;
	}

	lcore = lower(lit->first);
	rcore = upper(rit->first);
	assert(rcore > lcore);

	if(lcore > lpos) ltype = START_BOUNDARY;
	if(rcore < rpos) rtype = END_BOUNDARY;

	empty = false;
	int32_t m = compute_max_overlap(*imap, lit, rit);
	int32_t s = compute_sum_overlap(*imap, lit, rit);
	if(m < min_max_region_overlap) empty = true;
	if(1.0 * s / (rcore - lcore) < min_average_overlap) empty = true;

	return 0;
}


int region::evaluate_rectangle(int ll, int rr, double &ave, double &dev)
{
	ave = 0;
	dev = 1.0;

	if(empty == true) return 0;

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, ll, rr);

	if(lit == imap->end()) return 0;
	if(rit == imap->end()) return 0;

	ave = 1.0 * compute_sum_overlap(*imap, lit, rit) / (rr - ll);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave) * (it->second - ave) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	if(dev < 1.0) dev = 1.0;

	return 0;
}

int region::evaluate_triangle(int ll, int rr, double &ave, double &dev)
{
	ave = 0;
	dev = 1.0;

	if(empty == true) return 0;

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, ll, rr);

	if(lit == imap->end()) return 0;
	if(rit == imap->end()) return 0;

	vector<double> xv;
	vector<double> yv;
	double xm = 0;
	double ym = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (lower(it->first) + upper(it->first)) / 2.0;
		double yi = it->second;
		xv.push_back(xi);
		yv.push_back(yi);
		xm += xi;
		ym += yi;
		if(it == rit) break;
	}

	xm /= xv.size();
	ym /= yv.size();

	double f1 = 0;
	double f2 = 0;
	for(int i = 0; i < xv.size(); i++)
	{
		f1 += (xv[i] - xm) * (yv[i] - ym);
		f2 += (xv[i] - xm) * (xv[i] - xm);
	}

	double b1 = f1 / f2;
	double b0 = ym - b1 * xm;

	double a1 = b1 * rr + b0;
	double a0 = b1 * ll + b0;
	ave = (a1 > a0) ? a1 : a0;

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (upper(it->first) + lower(it->first)) / 2.0;
		double yi = b1 * xi + b0;
		var += (it->second - yi) * (it->second - yi) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	if(dev < 1.0) dev = 1.0;

	return 0;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*imap, lcore);
	int32_t rc = compute_overlap(*imap, rcore - 1);
	printf("region %d: empty = %c, pos = [%d, %d), core = [%d, %d), coverage = (%d, %d)\n", index, empty ? 'T' : 'F', lpos, rpos, lcore, rcore, lc, rc);

	return 0;
}
