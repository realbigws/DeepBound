#include "sample.h"
#include "config.h"
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <set>
#include <cmath>

int sample::add_position(const string &s)
{
	position p(s);
	positions.push_back(p);
	return 0;
}

int sample::clear()
{
	positions.clear();
	blocks0.clear();
	blocks1.clear();
	blocks2.clear();
	correct0 = correct1 = correct2 = 0;
	label0 = label1 = label2 = 0;
	accept = false;
	return 0;
}

int sample::process()
{
	build_blocks(0, blocks0);
	build_blocks(2, blocks2);
	align_blocks(0, blocks0, correct0, label0);
	align_blocks(2, blocks2, correct2, label2);
	predict0 = blocks0.size();
	predict2 = blocks2.size();
	build_blocks1();

	build_splice_positions();
	build_abundance();
	assess_abundance();

	qualify();
	return 0;
}

int sample::qualify()
{
	accept = true;
	return 0;

	if(label0 == 0 && label2 == 0) accept = false;

	for(int i = 0; i < positions.size(); i++)
	{
		assert(positions[i].tabd >= 0);
		if(positions[i].tabd < min_accept_expression) accept = false;
		if(accept == false) return 0;
	}
	
	return 0;
}

int sample::build_blocks1()
{
	blocks1.assign(positions.size(), 1);
	for(int i = 0; i < positions.size(); i++)
	{
		double pr = positions[i].xyz[1];
		if(pr >= min_prob) continue;
		int k1 = i - block_size / 2;
		int k2 = i + block_size / 2;
		if(k1 <= 0) k1 = 0;
		if(k2 >= positions.size() - 1) k2 = positions.size() - 1;
		for(int k = k1; k < k2; k++) blocks1[k] = 0;
	}

	correct1 = 0;
	predict1 = 0;
	label1 = 0;

	for(int i = 0; i < positions.size(); i++)
	{
		if(positions[i].tlab == 1) label1++;
		if(blocks1[i] == 1) predict1++;
		if(positions[i].tlab == 1 && blocks1[i] == 1) correct1++;
	}
	return 0;
}

int sample::build_blocks(int ff, vector<block> &blocks)
{
	blocks.clear();
	if(positions.size() < block_size) return 0;

	double pr = 0;
	int cnt = 0;
	for(int i = 0; i < block_size; i++)
	{
		position &p = positions[i];
		pr += p.xyz[ff];
		if(p.pred == ff) cnt++;
	}

	for(int i = 0; i < positions.size() - block_size; i++)
	{
		position &p1 = positions[i];
		position &p2 = positions[i + block_size];
		pr += p2.xyz[ff];
		pr -= p1.xyz[ff];
		if(p2.pred == ff) cnt++;
		if(p1.pred == ff) cnt--;
		assert(cnt <= block_size);

		double ave = pr * 1.0 / block_size;
		if(ave < min_prob) continue;
		//if(cnt < block_size) continue;

		blocks.push_back(block(i, ave));
	}

	sort(blocks.begin(), blocks.end());

	vector<block> vv;
	for(int i = 0; i < blocks.size(); i++)
	{
		block &b = blocks[i];
		bool f = false;
		for(int k = 0; k < vv.size(); k++)
		{
			if(b.distance(vv[k]) < min_block_distance) f = true;
			if(f == true) break;
		}
		if(f == true) continue;
		vv.push_back(b);
	}

	blocks = vv;
	return 0;
}

int sample::align_blocks(int ff, vector<block> &blocks, int &ncorrect, int &nlabel)
{
	vector<PI> vv;
	for(int i = 0; i < positions.size(); i++)
	{
		position &p = positions[i];
		if(p.tlab != ff) continue;
		vv.push_back(PI(i, -1));
	}
	nlabel = vv.size();
	ncorrect = 0;

	for(int i = 0; i < blocks.size(); i++)
	{
		block &b = blocks[i];
		vv.push_back(PI(b.pos + block_size / 2, i));
	}

	sort(vv.begin(), vv.end());

	set<int> s;
	for(int i = 0; i < vv.size(); i++)
	{
		int p = vv[i].first;
		int k = vv[i].second;
		if(k == -1) continue;

		int p1 = max_correct_distance * 2;
		int p2 = max_correct_distance * 2;

		if(i >= 1 && vv[i - 1].second == -1 && s.find(i - 1) == s.end()) p1 = p - vv[i - 1].first;
		if(i < vv.size() - 1 && vv[i + 1].second == -1 && s.find(i + 1) == s.end()) p2 = vv[i + 1].first - p;

		if(p1 < max_correct_distance && p1 <= p2)
		{
			blocks[k].match = p1;
			s.insert(i - 1);
			ncorrect++;
		}
		if(p2 < max_correct_distance && p2 < p1)
		{
			blocks[k].match = p2;
			s.insert(i + 1);
			ncorrect++;
		}
	}
	return 0;
}

int sample::build_splice_positions()
{
	splist.clear();
	for(int i = 0; i < blocks0.size(); i++)
	{
		int p = blocks0[i].pos + block_size / 2;
		splist.push_back(p);
	}
	for(int i = 0; i < blocks2.size(); i++)
	{
		int p = blocks2[i].pos + block_size / 2;
		splist.push_back(p);
	}
	
	splist.push_back(positions.size());
	sort(splist.begin(), splist.end());

	return 0;
}

int sample::build_abundance()
{
	vabd.clear();
	int p1 = 0;
	for(int i = 0; i < splist.size(); i++)
	{
		int p2 = splist[i];
		assert(p2 >= p1);
		if(p1 == p2) continue;

		double ave = 0;
		for(int k = p1; k < p2; k++)
		{
			if(extended == false) ave += positions[k].rabd;
			else ave += positions[k].pabd;
		}
		ave = ave / (p2 - p1);

		for(int k = p1; k < p2; k++) vabd.push_back(ave);

		p1 = p2;
	}
	assert(positions.size() == vabd.size());
	return 0;
}

int sample::assess_abundance()
{
	abdratio = 0;
	for(int i = 0; i < positions.size(); i++)
	{
		double tabd = positions[i].tabd + 1.0;
		double pabd = vabd[i] + 1.0;
		abdratio += fabs(tabd - pabd) / tabd;
	}
	abdratio /= vabd.size();
	return 0;
}

int sample::print(int index)
{
	printf("sample %d: positions = %lu, abd-ratio = %.3lf, label0 = %d / %d / %d, label2 = %d / %d / %d, label1 = %d / %d / %d\n", 
			index,
			positions.size(), 
			abdratio,
			correct0, predict0, label0, 
			correct2, predict2, label2, 
			correct1, predict1, label1);
	return 0;

	// print blocks
	for(int i = 0; i < blocks0.size(); i++)
	{
		printf("block0 %d (size = %lu): ", i, positions.size());
		blocks0[i].print();
		int p = blocks0[i].pos;
		positions[p].print();
	}
	for(int i = 0; i < blocks2.size(); i++)
	{
		printf("block2 %d (size = %lu): ", i, positions.size());
		blocks2[i].print();
		int p = blocks2[i].pos;
		positions[p].print();
	}
	return 0;
}

int sample::print_result(int index)
{
	for(int i = 0; i < blocks0.size(); i++)
	{
		int p = blocks0[i].pos;
		printf("sample %d: end boundary at position %d\n", index, p);
	}

	for(int i = 0; i < blocks2.size(); i++)
	{
		int p = blocks2[i].pos;
		printf("sample %d: start boundary at position %d\n", index, p);
	}

	//printf("sample %d: summary end boundaries = %d / %d / %d\n", index, correct0, predict0, label0); 
	//printf("sample %d: summary start boundaries = %d / %d / %d\n", index, correct2, predict2, label2); 

	return 0;
}

int sample::write(ofstream &fout)
{
	set<int> s0;
	set<int> s2;
	for(int i = 0; i < blocks0.size(); i++)
	{
		int p = blocks0[i].pos;
		s0.insert(p);
	}

	for(int i = 0; i < blocks2.size(); i++)
	{
		int p = blocks2[i].pos;
		s2.insert(p);
	}

	fout << header.c_str() << endl;
	for(int i = 0; i < positions.size(); i++)
	{
		if(s0.find(i) != s0.end()) fout << "0 ";
		else if(s2.find(i) != s2.end()) fout << "2 ";
		else fout << "1 ";
	}
	fout << endl;

	//printf("sample %d: summary end boundaries = %d / %d / %d\n", index, correct0, predict0, label0); 
	//printf("sample %d: summary start boundaries = %d / %d / %d\n", index, correct2, predict2, label2); 

	return 0;
}
