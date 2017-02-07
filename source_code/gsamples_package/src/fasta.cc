#include "fasta.h"
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

fasta::fasta(const string &d)
	: dir(d)
{
	chrm = "shaomingfu";
}

string fasta::get_seq(const string &chr, int s, int l)
{
	if(dir == "") return "";

	if(chr != chrm)
	{
		load_chrm(chr);
		return get_seq(chr, s, l);
	}

	string ss;
	
	if(seq.size() == 0) return ss;

	int k = s / seq[0].size(); 
	int a = k * seq[0].size();
	int t = s + l;
	for(int i = k; i < seq.size(); i++)
	{
		int b = a + seq[i].size();

		if(a >= t) break;

		if(b <= s) 
		{
			a = b;
			continue;
		}

		if(a >= s && t >= b)
		{
			ss += seq[i];
			a = b;
			continue;
		}

		int m = (a > s ? a : s) - a;
		int n = (t < b ? t : b) - a - m;

		assert(m >= 0);
		assert(n >= 0);

		ss += seq[i].substr(m, n);

		a = b;
	}

	return ss;
}

int fasta::load_chrm(const string &chr)
{
	if(dir == "") return 0;
	if(chr == chrm) return 0;

	string file = dir + "/" + chr + ".fa";

	ifstream fin(file.c_str());

	if(fin.fail())
	{
		printf("open fasta file error %s\n", file.c_str());
		assert(false);
	}

	seq.clear();
	string s;
	while(getline(fin, s, '\n'))
	{
		if(s.size() == 0) continue;
		if(s[0] == '>') continue;
		seq.push_back(s);
	}
	fin.close();

	for(int k = 1; k < seq.size() - 1; k++) assert(seq[k].size() == seq[k - 1].size());

	chrm = chr;
	return 0;
}
