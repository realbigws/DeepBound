#ifndef __FASTA_H__
#define __FASTA_H__

#include <string>
#include <vector>

using namespace std;

class fasta
{
public:
	fasta(const string &d);

public:
	string dir;
	string chrm;
	vector<string> seq;

public:
	string get_seq(const string &chr, int s, int l);
	int load_chrm(const string &chr);
};

#endif
