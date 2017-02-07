#include "fscore.h"
#include <cstdio>
#include <iomanip>

int fscore::init(int k)
{
	sa.assign(k, -1);
	sb.assign(k, -1);
	va.assign(k, -1);
	vb.assign(k, -1);
	return 0;
}

int fscore::write(ofstream &fout)
{
	fout << setiosflags(ios::fixed) << setprecision(3);
	for(int i = 0; i < sa.size(); i++) fout << sa[i] << " "; fout << "\n";
	for(int i = 0; i < sb.size(); i++) fout << sb[i] << " "; fout << "\n";
	for(int i = 0; i < va.size(); i++) fout << va[i] << " "; fout << "\n";
	for(int i = 0; i < vb.size(); i++) fout << vb[i] << " "; fout << "\n";
	return 0;
}
