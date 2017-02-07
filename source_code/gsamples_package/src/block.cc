#include "block.h"
#include "config.h"
#include <cassert>
#include <binomial.h>

int block::index = 0;

int block::clear()
{
	chrm = "";
	pos = -1;
	s.clear();
	q.clear();
	ss.clear();
	tt.clear();
	return 0;
}

int block::evaluate(int a, int b, double &ave, double &dev)
{
	ave = dev = 0;
	if(a >= b) return 0;

	double sum = 0;
	for(int i = a; i < b; i++)
	{
		//printf("a = %d, b = %d, i = %d, s.size() = %lu\n", a, b, i, s.size());
		assert(i >= 0 && i < s.size());
		sum += s[i];
	}

	ave = sum / (b - a);

	double var = 0;
	for(int i = a; i < b; i++) var += (s[i] - ave) * (s[i] - ave);
	dev = sqrt(var / (b - a));

	return 0;
}

int block::build_feature_score(fscore &fs)
{
	fs.init(s.size());
	for(int i = fs.w; i < s.size() - fs.w; i++)
	{
		double ave1, ave2, dev1, dev2;
		evaluate(i - fs.w, i, ave1, dev1);
		evaluate(i, i + fs.w, ave2, dev2);

		int n1 = ave1 * fs.w / 75.0 + 1;
		int n2 = ave2 * fs.w / 75.0 + 1;
		fs.sa[i] = compute_binomial_score(n1 + n2, 0.5, n1);
		fs.sb[i] = compute_binomial_score(n1 + n2, 0.5, n2);
		fs.va[i] = dev1 / (ave1 + 1.0);
		fs.vb[i] = dev2 / (ave2 + 1.0);
	}
	return 0;
}

int block::build_abundance(const join_interval_map &jmap)
{
	if(s.size() < min_sample_length) return 0;

	abd.clear();
	for(int i = 0; i < s.size(); i++)
	{
		int32_t p = i + pos;
		JIMI it = jmap.find(p);
		if(it == jmap.end()) abd.push_back(0);
		else abd.push_back(it->second);
		int l = locate_label(abd[i]);
		abl.push_back(l);
	}
	return 0;
}

int block::build_labels()
{
	if(s.size() < min_sample_length) return 0;
	labels.assign(s.size(), 1);
	for(int i = 0; i < ss.size(); i++)
	{
		int k = ss[i] - pos;
		if(k < 0) labels[0] = 2;
		else if(k < s.size()) labels[k] = 2;
	}
	for(int i = 0; i < tt.size(); i++)
	{
		int k = tt[i] - pos;
		if(k >= s.size()) labels[s.size() - 1] = 0;
		else if(k >= 0) labels[k] = 0;
	}
	return 0;
}

int block::build_labels(const set<int32_t> &ss, const set<int32_t> &tt)
{
	if(s.size() < min_sample_length) return 0;

	labels.clear();
	for(int i = 0; i < s.size(); i++)
	{
		int32_t p = i + pos;
		int label = 1;
		if(ss.find(p) != ss.end()) label = 2;
		if(tt.find(p) != tt.end()) label = 0;
		/* TODO important
		if(i == 0 && ltype == true) label = 2;
		if(i == s.size() - 1 && rtype == true) label = 0;
		*/
		labels.push_back(label);
	}
	return 0;
}

int block::build_features()
{
	if(s.size() < min_sample_length) return 0;

	fs20.w = 20;
	fs50.w = 50;
	fs100.w = 100;
	build_feature_score(fs20);
	build_feature_score(fs50);
	build_feature_score(fs100);

	return 0;
}

bool block::qualify_boundary_training()
{
	if(s.size() < min_sample_length) return false;
	else return true;
}

bool block::qualify_abundance_training()
{
	if(qualify_boundary_training() == false) return false;

	bool b = false;
	for(int i = 0; i < labels.size(); i++)
	{
		if(labels[i] == 0) b = true;
		if(labels[i] == 2) b = true;
		if(b == true) break;
	}

	if(b == false) return false;

	int cnt = 0;
	for(int i = 0; i < abd.size(); i++)
	{
		if(abd[i] >= min_transcript_expression) cnt++;
	}
	
	if(cnt * 1.0 / abd.size() < 0.8) return false;
	else return true;
}

int block::write_index()
{
	assert(s.size() == q.size());

	if(qualify_boundary_training() == false) return 0;

	block::index++;

	if(qualify_abundance_training() == false) return 0;

	printf("abundance-training-index %d\n", index);

	return 0;
}

int block::write_samples(ofstream &fout)
{
	assert(s.size() == q.size());

	if(qualify_boundary_training() == false) return 0;

	block::index++;
	//printf("# sample-id = %d, length = %lu, location = %s:%d-%lu\n", block::index, s.size(), chrm.c_str(), pos, pos + s.size());

	fout << setiosflags(ios::fixed) << setprecision(3);
	fout << "# sample-id = " << block::index <<", length = " << s.size();
	fout << ", location = " << chrm.c_str() << ":" << pos << "-" << pos + s.size() << "\n";

	// print label
	for(int i = 0; i < labels.size(); i++)  fout << labels[i] <<" "; fout << "\n";

	// print abundance
	for(int i = 0; i < abd.size(); i++) fout<< abd[i] << " "; fout<<"\n";
	//for(int i = 0; i < labels.size(); i++)  fout << abl[i] <<" "; fout << "\n";

	// print features
	for(int i = 0; i < s.size(); i++) fout<< s[i] << " "; fout<<"\n";
	for(int i = 0; i < q.size(); i++) fout<< q[i] << " "; fout<<"\n";

	fs20.write(fout);
	fs50.write(fout);
	fs100.write(fout);

	// print abundance
	//for(int i = 0; i < abd.size(); i++) printf("%d ", abd[i]); printf("\n");

	// write sequence
	//assert(seq.size() == s.size());
	for(int i = 0; i < s.size(); i++)
	{
		if(i < seq.size() && seq[i] == 'A') fout << "1 "; 
		else fout << "0 ";
	}
	fout << "\n";

	for(int i = 0; i < s.size(); i++)
	{
		if(i < seq.size() && seq[i] == 'C') fout << "1 "; 
		else fout << "0 ";
	}
	fout << "\n";

	for(int i = 0; i < s.size(); i++)
	{
		if(i < seq.size() && seq[i] == 'G') fout << "1 "; 
		else fout << "0 ";
	}
	fout << "\n";

	for(int i = 0; i < s.size(); i++)
	{
		if(i < seq.size() && seq[i] == 'T') fout << "1 "; 
		else fout << "0 ";
	}
	fout << "\n";

	// write boundaries features
	if(ltype == START_BOUNDARY) fout<< "1 ";
	else fout<<"0 ";
	for(int i = 1; i < s.size(); i++) fout<<"0 ";
	fout << "\n";

	for(int i = 1; i < s.size(); i++) fout<<"0 ";
	if(rtype == END_BOUNDARY) fout<< "1 ";
	else fout<<"0 ";
	fout << "\n";

	// write abundance training index
	if(qualify_abundance_training() == false) return 0;
	//printf("abundance-training-index %d\n", index);

	return 0;
}

int block::print_labels()
{
	for(int i = 0; i < labels.size(); i++)
	{
		if(labels[i] == 0) printf("block %s:%d-%lu, end boundary at %d\n", chrm.c_str(), pos, pos + s.size(), i + pos);
		if(labels[i] == 2) printf("block %s:%d-%lu, start boundary at %d\n", chrm.c_str(), pos, pos + s.size(), i + pos);
	}
	return 0;
}
