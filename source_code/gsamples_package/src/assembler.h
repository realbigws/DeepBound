#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>
#include <set>
#include <fstream>

#include "gtf.h"
#include "block.h"
#include "bundle_base.h"
#include "fasta.h"

using namespace std;

class assembler
{
public:
	assembler(const string &gtf_file, const string &fasta_dir);
	~assembler();

public:
	fasta fa;
	gtf gf;
	ofstream sample_fout;

public:
	int build_boundary_positions(const string &file);
	int process_bundle(bundle_base &bb, bam_hdr_t *h);
	int process(const string &file, const string &sample_file);
};

#endif
