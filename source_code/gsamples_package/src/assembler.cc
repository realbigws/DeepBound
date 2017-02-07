#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "assembler.h"
#include "bundle.h"
#include "genome.h"

assembler::assembler(const string &file, const string &fasta_dir)
	: gf(file), fa(fasta_dir)
{
}

assembler::~assembler()
{
}

int assembler::process(const string &file, const string &sample_file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	sample_fout.open(sample_file.c_str());

	bundle_base bb1;		// for + reads
	bundle_base bb2;		// for - reads
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;
		if((p.flag & 0x4) >= 1) continue;			// read is not mapped, TODO
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(p.n_cigar < 1) continue;				// should never happen
		if(p.n_cigar > MAX_NUM_CIGAR) continue;	// ignore hits with more than 7 cigar types
		//if(p.qual <= 4) continue;				// ignore hits with quality-score < 5
		
		if(bb1.get_num_hits() > 0 && (bb1.get_rpos() + min_bundle_gap < p.pos || p.tid != bb1.get_tid()))
		{
			process_bundle(bb1, h);
		}

		if(bb2.get_num_hits() > 0 && (bb2.get_rpos() + min_bundle_gap < p.pos || p.tid != bb2.get_tid()))
		{
			process_bundle(bb2, h);
		}

		hit ht(b);
		if(ht.xs == '+') bb1.add_hit(ht);
		if(ht.xs == '-') bb2.add_hit(ht);
    }
	process_bundle(bb1, h);
	process_bundle(bb2, h);

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);

	sample_fout.close();

	return 0;
}

int assembler::process_bundle(bundle_base &bb, bam_hdr_t *h)
{
	if(bb.get_num_hits() < min_num_hits_in_bundle) 
	{
		bb.clear();
		return 0;
	}

	char buf[1024];
	int tid = bb.get_tid();
	assert(tid >= 0);
	strcpy(buf, h->target_name[tid]);
	string chrm(buf);
	bb.set_chrm(chrm);

	bundle bd(bb, fa);
	bd.build();

	set<int32_t> ss;
	set<int32_t> tt;
	join_interval_map jmap0;

	set<int32_t> &sss = ss;
	set<int32_t> &ttt = tt;
	join_interval_map &jmap = jmap0;

	//printf("bundle for chrm of %s\n", chrm.c_str());

	if(bd.strand == '+' && gf.mss1.find(chrm) != gf.mss1.end()) sss = gf.mss1[chrm];
	if(bd.strand == '-' && gf.mss2.find(chrm) != gf.mss2.end()) sss = gf.mss2[chrm];
	if(bd.strand == '+' && gf.mtt1.find(chrm) != gf.mtt1.end()) ttt = gf.mtt1[chrm];
	if(bd.strand == '-' && gf.mtt2.find(chrm) != gf.mtt2.end()) ttt = gf.mtt2[chrm];
	if(gf.jmap.find(chrm) != gf.jmap.end()) jmap = gf.jmap[chrm];

	bd.build_boundaries(sss, ttt);
	bd.assign_boundaries();

	for(int i = 0; i < bd.blocks.size(); i++)
	{
		block &b = bd.blocks[i];
		b.build_labels();
		b.build_abundance(jmap);
		b.build_features();
		b.write_samples(sample_fout);
		b.print_labels();
	}

	bb.clear();
	return 0;
}
