#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>

#include "assembler.h"
#include "config.h"
#include "block.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		printf("usage: %s <in-bam-file> <out-sample-file> [-f fasta-dir] [-g gtf-file] [-e min-expression]\n", argv[0]);
		return 0;
	}

	parse_arguments(argc, argv);
	block::index = index_prefix;

	assembler asmbl(input_gtf_file, input_fasta_dir);
	asmbl.process(argv[1], argv[2]);

    return 0;
}
