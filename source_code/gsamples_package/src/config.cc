#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters

int32_t max_correct_distance = 50;
int32_t min_bundle_gap = 50;
int min_num_hits_in_bundle = 20;
int32_t min_splice_boundary_hits = 1;
uint32_t min_max_splice_boundary_qual = 0;
double min_average_overlap = 2;
int min_max_region_overlap = 5;
int min_sample_length = 100;
double min_region_coverage = 5.0;

string input_gtf_file = "";
string input_fasta_dir = "";
double min_transcript_expression = 1.0;
int index_prefix = 0;

//vector<double> abundance_labels = {15, 22, 33, 48, 68, 95, 132, 181, 255, 362, 530, 932};		// min_transcript_expression = 10
vector<double> abundance_labels = {1, 2, 3, 4, 6, 9, 15, 22, 33, 47, 68, 94, 130, 179, 252, 358, 524, 859};	// min_transcript_expression = 1.0

int locate_label(double abd)
{
	for(int k = 0; k < abundance_labels.size(); k++)
	{
		if(abd <= abundance_labels[k]) return k;
	}
	return abundance_labels.size();
}

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-g")
		{
			input_gtf_file = string(argv[i + 1]);
			i++;
		}
		if(string(argv[i]) == "-f")
		{
			input_fasta_dir = string(argv[i + 1]);
			i++;
		}
		if(string(argv[i]) == "-e")
		{
			min_transcript_expression = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--index_prefix")
		{
			index_prefix = atoi(argv[i + 1]);
			i++;
		}
	}
	return 0;
}
