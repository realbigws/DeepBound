#!/bin/bash

# ----- usage ------ #
usage()
{
	echo "DeepBound_Package v0.01 [Feb-10-2017] "
	echo "    Accurate Identification of Transcript Boundaries via Deep Convolutional Neural Fields (DeepCNF)"
	echo ""
	echo "USAGE: ./oneline_command.sh <-i input_BAM_file> [-o output_boundary] "
	echo "        [-f chromosomes_dir] [-g annotation_file] [-m min_expression] "
	echo "        [-O premature_bound] [-p prob_threshold] [-w window_size]  "
	echo ""
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_BAM_file  : specifies the input reads alignemnt file in BAM format. "
	echo "                     It has to be sorted (for example, by samtools). "
	echo "-o output_boundary : specifies the predicted boundaries in FASTA format. "
	echo "                     0 for end-boundary, 1 for non-boundary, and 2 for start-boundary. "
	echo ""
	echo "***** optional arguments *****"
	echo "--| Step 1. from BAM to generate ReadCount by bin/gsamples"
	echo "-f chromosomes_dir : is optional but strongly recommended, "
	echo "                     which specifies directory of chromosome sequences in FASTA format. "
	echo "-g annotation_file : is optional, which specifies the ground-truth expression abundance in GTF file. "
	echo "-m min_expression  : is together with '-g' option, which will make gsamples "
	echo "                     ignore these transcripts in GTF file with expression abundance under this value."
	echo ""
	echo "--| Step 2. predict premature boundary by DeepBound.sh"
	echo "-O premature_bound : is the directory containing predicted premature boundary, "
	echo "                     by DeepBound.sh which is based on AUC-maximizd DeepCNF model. "
	echo ""
	echo "--| Step 3. process premature boundary by bin/pinpoint"
	echo "-p prob_threshold  : is optional (default 0.6), which specifies the minimum average probability "
	echo "                     over a window that will be considerred as a true boundary. "
	echo "-w window_size     : is optional (default 10), which gives the window size for average calculation. "
	exit 1
}

if [ $# -lt 1 ];
then
	usage
fi
curdir="$(pwd)"

# ----- get arguments ----- #
#-> required arguments
input_BAM_file=""
output_boundary_prediction=""
#-> optional arguments
#--| for bin/gsamples
chromosomes_fasta_dir=""
annotation_gtf_file=""
min_expression=""
#--| for DeepBound.sh
premature_bound=""
#--| for bin/pinpoint
prob_threshold=0.6
window_size=10

#-> parse arguments
while getopts ":i:o:f:g:m:O:p:w:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input_BAM_file=$OPTARG
		;;
	o)
		output_boundary_prediction=$OPTARG
		;;
	#-> optional arguments
	#--| for bin/gsamples
	f)
		chromosomes_fasta_dir=$OPTARG
		;;
	g)
		annotation_gtf_file=$OPTARG
		;;
	m)
		min_expression=$OPTARG
		;;
	#--| for DeepBound.sh
	O)
		premature_bound=$OPTARG
		;;
	#--| for bin/pinpoint
	p)
		prob_threshold=$OPTARG
		;;
	w)
		window_size=$OPTARG
		;;
	#-> exception
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done


#----- check input BAM file ---#
if [ -z "$input_BAM_file" ]
then
	echo "input_BAM_file is null !!"
	exit 1
fi
if [ ! -f "$curdir/$input_BAM_file" ]
then
	if [ ! -f "$input_BAM_file" ]
	then
		echo "input_BAM_file $input_BAM_file not found !!"
		exit 1
	fi
else
	input_BAM_file=$curdir/$input_BAM_file
fi

# ------ get name ----- #
fulnam=`basename $input_BAM_file`
relnam=${fulnam%.*}
if [ "$premature_bound" == "" ]
then
	premature_bound=${relnam}_out
fi
if [ "$output_boundary_prediction" == "" ]
then
	output_boundary_prediction=${relnam}.boundary
fi
mkdir -p $premature_bound


#======== Section I: from BAM, generate ReadCount by bin/gsamples ==============#
#-> determine parameters
chromosomes_fasta_dir_=""
if [ "$chromosomes_fasta_dir" != "" ]
then
	chromosomes_fasta_dir_=" -f "$chromosomes_fasta_dir
fi
annotation_gtf_file_=""
if [ "$annotation_gtf_file" != "" ]
then
	annotation_gtf_file_=" -g "$annotation_gtf_file
fi
min_expression_=""
if [ "$min_expression" != "" ]
then
	min_expression_=" -e "$min_expression
fi

#-> run bin/gsamples
bin/gsamples $input_BAM_file $premature_bound/${relnam}.ReadCount $chromosomes_fasta_dir_ $annotation_gtf_file_ $min_expression_


#======== Section II: predict premature boundary by DeepBound.sh ==================#
./DeepBound.sh $premature_bound/${relnam}.ReadCount $premature_bound


#======== Section III: process premature boundary prediction by bin/pinpoint ===========#
bin/pinpoint $premature_bound/${relnam}.ReadCount $premature_bound/${relnam}.premature $output_boundary_prediction -p $prob_threshold -w $window_size

