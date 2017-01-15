#!/bin/bash

if [  $# -lt  2 ]
then
	echo "Usage: ./oneline_command.sh <original_input> <output_folder> "
	echo "[Note]:    <original_input> should be concatenated into one file"
	echo "           <output_folder> shall contain output folders"
	echo "           see README file for detail description of output folders"
	exit 1
fi

#---- input arguments ---#
original_input=$1
output_folder=$2
fulnam=`basename $original_input`
relnam=${fulnam%.*}

#---- additional parameters ---#
#-> boundary label model
boundary_model_param="-w 75 -d 100 -s 3 -l 20"
boudanry_model_file="bou_model_MaxAUC_3lab"

#---- oneline command ---#
bin=bin/
mkdir -p $output_folder/feat_out/
mkdir -p $output_folder/cov_out/
mkdir -p $output_folder/bou_out/
mkdir -p $output_folder/abu_out/
mkdir -p $output_folder/abu_val/
mkdir -p $output_folder/bou_pred_out/


#-> step 1. predict boundary
#--> generate feature
feat_dim=20
$bin/Generate_Feature $original_input $output_folder/feat_out/ $output_folder/cov_out/ $output_folder/bou_out/ $output_folder/abu_val/ $output_folder/abu_out/ $feat_dim
totnum=`ls $output_folder/feat_out/ | wc | awk '{print $1}'`
rm -f ${relnam}_sample_cat
rm -f $output_folder/${relnam}_sample_list
for ((i=0;i<$totnum;i++))
do
	echo "$output_folder/feat_out/sample_$i" >> ${relnam}_sample_cat
	echo "$output_folder/bou_out/sample_$i" >> ${relnam}_sample_cat
	echo "sample_$i" >> $output_folder/${relnam}_sample_list
done
$bin/Fast_CAT_v2 ${relnam}_sample_cat $output_folder/$relnam.bou_feat
rm -f ${relnam}_sample_cat


#--> predict feature
$bin/DeepCNF_Pred -i $output_folder/$relnam.bou_feat $boundary_model_param -m parameters/$boudanry_model_file > $output_folder/$relnam.bou_pred
$bin/Decompose_Commented_File $output_folder/$relnam.bou_pred $output_folder/bou_pred_out/ sample
rm -f $output_folder/$relnam.bou_feat
rm -f $output_folder/$relnam.bou_pred

#--> remove temporary folders
rm -rf $output_folder/cov_out/ $output_folder/abu_out/ $output_folder/abu_val/ $output_folder/feat_out/ $output_folder/bou_out/

