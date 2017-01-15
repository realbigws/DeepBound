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
#-> abundance label model
abundance_model_param="-w 75 -d 20 -s 15 -l 21"
abundance_model_file="abu_model_MaxLik_15lab"
#-> abundance value model
abundance_val_model="abu_value_model_MSE_relval"

#---- oneline command ---#
bin=bin/
mkdir -p $output_folder/feat_out/
mkdir -p $output_folder/cov_out/
mkdir -p $output_folder/bou_out/
mkdir -p $output_folder/abu_out/
mkdir -p $output_folder/abu_val/
mkdir -p $output_folder/bou_pred_out/
mkdir -p $output_folder/abu_feat_out/
mkdir -p $output_folder/abu_pred_out/
mkdir -p $output_folder/abu_feat_val/
mkdir -p $output_folder/abu_pred_val/
mkdir -p $output_folder/abu_pred_fin/


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


exit

#-> step 2. predict abundance label
#--> generate feature
rm -f ${relnam}_abu_feat_proc
for ((i=0;i<$totnum;i++))
do
	echo "$bin/Feature_Merge $output_folder/feat_out/sample_$i $output_folder/bou_pred_out/sample_$i $output_folder/abu_out/sample_$i > $output_folder/abu_feat_out/sample_$i" >> ${relnam}_abu_feat_proc
done
$bin/distribute_ii_openmp ${relnam}_abu_feat_proc
rm -f ${relnam}_abu_feat_proc
$bin/Fast_CAT $output_folder/${relnam}_sample_list $output_folder/abu_feat_out $output_folder/$relnam.abu_feat

#--> predict feature
$bin/DeepCNF_Pred -i $output_folder/$relnam.abu_feat $abundance_model_param -m parameters/$abundance_model_file > $output_folder/$relnam.abu_pred
$bin/Decompose_Commented_File $output_folder/$relnam.abu_pred $output_folder/abu_pred_out/ sample
rm -f $output_folder/$relnam.abu_feat
rm -f $output_folder/$relnam.abu_pred


#-> step 3. predict abundance value
#--> generate feature
rm -f ${relnam}_abu_feat_val_proc
for ((i=0;i<$totnum;i++))
do
	echo "$bin/Feature_Merge_wind $output_folder/bou_pred_out/sample_$i $output_folder/abu_pred_out/sample_$i $output_folder/cov_out/sample_$i $output_folder/abu_val/sample_$i > $output_folder/abu_feat_val/sample_$i" >> ${relnam}_abu_feat_val_proc
done
$bin/distribute_ii_openmp ${relnam}_abu_feat_val_proc
rm -f ${relnam}_abu_feat_val_proc

#--> predict 
rm -f ${relnam}_abu_pred_val_proc
for ((i=0;i<$totnum;i++))
do
	echo "$bin/nnreg -iFile $output_folder/abu_feat_val/sample_$i -model parameters/$abundance_val_model -act predict 1> $output_folder/abu_pred_val/sample_$i 2> ws2" >> ${relnam}_abu_pred_val_proc
done
$bin/distribute_ii_openmp ${relnam}_abu_pred_val_proc
rm -f ${relnam}_abu_pred_val_proc

#--> final predict
rm -f ${relnam}_abu_pred_fin_proc
for ((i=0;i<$totnum;i++))
do
	echo "$bin/Abundance_Estimate $output_folder/abu_pred_out/sample_$i $output_folder/abu_pred_val/sample_$i $output_folder/cov_out/sample_$i > $output_folder/abu_pred_fin/sample_$i" >> ${relnam}_abu_pred_fin_proc
done
$bin/distribute_ii_openmp ${relnam}_abu_pred_fin_proc
rm -f ${relnam}_abu_pred_fin_proc


#-> step 4. remove unnesessary files/folders
rm -rf $output_folder/abu_feat_val/
rm -rf $output_folder/abu_feat_out/
rm -rf $output_folder/abu_pred_val/
mv $output_folder/abu_pred_fin $output_folder/abu_pred_val
rm -f ws1 ws2

