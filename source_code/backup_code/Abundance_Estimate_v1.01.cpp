#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
using namespace std;


//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//---- parse feat_num ---//
int Parse_FeatNum(string &in_str)
{
	int feat_num=0;
	istringstream www(in_str);
	for(;;)
	{
		string temp;
		if(! (www>>temp) )break;
		feat_num++;
	}
	return feat_num;
}

//---------- input a string, output a vector -----//
int String_To_Vector(string &input,vector <double> &output)
{
	istringstream www(input);
	output.clear();
	int count=0;
	double value;
	for(;;)
	{
		if(! (www>>value) )break;
		output.push_back(value);
		count++;
	}
	return count;
}

//=========== load NNreg predicted abundance value ============//
int Load_NNreg_Abudance(string &input_file, vector <double> &output)
{
	ifstream fin;
	string buf,temp;
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	output.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		double value=atof(buf.c_str());
		output.push_back(value);
		count++;
	}
	//return
	return count;
}



//-------- determine abundance label -------//-> modified on 2016_10_20
int abundance_labels[15] = {1, 8, 15, 22, 33, 47, 68, 94, 130, 179, 252, 358, 524, 859, 1500};     // min_transcript_expression = 1.0
double abundance_labels_log[15]={0, 2.08, 2.71, 3.09, 3.50, 3.85, 4.22, 4.54, 4.87, 5.19, 5.53, 5.88, 6.26, 6.76, 7.31};
int locate_label(int abd)
{
	for(int k = 0; k < 14; k++)
	{
		if(abd < abundance_labels[k]) return k;
	}
	return 14;
}
int locate_label_log(double abd)
{
	for(int k = 0; k < 14; k++)
	{
		if(abd < abundance_labels_log[k]) return k;
	}
	return 14;
}

//================= load START/END prediction file ===============//
//-> format
/*
#-> 225
 0 -> 0.745086 0.019292 0.006568 0.005485 0.004362 0.007578 0.004736 0.006087 0.015438 0.012847 0.017210 0.018575 0.027046 0.041996 0.067693 -> 0.745086  0
 0 -> 0.862038 0.009437 0.002133 0.001682 0.001282 0.002516 0.001402 0.001866 0.006455 0.004903 0.007467 0.008011 0.013877 0.025433 0.051499 -> 0.862038  0
 0 -> 0.913919 0.005420 0.001236 0.000986 0.000752 0.001380 0.000791 0.001041 0.003440 0.002595 0.003975 0.004266 0.007770 0.015501 0.036926 -> 0.913919  0
 0 -> 0.939699 0.003742 0.000941 0.000753 0.000571 0.001004 0.000588 0.000771 0.002346 0.001809 0.002678 0.002884 0.005115 0.010315 0.026784 -> 0.939699  0
...
*/

//---- read between "->" ------//
void Read_Between_Arrow(string &input, string &output)
{
	int i;
	int len=(int)input.length();
	int start,end;
	//start
	start=-1;
	for(i=0;i<len-1;i++)
	{
		if(input[i]=='-' && input[i+1]=='>') //first arrow
		{
			start=i+3;
			break;
		}
	}
	//end
	end=-1;
	for(i=len-1;i>=1;i--)
	{
		if(input[i]=='>' && input[i-1]=='-') //second arrow
		{
			end=i-3;
			break;
		}
	}
	//judge
	if(start==-1 || end==-1)
	{
		fprintf(stderr,"bad start %d or end %d \n",start,end);
		exit(-1);
	}
	//final
	output=input.substr(start,end-start+1);
}

int Load_Pred_File(string &input_file,vector <string> &feat_rec)
{
	ifstream fin;
	string buf,temp;
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	int first=1;
	int feat_num_ori=-1;
	int length_s;
	for(;;)
	{
		//-> construct a new sequence
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
			exit(-1);
		}
		if(buf=="")continue;
		if(buf[0]=='#')
		{
			istringstream www(buf);
			string temp;
			www>>temp>>temp;
			length_s=atoi(temp.c_str());
			break;
		}
	}
	//-> read in features
	feat_rec.clear();
	for(int j=0;j<length_s;j++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
			exit(-1);
		}
		//--> check feature number
		int feat_num=Parse_FeatNum(buf);
		if(first==1)
		{
			first=0;
			feat_num_ori=feat_num;
		}
		else
		{
			if(feat_num!=feat_num_ori)
			{
				fprintf(stderr,"current feat_num %d not equal to first feat_num %d at line %s \n",
					feat_num,feat_num_ori,buf.c_str());
				exit(-1);
			}
		}
		//--> process
		string temp;
		Read_Between_Arrow(buf,temp);
		feat_rec.push_back(temp);
	}
	//return
	return (int)feat_rec.size();
}

//============== estimate abundance_value using predicted probability ===========//
void Estimate_Abundance_Value(vector < vector < double > > &abun_prob,vector <double> &nnreg_value)
{
	//-> estimate distance
	int length=(int)abun_prob.size();
	for(int i=0;i<length;i++)
	{
		//check label_0
		if(abun_prob[i][0]>0.8)
		{
			printf("0\n");
			continue;
		}
	
		//output estimated distance
		double estimated_abun=0;
		for(int k=0;k<15;k++)
		{
			double cur_abun=abundance_labels_log[k];
			estimated_abun+=cur_abun*abun_prob[i][k];
		}
		//output estimated variance
		double estimated_vari=0;
		for(int k=0;k<15;k++)
		{
			double cur_abun=abundance_labels_log[k];
			estimated_vari+=(cur_abun-estimated_abun)*(cur_abun-estimated_abun)*abun_prob[i][k];
		}
		double estimated_std=sqrt(estimated_vari);

		//---- upper and lower estimated variance -----//
		int end_lab=locate_label_log(estimated_abun);
		if(end_lab<0)end_lab=0;
		if(end_lab>13)end_lab=13;
		double end_dist=(estimated_abun-locate_label_log(end_lab))/(locate_label_log(end_lab+1) - locate_label_log(end_lab));
		if(end_dist<0)end_dist=0;
		if(end_dist>1)end_dist=1;
		//-> lower
		double estimated_vari1=0;
		double inner_vari1=0;
		for(int k=0;k<end_lab;k++)
		{
			double cur_abun=abundance_labels_log[k];
			estimated_vari1+=(cur_abun-estimated_abun)*(cur_abun-estimated_abun)*abun_prob[i][k];
		}
		{
			double cur_abun=abundance_labels_log[end_lab];
			estimated_vari1+=(cur_abun-estimated_abun)*(cur_abun-estimated_abun)*end_dist*abun_prob[i][end_lab];
		}
		double estimated_std1=sqrt(estimated_vari1);
		//-> upper
		double estimated_vari2=0;
		double inner_vari2=0;
		{
			double cur_abun=abundance_labels_log[end_lab];
			estimated_vari2+=(cur_abun-estimated_abun)*(cur_abun-estimated_abun)*(1-end_dist)*abun_prob[i][end_lab];
		}
		for(int k=end_lab+1;k<15;k++)
		{
			double cur_abun=abundance_labels_log[k];
			estimated_vari2+=(cur_abun-estimated_abun)*(cur_abun-estimated_abun)*abun_prob[i][k];
		}
		double estimated_std2=sqrt(estimated_vari2);

		//check label_0
		if(abun_prob[i][14]>0.45)
		{
			double lower_bound=estimated_abun-2.0*estimated_std1;
			if(lower_bound<0)lower_bound=0;
			if(nnreg_value[i]>lower_bound)printf("%lf\n",nnreg_value[i]);
			else printf("%lf\n",lower_bound);
			continue;
		}

		//check other labels
		{
			double lower_bound=estimated_abun-2.0*estimated_std1;
			double upper_bound=estimated_abun+2.0*estimated_std2;
			if(lower_bound<0)lower_bound=0;
			if(nnreg_value[i]>lower_bound && nnreg_value[i]<upper_bound)printf("%lf\n",nnreg_value[i]);
			else
			{
				if(nnreg_value[i]>upper_bound)printf("%lf\n",upper_bound);
				else printf("%lf\n",lower_bound);
			}
		}

	}
}


//----------- main -------------//
int main(int argc,char **argv)
{
	//---- Abundance_Estimate ----//__170101__//
	{
		if(argc<3)
		{
			fprintf(stderr,"Abundance_Estimate <abundance_prob> <nnreg_value> \n");
			exit(-1);
		}
		//argument
		string abundance_prob=argv[1];
		string nnreg_value=argv[2];
		//load feature
		vector <string> raw_feat;
		int len=Load_Pred_File(abundance_prob,raw_feat);
		vector <vector <double> > abun_prob;
		for(int i=0;i<len;i++)
		{
			vector <double> temp;
			String_To_Vector(raw_feat[i],temp);
			abun_prob.push_back(temp);
		}
		//load nnreg
		vector <double> nnreg_val;
		int ll=Load_NNreg_Abudance(nnreg_value, nnreg_val);
		//check length
		if(len!=ll)
		{
			fprintf(stderr,"len %d not equal to ll %d \n",len,ll);
			exit(-1);
		}
		//abundance estimate
		Estimate_Abundance_Value(abun_prob,nnreg_val);
		exit(0);
	}
}
