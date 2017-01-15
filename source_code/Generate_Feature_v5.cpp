#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <omp.h>
using namespace std;

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

//---- parse string -----//
int Parse_Str(string &in,vector <int> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		int value=atoi(buf.c_str());
		out.push_back(value);
		count++;
	}
	return count;
}
int Parse_Str_Double(string &in,vector <double> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		double value=atof(buf.c_str());
		out.push_back(value);
		count++;
	}
	return count;
}
int Parse_Str_String(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		out.push_back(buf);
		count++;
	}
	return count;
}


//----- calculate mean, vari ------//
void Calculate_Mean_Vari(vector <double> &input, double &mean, double &vari)
{
	int i;
	int size=(int)input.size();
	mean=0;
	vari=1;
	if(size==0)return;
	for(i=0;i<size;i++)
	{
		mean+=input[i];
	}
	mean/=size;
	vari=0;
	for(i=0;i<size;i++)
	{
		vari+=(input[i]-mean)*(input[i]-mean);
	}
	vari=sqrt(1.0*vari/size);
}


//-------- determine abundance label -------//-> modified on 2016_10_20
//int abundance_labels[12] = {15, 22, 33, 48, 68, 95, 132, 181, 255, 362, 530, 932};  //-> old bin
int abundance_labels[14] = {1, 8, 15, 22, 33, 47, 68, 94, 130, 179, 252, 358, 524, 859};     // min_transcript_expression = 1.0
int locate_label(int abd)
{
	for(int k = 0; k < 14; k++)
	{
		if(abd < abundance_labels[k]) return k;
	}
	return 14;
}

//--------- ATCG label -------//-> modified on 2016_10_23
void ATCG_Label(char c,string &out_str)
{
	out_str="";
	if(c=='A')out_str="1 0 0 0 0 ";
	else if(c=='T')out_str="0 1 0 0 0 ";
	else if(c=='C')out_str="0 0 1 0 0 ";
	else if(c=='G')out_str="0 0 0 1 0 ";
	else out_str="0 0 0 0 1 ";
}


//------- read Mingfu_Data ---------//-> modified on 2016_10_20
/*
As we discussed, one line (true abundance, under the label line) has been added to sample file.
I also changed labels from {-1, 0, 1} to {0, 1, 2} to make everything consistent.
*/
//-> data format
/*
# sample-id = 4, length = 636, location = chr1:4342282-4342918
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 
68 67 68 69 69 70 70 69 70 72 72 72 75 75 73 73 72 72 74 74 75 75 75 76 75 74 76 76 75 75 75 75 73 74 74 74 74 72 73 72 70 70 71 69 69 67 70 70 71 72 72 72 77 76 76 75 79 79 80 80 80 77 77 77 77 80 82 84 84 86 88 88 86 86 86 86 86 87 85 85 84 85 87 87 87 86 85 87 90 88 87 90 91 91 91 90 90 91 89 90 91 92 91 90 92 90 90 92 91 90 92 93 90 89 90 92 94 93 91 91 90 90 91 90 90 91 88 88 88 90 90 89 92 92 92 92 93 93 93 93 100 101 100 100 99 99 97 97 96 94 95 96 93 93 94 94 90 91 90 90 91 93 92 93 93 91 89 88 88 86 84 84 84 83 83 84 84 83 85 85 85 85 86 88 87 88 89 87 84 84 83 82 81 80 81 83 82 83 85 85 84 85 84 83 81 82 82 81 81 82 81 81 81 83 83 83 82 84 85 85 85 86 85 87 87 87 87 87 88 87 89 89 86 85 85 84 83 85 84 84 80 80 82 82 82 83 83 84 84 84 83 83 83 82 83 83 83 82 82 82 83 81 81 82 83 83 83 83 83 84 85 85 85 87 87 86 86 88 86 86 85 84 85 85 85 85 86 87 88 88 89 88 88 89 88 86 88 88 86 85 85 84 84 84 86 86 86 86 90 89 88 87 87 86 87 88 87 86 86 86 88 87 87 85 85 84 84 85 84 83 81 82 82 83 82 82 83 81 83 82 81 80 79 79 79 78 77 76 77 77 76 76 74 74 76 76 77 76 76 78 76 77 79 78 77 78 78 77 78 76 75 75 75 74 74 73 74 73 74 75 76 76 73 72 72 74 73 73 73 74 75 74 75 74 75 75 75 75 75 76 76 75 75 77 76 76 76 75 71 72 73 74 75 74 72 69 69 69 68 68 66 66 66 67 70 73 73 72 73 74 74 75 76 75 76 76 75 76 75 77 78 78 77 78 78 78 79 79 77 76 76 77 78 80 76 77 76 77 79 77 77 78 77 77 76 74 74 74 73 74 76 76 77 76 76 79 79 80 80 80 79 81 82 82 83 80 79 79 78 79 78 79 80 81 81 82 80 79 79 77 77 77 77 75 74 73 75 75 75 74 75 74 75 76 77 77 78 79 80 82 82 83 83 84 82 82 82 82 83 86 85 84 83 83 83 85 87 86 86 84 84 85 87 88 88 88 89 92 92 92 92 90 91 91 91 90 90 90 89 90 91 91 90 90 91 92 92 93 94 93 91 93 93 95 96 97 95 94 95 95 95 94 93 93 93 94 95 94 95 93 92 92 92 91 90 90 91 90 90 94 95 95 96 97 97 97 97 98 98 97 96 96 95 94 95 95 94 92 91 90 90 90 91 89 90 88 89 91 89 86 86 85 85 85 84 83
25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 32 32 32 35 35 35 35 35 35 35 39 39 39
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 22 22 22 20 20 20 20 20 20 20 18 18 18 18 18 18 20 20
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00
-1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00

*/
void Read_Single(vector <string> &input, 
	vector <int> &label_boundary, 
	vector <int> &label_abundance,
	vector <vector <double> > &feature_)
{
	Parse_Str(input[1],label_boundary,' ');
	Parse_Str(input[2],label_abundance,' ');
	//check length
	if(label_boundary.size() != label_abundance.size() )
	{
		fprintf(stderr,"label_boudanry size %d not equal to label_abundance size %d \n",
			label_boundary.size(),label_abundance.size());
		exit(-1);
	}
	//process
	int canon_len=(int)label_boundary.size();
	vector <vector <double> > feature;
	feature.clear();
	int dim=0;
	for(int i=3;i<(int)input.size();i++)
	{
		vector <double> tmp_feat;
		Parse_Str_Double(input[i],tmp_feat,' ');
		//check length
		if(tmp_feat.size()!=canon_len)
		{
			fprintf(stderr,"tmp_feat_size %d not equal to canon_len %d \n",
				tmp_feat.size(),canon_len);
			fprintf(stderr,"%s\n",input[0].c_str());
			exit(-1);
		}
		feature.push_back(tmp_feat);
		dim++;
	}
	//just perform log for the main feature
	for(int k=0;k<(int)feature[dim-1].size();k++)
	{
		feature[0][k]=log(1.0*feature[0][k]+1);
	}
	//final assign
	feature_.resize(canon_len);
	for(int i=0;i<canon_len;i++)
	{
		feature_[i].resize(dim);
		for(int k=0;k<dim;k++)feature_[i][k]=feature[k][i];
	}

/*
	//process final ATCG string
	{
		vector <string> tmp_feat;
		Parse_Str_String(input[(int)input.size()-1],tmp_feat,' ');
		//check length
		if(tmp_feat.size()!=canon_len)
		{
			fprintf(stderr,"tmp_feat_size %d not equal to canon_len %d \n",
				tmp_feat.size(),canon_len);
			fprintf(stderr,"%s\n",input[0].c_str());
			exit(-1);
		}
		//--- add one-hot vector ---//
		for(int i=0;i<canon_len;i++)
		{
			string out_str;
			ATCG_Label(tmp_feat[i][0],out_str);
			vector <double> tmp_feat;
			Parse_Str_Double(out_str,tmp_feat,' ');
			for(int j=0;j<(int)tmp_feat.size();j++)feature_[i].push_back(tmp_feat[j]);
		}
	}
*/
}

void Read_All(string &in_file,
	vector <vector <int> > &label_tot, 
	vector <vector <int> > &label_abu,
	vector <vector <vector <double> > > &feature_tot, int feat_dim)
{
	ifstream fin;
	string buf;
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input file not found [%s] !!!\n",in_file.c_str());
		exit(-1);
	}
	//read
	vector <string> tmp_rec;
	label_tot.clear();
	label_abu.clear();
	feature_tot.clear();
	int count=0;
	for(;;)
	{
		tmp_rec.clear();
		for(int i=0;i<feat_dim+3;i++)  //-> now we have 23 lines (20 features, and the remaining 3 lines are '#', bou_label, abu_label)
		{
			if(!getline(fin,buf,'\n'))goto end;
			tmp_rec.push_back(buf);
		}
		//check and process
		if(tmp_rec[0][0]!='#')
		{
			fprintf(stderr,"bad format ! \n");
			exit(-1);
		}
		vector <int> label_1;
		vector <int> label_2;
		vector <vector <double> > feature;
		Read_Single(tmp_rec,label_1,label_2,feature);
		label_tot.push_back(label_1);
		label_abu.push_back(label_2);
		feature_tot.push_back(feature);
		//count
		count++;
		printf("cur=%d\r",count);
	}
end:
	;
}

//------ output feature --------//
void Output_Feature(vector <vector <int> > &label_tot, vector <vector <double> > &feature_tot, 
	int label, int window_size,double samp_ratio,FILE *fp)
{
	int i,j,k;
	for(i=0;i<(int)label_tot.size();i++)
	{
		for(j=0;j<(int)label_tot[i].size();j++)
		{
			if(label_tot[i][j]==label) //target label
			{
				//random number
				double rand_num=1.0* rand() / RAND_MAX;
				if(rand_num<0)rand_num=0;
				else if(rand_num>1)rand_num=1;
				if(rand_num>samp_ratio)continue;
				//get window size
				vector <double> value;
				value.clear();
				for(k=j-window_size;k<=j+window_size;k++)
				{
					if(k<0)
					{
						for(int w=0;w<3;w++)value.push_back(0);
					}
					else if(k>=(int)label_tot[i].size())
					{
						for(int w=0;w<3;w++)value.push_back(0);
					}
					else
					{
						//center
						value.push_back(feature_tot[i][k]);
						//left
						if(k>0)value.push_back(feature_tot[i][k-1]-feature_tot[i][k]);
						else value.push_back(0);
						//right
						if(k<(int)label_tot[i].size()-1)value.push_back(feature_tot[i][k+1]-feature_tot[i][k]);
						else value.push_back(0);
					}
				}
				//output
				fprintf(fp,"1|");
				for(k=0;k<(int)value.size();k++)fprintf(fp,"%5.3f ",value[k]);
				fprintf(fp,"|");
				fprintf(fp,"%d\n",label);
			}
		}
	}
}

//----------- main -------------//
int main(int argc,char **argv)
{
	//------- Generate_Feature --------// 
	{
		if(argc<8)
		{
			fprintf(stderr,"Generate_Feature_v5 <raw_input> <root_feat> <root_cov> <root_l1> <root_l2> <root_l3> <feat_dim> \n");
			fprintf(stderr,"[note1]: root_l1 is for 'End Mid Start' (0,1,2) labels \n");
			fprintf(stderr,"         root_l2 is for abundance original value \n");
			fprintf(stderr,"         root_l3 is for 15 abundance labels \n"); 
			fprintf(stderr,"[note2]: feat_dim is the feature dimention. [default=20] \n");
			exit(-1);
		}
		string raw_input=argv[1];
		string feat_root=argv[2];
		string cov_root=argv[3];
		string lab_root_1=argv[4];
		string lab_root_2=argv[5];
		string lab_root_3=argv[6];
		int feat_dim=atoi(argv[7]);
		//mkdir -p
		char command_[30000];
		sprintf(command_,"mkdir -p %s",feat_root.c_str());
		system(command_);
		sprintf(command_,"mkdir -p %s",cov_root.c_str());
		system(command_);
		sprintf(command_,"mkdir -p %s",lab_root_1.c_str());
		system(command_);
		sprintf(command_,"mkdir -p %s",lab_root_2.c_str());
		system(command_);
		sprintf(command_,"mkdir -p %s",lab_root_3.c_str());
		system(command_);
		//process
		vector <vector <int> > label_tot;
		vector <vector <int> > label_abu;
		vector <vector <vector <double> > > feature_tot;
		Read_All(raw_input,label_tot,label_abu,feature_tot,feat_dim);
		//output
		#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<(int)label_tot.size();i++)
		{
			FILE *fp;
			char command[30000];
			//-> features
			sprintf(command,"%s/sample_%d",feat_root.c_str(),i);
			fp=fopen(command,"wb");
			fprintf(fp,"%d\n",(int)label_tot[i].size());
			for(int k=0;k<(int)label_tot[i].size();k++)
			{
				int featdim=(int)feature_tot[i][k].size();
				stringstream oss;
				for(int l=0;l<featdim;l++)
				{
					int wsiii=(int)feature_tot[i][k][l];
					if(wsiii!=feature_tot[i][k][l])oss << feature_tot[i][k][l] << " ";
					else oss << wsiii << " ";
				}
				string wsbuf=oss.str();
				fprintf(fp,"%s\n",wsbuf.c_str());
			}
			fclose(fp);
			//-> coverage
			sprintf(command,"%s/sample_%d",cov_root.c_str(),i);
			fp=fopen(command,"wb");
			for(int k=0;k<(int)label_tot[i].size();k++)fprintf(fp,"%lf\n",feature_tot[i][k][0]);
			fclose(fp);
			//-> labels (boudanry and abudance)
			//--> boundary label
			sprintf(command,"%s/sample_%d",lab_root_1.c_str(),i);
			fp=fopen(command,"wb");
			for(int k=0;k<(int)label_tot[i].size();k++)fprintf(fp,"%d\n",label_tot[i][k]);
			fclose(fp);
			//--> abundance real value
			sprintf(command,"%s/sample_%d",lab_root_2.c_str(),i);
			fp=fopen(command,"wb");
			for(int k=0;k<(int)label_abu[i].size();k++)fprintf(fp,"%d\n",label_abu[i][k]);
			fclose(fp);
			//--> abundance 15 label
			sprintf(command,"%s/sample_%d",lab_root_3.c_str(),i);
			fp=fopen(command,"wb");
			for(int k=0;k<(int)label_abu[i].size();k++)
			{
				int label=locate_label(label_abu[i][k]);
				fprintf(fp,"%d\n",label);
			}
			fclose(fp);
		}
		exit(0);
	}
}
