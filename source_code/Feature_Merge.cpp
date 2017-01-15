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

//================= load feature file ===========//
//-> format
/*
105
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
4.94876 14.2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
...
*/

int Load_Feat_File(string &input_file,vector <string> &feat_rec)
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
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		int length_s=atoi(buf.c_str());
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
			feat_rec.push_back(buf);
		}
		break;
	}
	//return
	return (int)feat_rec.size();
}

//================= load START/END prediction file ===============//
//-> format
/*
#-> 105
13 -> 0.169958 0.719566 0.110477 -> 0.719566  1
13 -> 0.171116 0.472000 0.356884 -> 0.472000  1
13 -> 0.250478 0.529205 0.220316 -> 0.529205  1
13 -> 0.174095 0.690448 0.135458 -> 0.690448  1
13 -> 0.249481 0.488221 0.262297 -> 0.488221  1
13 -> 0.199086 0.669642 0.131272 -> 0.669642  1
13 -> 0.103673 0.817112 0.079215 -> 0.817112  1
13 -> 0.091546 0.827624 0.080830 -> 0.827624  1
...
*/

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
		istringstream www(buf);
		string temp;
		string f1,f2,f3;
		www>>temp>>temp>>f1>>f2>>f3;
		temp=f1+" "+f2+" "+f3;
		feat_rec.push_back(temp);
	}
	//return
	return (int)feat_rec.size();
}

//---------- load label file -------//
int Load_Label_File(string &input_file,vector <int> &label_rec)
{
	ifstream fin;
	string buf,temp;
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//process
	label_rec.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		int label=atoi(buf.c_str());
		label_rec.push_back(label);
	}
	//return
	return (int)label_rec.size();
}


//--------- main process ----------//
void Main_Process(string &feat_file, string &pred_file, string &label_file)
{
	//load feat
	vector <string> feat_rec;
	int feat_len=Load_Feat_File(feat_file,feat_rec);
	//load pred
	vector <string> pred_rec;
	int pred_len=Load_Pred_File(pred_file,pred_rec);
	//load label
	vector <int> label_rec;
	int label_len=Load_Label_File(label_file,label_rec);
	//check length
	if(feat_len!=pred_len || feat_len!=label_len)
	{
		fprintf(stderr,"feat_len %d not equal to pred_len %d or label_len %d\n",
			feat_len,pred_len,label_len);
		exit(-1);
	}
	//output
	printf("%d\n",feat_len);
	for(int i=0;i<feat_len;i++)
	{
		printf("%s %s\n",feat_rec[i].c_str(),pred_rec[i].c_str());
	}
	for(int i=0;i<feat_len;i++)
	{
		printf("%d\n",label_rec[i]);
	}
	
}


//----------- main -------------//
int main(int argc,char **argv)
{
	//---- Mingfu_Feat_Merge ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Mingfu_Feat_Merge <feat_file> <pred_file> <label_file> \n");
			exit(-1);
		}
		//argument
		string feat_file=argv[1];
		string pred_file=argv[2];
		string label_file=argv[3];
		//process
		Main_Process(feat_file,pred_file,label_file);
		exit(0);		
	}
}
