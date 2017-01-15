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
		//-> break
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

//------- read in real_value label ---------//
int Parse_Str_Double(string &in,vector <double> &out, char separator)
{
	istringstream www(in);
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
void Double_To_String(vector <double> &features, string &out)
{
	int i;
	int size=(int)features.size();
	stringstream oss;
	for(i=0;i<size;i++)
	{
		int wsiii=(int)features[i];
		if(wsiii!=features[i])oss << features[i] << " ";
		else oss << wsiii << " ";
	}
	out=oss.str();
}

//---------- load label file -------//
int Load_Label_File(string &input_file,vector <double> &label_rec)
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
		double label=atof(buf.c_str());
		int wsiii=(int)label;
		if(wsiii!=label)label_rec.push_back(label);
		else label_rec.push_back(wsiii);
	}
	//return
	return (int)label_rec.size();
}


//--------- calculate average value of a vector ---------//
void Average_Value_Vector_Single(vector < vector <double> > &input, int start,int end,string &out_str)
{
	int i,j;
	int len=(int)input[0].size();
	vector <double> output;
	output.resize(len,0);
	for(i=start;i<=end;i++)
	{
		for(j=0;j<len;j++)output[j]+=input[i][j];
	}
	int size=end-start+1;
	for(j=0;j<len;j++)output[j]/=size;
	//transfer to string
	Double_To_String(output, out_str);
}
void Average_Value_Vector(vector < vector <double> > &feat_mat, vector <string> &output, int half_wind)
{
	int i;
	int size=(int)feat_mat.size();
	output.clear();
	for(i=0;i<size;i++)
	{
		//determine start and end
		int start,end;
		if(half_wind<0)      //left
		{
			start=i+half_wind;
			end=i;
			if(start<0)start=0;
		}
		else if(half_wind>0) //right
		{
			start=i;
			end=i+half_wind;
			if(end>=size)end=size-1;
		}
		else                 //center
		{
			start=i;
			end=i;
		}
		//process
		string str_out;
		Average_Value_Vector_Single(feat_mat, start,end,str_out);
		//push_back
		output.push_back(str_out);
	}
}


//--------- main process ----------//
//-> we consider the average values in different half-windows
/*
   -50   -25    -10    -5   0    +5    +10    +25      +50
*/
void Main_Process(string &pred_file1, string &pred_file2, string &cov_file, string &lab_file)
{
	//load pred1 (START,MID,END)
	vector <string> feat_rec_1;
	int pred_len1=Load_Pred_File(pred_file1,feat_rec_1);
	//load pred2 (abundance label)
	vector <string> feat_rec_2;
	int pred_len2=Load_Pred_File(pred_file2,feat_rec_2);
	//load coverage
	vector <double> feat;
	int cov_len=Load_Label_File(cov_file,feat);
	//load label
	vector <double> label;
	int lab_len=Load_Label_File(lab_file,label);
	//check length
	if(pred_len1!=pred_len2 || pred_len1!=cov_len || pred_len1!=lab_len )
	{
		fprintf(stderr,"pred_len1 %d not equal to pred_len2 %d or cov_len %d or lab_len %d \n",
			pred_len1,pred_len2,cov_len,lab_len);
		exit(-1);
	}
	//--- create feature matrix ---//
	int first=1;
	int rec_len;
	vector < vector <double> > feat_mat;
	feat_mat.clear();
	for(int i=0;i<pred_len1;i++)
	{
		//feat_vec
		vector <double> feat_vec;
		feat_vec.clear();
		//l1
		int l1=Parse_Str_Double(feat_rec_1[i],feat_vec,' ');
		//l2
		int l2=Parse_Str_Double(feat_rec_2[i],feat_vec,' ');
		//l3
		feat_vec.push_back(feat[i]);
		//fin_len
		int fin_len=l1+l2+1;
		//check
		if(first==1)
		{
			first=0;
			rec_len=fin_len;
		}
		else
		{
			if(fin_len!=rec_len)
			{
				fprintf(stderr,"fin_len %d not equal to rec_len %d \n",fin_len,rec_len);
				exit(-1);
			}
		}
		//push_back
		feat_mat.push_back(feat_vec);
	}
	//----- generate feature string with windowns ------//
	vector < vector <string> > feat_str_tot;
	vector <string> feat_str_cur;
	feat_str_tot.clear();
	//-> add half_wind = -50
	Average_Value_Vector(feat_mat, feat_str_cur, -50);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = -25
	Average_Value_Vector(feat_mat, feat_str_cur, -25);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = -10
	Average_Value_Vector(feat_mat, feat_str_cur, -10);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = -5
	Average_Value_Vector(feat_mat, feat_str_cur, -5);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = 0
	Average_Value_Vector(feat_mat, feat_str_cur, 0);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = +5
	Average_Value_Vector(feat_mat, feat_str_cur,  5);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = +10
	Average_Value_Vector(feat_mat, feat_str_cur,  10);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = +25
	Average_Value_Vector(feat_mat, feat_str_cur,  25);
	feat_str_tot.push_back(feat_str_cur);
	//-> add half_wind = +50
	Average_Value_Vector(feat_mat, feat_str_cur,  50);
	feat_str_tot.push_back(feat_str_cur);

	//output
	for(int i=0;i<pred_len1;i++)
	{
		printf("%lf ",log(label[i]+1));
		for(int j=0;j<(int)feat_str_tot.size();j++)printf("%s ",feat_str_tot[j][i].c_str());
		printf("\n");
	}
}


//----------- main -------------//
int main(int argc,char **argv)
{
	//---- Mingfu_Feat_Merge_wind ----//
	{
		if(argc<5)
		{
			fprintf(stderr,"Mingfu_Feat_Merge_wind <pred_file_bou> <pred_file_abu> <cov_file> <label_file> \n");
			fprintf(stderr,"[note]: pred_file_bou should be (START,MID,END) boundary prediction \n");
			fprintf(stderr,"        pred_file_abu should be (abundance label) prediction \n");
			fprintf(stderr,"        cov_file should be coverage file\n");
			fprintf(stderr,"        label_file should contain ground_truth abudance \n");
			fprintf(stderr,"[wind]: the output will use -50,-25,-10,-5,0,+5,+10,+25,+50 hafl_wind \n");
			exit(-1);
		}
		//argument
		string pred_file_1=argv[1];
		string pred_file_2=argv[2];
		string cov_file=argv[3];
		string label_file=argv[4];
		//process
		Main_Process(pred_file_1,pred_file_2,cov_file,label_file);
		exit(0);		
	}
}
