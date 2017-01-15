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

//================= decompose '#' commented file into sub files =============//
int Decompose_Commented_File(string &input_file, string &out_root, string &out_name)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//load
	int first=1;
	int count=0;
	char command[30000];
	vector <string> tmp_rec;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')
		{
			if(first==1)
			{
				first=0;
			}
			else
			{
				//output file
				sprintf(command,"%s/%s_%d",out_root.c_str(),out_name.c_str(),count);
				string out_file=command;
				FILE *fp=fopen(out_file.c_str(),"wb");
				for(int i=0;i<(int)tmp_rec.size();i++)fprintf(fp,"%s\n",tmp_rec[i].c_str());
				fclose(fp);
				//count++
				count++;
			}
			//clear
			tmp_rec.clear();
			tmp_rec.push_back(buf);
		}
		else
		{
			tmp_rec.push_back(buf);
		}
	}
	//final check
	if(tmp_rec.size()>0)
	{
		//output file
		sprintf(command,"%s/%s_%d",out_root.c_str(),out_name.c_str(),count);
		string out_file=command;
		FILE *fp=fopen(out_file.c_str(),"wb");
		for(int i=0;i<(int)tmp_rec.size();i++)fprintf(fp,"%s\n",tmp_rec[i].c_str());
		fclose(fp);
		//count++
		count++;
	}
	//return
	return count;
}

//----------- main -------------//
int main(int argc,char **argv)
{
	//---- Decompose_Commented_File ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Decompose_Commented_File <input_commented_file> <output_root> <output_name> \n");
			fprintf(stderr,"[note]:      output files would be put to output_root/output_name_X \n");
			exit(-1);
		}
		//argument
		string input_commented_file=argv[1];
		string output_root=argv[2];
		string output_name=argv[3];
		//process
		Decompose_Commented_File(input_commented_file,output_root,output_name);
		exit(0);		
	}
}
