#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
using namespace std;


//----- just utility ----//
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}

//================ WS_Basic_Script for List_Differ ==============//__
//[1] given list cur + nxt, return the difference list that only appeared in nxt
void WS_Return_Difference(string &list1,string &list2)
{
	map<string, int > ws_mapping;   //M, mapping the PDB's name
	ws_mapping.clear();
	ifstream fin;
	string buf;
	//read list1 for basic mapping
	fin.open(list1.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list1 %s not found!!\n",list1.c_str());
		exit(-1);
	}
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		count++;
		ws_mapping.insert(map < string, int >::value_type(buf, count));
	}
	fin.close();
	fin.clear();

	//read list2 for checking
	fin.open(list2.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list2 %s not found!!\n",list2.c_str());
		exit(-1);
	}
	int key;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		key=ws_mapping[buf];
		if(key<1 || key>count)printf("%s\n",buf.c_str());
	}
	fin.close();
	fin.clear();
}
//[2] given list1 and list2, return the intersection
void WS_Return_Intersection(string &list1,string &list2)
{
	map<string, int > ws_mapping;   //M, mapping the PDB's name
	ws_mapping.clear();
	ifstream fin;
	string buf;
	//read list1 for basic mapping
	fin.open(list1.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list1 %s not found!!\n",list1.c_str());
		exit(-1);
	}
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		count++;
		ws_mapping.insert(map < string, int >::value_type(buf, count));
	}
	fin.close();
	fin.clear();

	//read list2 for checking
	fin.open(list2.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list2 %s not found!!\n",list2.c_str());
		exit(-1);
	}
	int key;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		key=ws_mapping[buf];
		if(key<1 || key>count)continue;
		else printf("%s\n",buf.c_str());
	}
	fin.close();
	fin.clear();
}

//[2] given list1 and list2, return the union
void WS_Return_Union(string &list1,string &list2)
{
	map<string, int > ws_mapping;   //M, mapping the PDB's name
	ws_mapping.clear();
	ifstream fin;
	string buf;
	//read list1 for basic mapping
	fin.open(list1.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list1 %s not found!!\n",list1.c_str());
		exit(-1);
	}
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		count++;
		ws_mapping.insert(map < string, int >::value_type(buf, count));
		printf("%s\n",buf.c_str());
	}
	fin.close();
	fin.clear();

	//read list2 for checking
	fin.open(list2.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list2 %s not found!!\n",list2.c_str());
		exit(-1);
	}
	int key;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		key=ws_mapping[buf];
		if(key<1 || key>count)printf("%s\n",buf.c_str());
	}
	fin.close();
	fin.clear();
}

//------- ws_fast_cat ---------//__120430__//
int WS_Read_In_Cat(string &list,vector <string> &rec)
{
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",list.c_str());
		return -1;
	}
	rec.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		rec.push_back(buf);
	}
	return (int)rec.size();
}
void WS_Fast_Cat(string &list,string &out)
{
	//read in list
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",list.c_str());
		exit(-1);
	}
	vector <string> ws_list_rec;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		ws_list_rec.push_back(buf);
	}
	//fast cat
	int i;
	int totnum=(int)ws_list_rec.size();
	vector <vector <string> > ws_total_rec;
	ws_total_rec.resize(totnum);
	//the following could be run under OpenMP
	#pragma omp parallel for schedule(dynamic)
	for(i=0;i<totnum;i++)
	{
		WS_Read_In_Cat(ws_list_rec[i],ws_total_rec[i]);
	}
	//======= final output =====//
	FILE *fp=fopen(out.c_str(),"wb");
	int j;
	for(i=0;i<totnum;i++)for(j=0;j<(int)ws_total_rec[i].size();j++)fprintf(fp,"%s\n",ws_total_rec[i][j].c_str());
	fclose(fp);
}





//-------------- main ---------------//
int main(int argc,char **argv)
{
	//--------- WS_FAST_Cat -----------//_V2_//
	{
		if(argc<3)
		{
			printf("Fast_CAT <list> <out> \n");
			exit(-1);
		}
		string list=argv[1];
		string out=argv[2];
		WS_Fast_Cat(list,out);
		exit(0);
	}
}
