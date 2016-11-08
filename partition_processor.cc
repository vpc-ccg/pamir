#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <set>
#include <string.h>
#include "common.h"

using namespace std;

int main(int argc, char **argv)
{
	int INSSIZE=450;
	//map token is from oea to list of orphan contigs that can be mapped
	map<string, vector<vector<string> > > map_token;
	map<string, int> map_sup;
	map<string, string> map_cont;

	//map<string, vector<string> >::iterator it;
	vector<string>::iterator it;	
	vector<vector<string> >::iterator nit;	
	
	ifstream fin;

	string cont;
	int sup;
	string rname,gname;

	// reading the conting content

	fin.open(argv[1]);
	while( fin >> gname >> cont )
	{
		rname = gname.substr(1, gname.length());
		map_sup[rname] = 0;
		map_cont[rname] = cont;
	}
	fin.close();
	cerr << "Loading orphan contig file"<<endl; 

	// orphan to orphan-contig mapping file; 
	fin.open(argv[2]);
	string  name, ref, cigar, nref, rstr, qstr, t1, t2;
	int flag, pos, qual, npos, tlen;  
	while (fin >> name >> flag >> ref >> pos >> qual >> cigar >> nref >> npos >> tlen >> rstr >> qstr >> t1 >> t2) {
		map_sup[ref]++;
	}

	fin.close();
	cerr << "Updating orphan contig support"<<endl; 

	// oea to orphan-contig mapping file
	fin.open(argv[3]);
	while (fin >> rname >> flag >> gname >> pos >> qual >> cigar >> nref >> npos >> tlen >> rstr >> qstr >> t1 >> t2)
	{
		if(rname[rname.length()-2]=='/'){
			rname.erase(rname.length()-2,2);
		}
		vector<string> tmpv;
		tmpv.push_back(gname);
		tmpv.push_back( ( (flag & 0x10) == 0x10 )?"-":"+" );
		
		int contiglen = map_cont[gname].length();
		int readlen   = rstr.length();
		if ( pos < INSSIZE)
			tmpv.push_back("l");
		else if(pos > contiglen-INSSIZE-readlen)
			tmpv.push_back("r");
		else
			tmpv.push_back("w");

		if (map_token.find(rname) != map_token.end() ) {
			map_token[rname].push_back(tmpv);
		}
		else
		{
			vector<vector<string> > tmp_vec;
			tmp_vec.push_back(tmpv);
			map_token[rname] = tmp_vec;
		}

	}
	fin.close();
	cerr << "Updating OEA contigs"<<endl; 


	// original oea partition file

	// converted oea mapped file

	fin.open(argv[4]);
	

	int lines, start, end, cluster_id, numOfReads;
	int i;
	vector <string> reads;
	vector <string> reads_content;
	vector <int> reads_support;
	vector <int> reads_pos;
	
	map<string, vector<int> > myset;
	map<string, vector<int> >::iterator mit;
	char *logfname = new char[1000];
	sprintf(logfname,"%s_partitionprocessor.log",argv[4]);

	FILE *fucluster = fopen(argv[5],"w");
	FILE *fidx	 = fopen((string(argv[5])+".idx").c_str(),"wb");
	FILE *flog 	 = fopen(logfname, "w");
	int tie 	 = 0;
	while (fin>> cluster_id >> lines >> start >> end >> ref){
		fprintf(flog,"CLUSTER ID: %d\n",cluster_id);
		numOfReads=lines;
		reads.clear();
		reads_content.clear();
		reads_support.clear();
		reads_pos.clear();

		for ( i=0; i<lines; i++) {
			fin >> rname >> rstr >> sup >> pos;
			reads.push_back(rname);
			reads_content.push_back(rstr);
			reads_support.push_back(sup);
			reads_pos.push_back(pos);
		}
		cerr << "Partition Loaded... "; 
		cerr << lines << " oea";
		myset.clear();
		char sign;
		for (it = reads.begin(); it != reads.end(); it++)
		{
			qstr = (*it).substr(0, (*it).size()-1);			
			sign = (*it)[(*it).size()-1];
			if (map_token.find(qstr) != map_token.end()) {
				for (nit=map_token[qstr].begin(); nit != map_token[qstr].end(); nit++){
					if (myset.find((*nit)[0]) == myset.end()) {
						vector<int> tmpev;
						tmpev.push_back(0);
						tmpev.push_back(0);
						tmpev.push_back(0);
						tmpev.push_back(0);
						myset[(*nit)[0]] = tmpev;
					}

					if (sign != (*nit)[1][0])
						myset[(*nit)[0]][1]++;
					else
						myset[(*nit)[0]][0]++;

					if ((*nit)[2] == "l")
						myset[(*nit)[0]][2]++;
					else
						myset[(*nit)[0]][3]++;
				}	
			}
		}
		cerr << "Processed..."; 
		int acceptedContigNum =0;
		for (mit=myset.begin(); mit!=myset.end(); mit++) 
		{
			fprintf(flog,"readNum: %d\tcontigname: %s\tleft: %d\tright: %d\tforward: %d\treverse: %d\n", numOfReads, (*mit).first.c_str(), myset[(*mit).first][2], myset[(*mit).first][3], myset[(*mit).first][0], myset[(*mit).first][1]);
			int decided =0;
			if (map_sup.find((*mit).first) != map_sup.end() ) 
			{
				decided = 1;
				int contiglen = map_cont[(*mit).first].length();
				if(contiglen <= INSSIZE)
				{
					if(myset[(*mit).first][2] ==0)
					{
						decided =0;
						fprintf(flog, "Short contig and no left or right flank supporters\n");
					}
				}
				else
				{
					if((myset[(*mit).first][2]==0|| myset[(*mit).first][3]==0) && myset[(*mit).first][2] +myset[(*mit).first][3] < 0.3*numOfReads )
					{
						decided =0;
						fprintf(flog, "either no OEAs mapping on left or right flank and not enough left + right flank support < 30%%\n");
					}
				}
				if ( myset[(*mit).first][0] == myset[(*mit).first][1])
				{
					tie++;
					fprintf(flog, "Forward and reverse are exactly same, we cannot decide which version of the contig to add\n");
				}

			}			
			if(decided == 1)
			{
				string revcontent;
				//if (strcmp( (*mit).first.c_str(),"contig-222")==0 && cluster_id ==2129)
				//	cerr<<"\nXXXXXXXXXXXXXXXXXXXX\t"<<myset[ (*mit).first][0]<<"\t"<<myset[(*mit).first][1]<<endl;
				if( myset[(*mit).first][0] ==  myset[(*mit).first][1] )
				{
				//	if (strcmp( (*mit).first.c_str(),"contig-222")==0)
				//		cerr<<"IF"<<endl;
					reads.push_back((*mit).first);
					reads_content.push_back(map_cont[(*mit).first]);
					reads_support.push_back(map_sup[(*mit).first]);
					reads_pos.push_back(-1);
					string content = map_cont[(*mit).first];
					revcontent = reverse_complement(content);
					reads.push_back((*mit).first + "_rv");
					reads_content.push_back(revcontent);
					reads_support.push_back(map_sup[(*mit).first]);
					reads_pos.push_back(-1);
					acceptedContigNum+=2;

				}
				else
				{
				//	if (strcmp( (*mit).first.c_str(),"contig-222")==0 && cluster_id == 2129)
				//		cerr<<"ELSE"<<endl;
					if(myset[(*mit).first][0] > myset[(*mit).first][1])
					{
						revcontent = map_cont[(*mit).first];
					}
					else if(myset[(*mit).first][1] > myset[(*mit).first][0])
					{
						string content = map_cont[(*mit).first];
						revcontent = reverse_complement(content);
					}
					reads.push_back((*mit).first);
					reads_content.push_back(revcontent);
					reads_support.push_back(map_sup[(*mit).first]);
					reads_pos.push_back(-1);
					acceptedContigNum++;
				}
			}
		}
		lines+=acceptedContigNum;
		
		size_t fout_size;
		fout_size= ftell(fucluster);
		fwrite(&fout_size,1,sizeof(size_t),fidx);
		fprintf(fucluster,"%d %d %d %d %s\n",cluster_id, lines, start, end, ref.c_str());

		for (i=0; i<lines; i++) {
			fprintf(fucluster,"%s %s %d %d\n",reads[i].c_str(),reads_content[i].c_str(),reads_support[i],reads_pos[i]);
		}

		if ( 0 < myset.size() )
		{	cerr << "Updated " <<cluster_id << "\t" << lines << "("<< myset.size()<<")." << endl;}
		else
		{   cerr << "Unchanged." << endl;}
	}
	fprintf(flog, "Forward and reverse were in tie condition : %d times.\n", tie );
	
	cout<<cluster_id<<endl;
	fclose(fucluster);
	fclose(fidx);
}
