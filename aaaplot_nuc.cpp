//aaaplot_nuc.cpp    Brian Peng He 2012 Dec
//This program plots the absolute density of reads near the boundaries of defined regions with one-nucleotide resolution
//strand considered, no curve smoothing used 
//screen output: how many lines of position file and region file are checked, which file reached the end first as well as the current progress
//
//
//input format must be:
//region file:                            pos file:('chr' is removed, X => 23, Y => 24)
//chrX + 198892 199382             ||     23 + 92738
//chrY - 129880 149790             ||     24 - 10480
//input files MUST be sorted 
// sort -k 1d,1 3n,3 RegionFile >> 
// sort -k 1n,1 3n,3 PosFile >>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

const int Gwidth = 100; // # of bins for gene region
#define Twidth (Gwidth*4)

using namespace std;  

  
int mystoi(string chr_c){   //function to convert string to integer
                int chr =0;
		if(chr_c=="chr1")chr=1;
		else if(chr_c=="chr2")chr=2;
		else if(chr_c=="chr3")chr=3;
		else if(chr_c=="chr4")chr=4;
		else if(chr_c=="chr5")chr=5;
		else if(chr_c=="chr6")chr=6;
		else if(chr_c=="chr7")chr=7;
		else if(chr_c=="chr8")chr=8;
		else if(chr_c=="chr9")chr=9;
		else if(chr_c=="chr10")chr=10;
		else if(chr_c=="chr11")chr=11;
		else if(chr_c=="chr12")chr=12;
		else if(chr_c=="chr13")chr=13;
		else if(chr_c=="chr14")chr=14;
		else if(chr_c=="chr15")chr=15;
		else if(chr_c=="chr16")chr=16;
		else if(chr_c=="chr17")chr=17;
		else if(chr_c=="chr18")chr=18;
		else if(chr_c=="chr19")chr=19;
		else if(chr_c=="chr20")chr=20;
		else if(chr_c=="chr21")chr=21;
		else if(chr_c=="chr22")chr=22;
		else if(chr_c=="chrX")chr=23;
		else if(chr_c=="chrY")chr=24;
		return chr;
			}


int coverage = 0;//a global variable for normalization

				
void iofile(string pf, string rf, unsigned long long int mCXX[Twidth]){//the core function
     //////////////////////////////////////////////////////////////////////////
     /////////                                                      ///////////
     /////////                                                      ///////////
     /////////          left    interesting      right              ///////////
     ///////// a-100---a-1|a---a+99-----b-100---b-1|b-------b+99     ///////////
     /////////                                                      ///////////
     //////////////////////////////////////////////////////////////////////////
     

    int c=0;   //chromosome int
    string C;// chromosome string
    char strand;
    int p;     //position in the POSITION file


    int chr,a,b,ab;    //chromosome,start,end in the REGION file. ab marks the boundary
    char str;  //strand in region file


    int clr = 0;//count the lines checked in region file
    int clp = 0;//count the lines checked in position file


	for(int round = 1;round < 3; round+=1){//two rounds: upstreamboundary ---- downstreamboundary
		coverage=0;
		clr=0; clp=0;

		ifstream pos;    //open the position file(mapped)
    		pos.open(pf.c_str());
    		if(pos.fail()){cout<<"position fail";exit(1);}
        
    		ifstream region;   //open the region file(from Andrey) 
    		region.open(rf.c_str());
    		if(region.fail()){cout<<"bed fail";exit(1);}


		    region>>chr>>str>>a>>b;
			coverage++;
			clr++;
    			pos>>C>>strand>>p;
    			clp++; 
			c = mystoi(C);  //initialize: get first lines of PositionFile and RegionFile		
		
		if(1==round){ab=a-Gwidth;}
		else if(2==round){ab=b-Gwidth;}
		else {cout<<"round error! impossible!";exit(1);}
		while(true){ 

	                if(c>chr||(c==chr&&p>=ab+Gwidth)){//if the position is on the right of the region
					if(region>>chr==false){cout<<"check2";break;}//reached the end of the region file 
					else	{
						region>>str>>a>>b;
						coverage++;
					if(1==round){ab=a-Gwidth;}
					else if(2==round){ab=b-Gwidth;}
					else {cout<<"round error! impossible!";exit(1);}
						clr++;
						}
					continue;
							}// go to the next region
			else if(c==chr&&ab<=p && p < ab+2*Gwidth){//if the position falls in
				unsigned int LOCTN = p-ab+(round-1)*2*Gwidth;      
					if(str=='+')mCXX[LOCTN]++;//increase the corresponding bin's signal strength 
        	                  	else if(str=='-')mCXX[Twidth-1-LOCTN]++;
        	                  	else    {cout<<"strand neither + nor -"<<endl;exit(1);}
									
								}
		if(pos>>C==false){cout<<"check1";break;}//get one more position
		pos>>strand>>p;
		clp++;
		c=mystoi(C);    
			}//end of while loop         	
    region.close();
    pos.close();
cout<<"went through "<<clr<<" r_lines "<<clp << " p_lines" << endl;                    
  						}//end of the second ronud
	}
	
int main(){
    cout<<"How many position files?(less than 100 plz)"<<endl;
    int np;
    string ofCXX;
    cin >> np;
    cout<<"show me the position files' names, please."<<endl;
    string pf[100];
    for(int i = 0;i<np;i++)cin>>pf[i];//get methylation file order

    cout<<"How many region files then?(less than 100 plz)"<<endl;
    int nr;
    cin >> nr;
    cout<<"And their names?"<<endl;
    string rf[100];
    for(int j = 0;j<nr;j++)cin>>rf[j];//get region files' list

    for(int j = 0;j<nr;j++){
	for(int k = 0;k<np;k++){
		unsigned long long int mCXX[Twidth];//define an array to store the signal
		for(int i=0;i<Twidth;i++)mCXX[i]=0;
		cout << k << '\t' << j << endl;		
		iofile(pf[k],rf[j], mCXX);//for each combination of region & file, plot 
		ofstream outfile;
		ofCXX=rf[j]+pf[k];
		outfile.open(ofCXX.c_str());
		for(int i=0;i<Twidth;i++) outfile <<i<<'\t'<<(mCXX[i]+0.000)/coverage <<endl;//the result is normalized		
		outfile.close();
				}

			}
    return 0;}
