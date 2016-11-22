//strand considered
//Dr. Feng's algorithm
//For formatted BSseeker output
//report progress
//use #define
///02	+	12583083	CHH	0	1
///02	+	12583084	CHH	0	1
///02	+	12583087	CHH	0	1
   
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

const int Gwidth = 20; // # of bins for gene region
#define Twidth (Gwidth*3)

using namespace std;//convert BS-seeker's chr coordinate
int mystoi(string chr_c){
                int chr =0;
			if(chr_c=="01")chr=1;
		else if(chr_c=="02")chr=2;
		else if(chr_c=="03")chr=3;
		else if(chr_c=="04")chr=4;
		else if(chr_c=="05")chr=5;
		else if(chr_c=="06")chr=6;
		else if(chr_c=="07")chr=7;
		else if(chr_c=="08")chr=8;
		else if(chr_c=="09")chr=9;
		else if(chr_c=="10")chr=10;
		else if(chr_c=="11")chr=11;
		else if(chr_c=="12")chr=12;
		else if(chr_c=="13")chr=13;
		else if(chr_c=="14")chr=14;
		else if(chr_c=="15")chr=15;
		else if(chr_c=="16")chr=16;
		else if(chr_c=="17")chr=17;
		else if(chr_c=="18")chr=18;
		else if(chr_c=="19")chr=19;
		else if(chr_c=="20")chr=20;
		else if(chr_c=="21")chr=21;
		else if(chr_c=="22")chr=22;
		else if(chr_c=="23")chr=23;
		else if(chr_c=="24")chr=24;
		return chr;
}

int STRAND(char strand,char str){//distinguish which strand, anti or sense
	if(strand=='+'&&str=='+'||strand=='-'&&str=='-')return 0; 
	else if(strand=='+'&&str=='-'||strand=='-'&&str=='+')return 1;
	else{
		cout<<"strand wrong or str wrong";
		exit(1);
	    }
				}
				
void iofile(string pf, string rf, double mCG[2][Twidth],double mCHG[2][Twidth],double mCHH[2][Twidth],int wCG[2][Twidth],int wCHG[2][Twidth], int wCHH[2][Twidth]){
     //////////////////////////////////////////////////////////////////////////
     /////////                                                      ///////////
     /////////                                                      ///////////
     /////////             left    interesting   right              ///////////
     /////////      2*a-b -------a-------------b----------2*b-a     ///////////
     /////////                                                      ///////////
     //////////////////////////////////////////////////////////////////////////
     
     //////////////////////////////////////////////////
     //////////part1: for the left-to-interest region (2a-b ~ a)
     ///////////////////////////////////////////////////
     
    ifstream pos;    //the position (mapped)
    pos.open(pf.c_str());
    if(pos.fail()){cout<<"position fail";exit(1);}
        
    ifstream region;   //the interested region 
    region.open(rf.c_str());
    if(region.fail()){cout<<"bed fail";exit(1);}


    int c=0;   //chromosome int
    	string C;// chromosome string
    char strand;
    int p;     //position in the position file
    string context;
    int w;//w is the # of reads acting as unmethylated ones
    double mc;

    int chr,a,b;    //chromosome,start,end in the REGION file
    char str;  //strand in region file
    const int Con=0;//the constant to take into account of the length
    int clr = 0;//count the lines checked in region file
    int clp = 0;//count the lines checked in position file

    region>>chr>>str>>a>>b;
    clr++;
    pos>>C>>strand>>p>>context>>mc>>w;
    clp++; 
    c = mystoi(C);
    while(true){
//		if(mc+unmc<4){
//                      if(pos>>C==false)break;
//                      pos>>strand>>p>>context>>mc>>unmc>>Creads>>Greads>>perc;
//                      continue;
//                      } 

                if(c>chr||(c==chr&&(p+Con)>=a)){//if the position is on the right of the upstream end
                                          if(region>>chr==false)break;//reached the end of the region file 
                                          else 			{
								region>>str>>a>>b;
								clr++;
								}
						}// go to the next region
                else if(c==chr&&(2*a-b-1)<=(p+Con)&&(p+Con)<a){//if the position falls into the upstream region
			unsigned int pUPSTREAM = ((p+Con-2*a+b)*Gwidth+10)/(b-a+1);//(p+Con-2*a+b+0.5)*Gwidth/(b-a+1);                     
                     if(context=="CG"){
                         if(str=='+'){
				                      mCG[STRAND(strand,str)][pUPSTREAM]+=mc;//increase the corresponding signal counts(0-19)
				                      wCG[STRAND(strand,str)][pUPSTREAM]+=w;
				                      }
                          else if(str=='-'){
				               mCG[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=mc;
				               wCG[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=w;
                                           }
                          else {cout<<"strand neither + nor -"<<endl;exit(1);}
				      }
                     else if(context=="CHG"){
                         if(str=='+'){
				mCHG[STRAND(strand,str)][pUPSTREAM]+=mc;//increase the corresponding signal counts(0-19)
				wCHG[STRAND(strand,str)][pUPSTREAM]+=w;
				}
			else if(str=='-'){
				mCHG[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=mc;
				wCHG[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=w;
				}
			else {cout<<"strand neither + nor -"<<endl;exit(1);}
					     }
                     else if(context=="CHH"){
                         if(str=='+'){
				mCHH[STRAND(strand,str)][pUPSTREAM]+=mc;//increase the corresponding signal counts(0-19)
				wCHH[STRAND(strand,str)][pUPSTREAM]+=w;
				     }
			else if(str=='-'){
				mCHH[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=mc;
				wCHH[STRAND(strand,str)][Twidth-1-(pUPSTREAM)]+=w;
					 }
			else {cout<<"strand neither + nor -"<<endl;exit(1);}
					    }
                     else {cout<<"error in context"<<context;}
                                          if(pos>>C==false)break;
                     pos>>strand>>p>>context>>mc>>w;
			clp++;
                     c=mystoi(C);    
								}          
                                      //     pos>>C;
					//   if(C=="M")break;//reached the end of the position file
                                         //  else pos>>p>>strand>>context>>mc>>reads>>id;}//go to the next position but still examine the already fitted region
                else {//if the position is on the left of the upstream region
                     if(pos>>C==false)break;
                     pos>>strand>>p>>context>>mc>>w;
			clp++;
                     c=mystoi(C);
                     }
//			if(C=="M")break; //reached the end of position file
//                     else pos>>p>>strand>>context>>mc>>reads>>id;}//go to the next position 
                }
	
    region.close();
    pos.close();
cout<<"went through "<<clr<<" r_lines "<<clp << " p_lines for upstream" << endl;

//////////////////////////////////////////////////////////////
/////next part la!!!//////the original region (a, b)//////////////////////////////////
//////////////////////////////////////////////////////////
    clr = 0;
    clp = 0;

    ifstream pos1;    //the position (mapped)
    pos1.open(pf.c_str());
    if(pos1.fail()){cout<<"position fail";exit(1);}
        
    ifstream region1;   //the interested region 
    region1.open(rf.c_str());
    if(region1.fail()){cout<<"bed fail";exit(1);}

    c=0;   //chromosome
//    int p;     //position in the position file
//    int chr,a,b;    //chromosome,start,end in the region file
//    const int Con=0;//the constant to take into account of the length
/*Above is a very important modification!!!*/

//    int temp1,temp2;
    region1>>chr>>str>>a>>b;
    pos1>>C>>strand>>p>>context>>mc>>w;
	clr++;clp++;
    c = mystoi(C);
    while(true){
//		if(mc+unmc<4){
//                      if(pos1>>C==false)break;
//                      pos1>>strand>>p>>context>>mc>>unmc>>Creads>>Greads>>perc;
//                      continue;
//                      }
                if(c>chr||(c==chr&&(p+Con)>b)){//if the position is on the right of the region
                                          if(region1>>chr==false)break; // reached the end of the region file
                                          else 		{
							region1>>str>>a>>b;
							clr++;
							}
						} // go to the next region
                else if(c==chr&&a<=(p+Con)&&(p+Con)<=b){//if the position falls into the region
			unsigned int pWITHIN = ((p+Con-2*a+b)*Gwidth+10)/(b-a+1);//(p+Con-2*a+b+0.5)*Gwidth/(b-a+1);                    
                     if(context=="CG"){
                     	if(str=='+'){
				mCG[STRAND(strand,str)][pWITHIN]+=mc;//increase the corresponding signal density(20-39)
				wCG[STRAND(strand,str)][pWITHIN]+=w;
				    }
			else if(str=='-'){
				mCG[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=mc;
				wCG[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=w;
					 }
			else {cout<<"strand neither + nor -"<<endl;exit(1);}
					}
                     else if(context=="CHG"){
                     	if(str=='+'){
				mCHG[STRAND(strand,str)][pWITHIN]+=mc;//increase the corresponding signal density(20-39)
				wCHG[STRAND(strand,str)][pWITHIN]+=w;
				    }
			else if(str=='-'){
				mCHG[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=mc;
				wCHG[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=w;
					 }
			else {cout<<"strand neither + nor -"<<endl;exit(1);}
					     }
                     else if(context=="CHH"){
                     	if(str=='+'){
				mCHH[STRAND(strand,str)][pWITHIN]+=mc;//increase the corresponding signal density(20-39)
				wCHH[STRAND(strand,str)][pWITHIN]+=w;
				    }
			else if(str=='-'){
				mCHH[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=mc;
				wCHH[STRAND(strand,str)][Twidth -1- (pWITHIN)]+=w;
					 }
			else {cout<<"strand neither + nor -"<<endl;exit(1);}
					    }
                     else {cout<<"error in context"<<context;}
                                          if(pos1>>C==false)break;
                     pos1>>strand>>p>>context>>mc>>w;
                     c=mystoi(C);
			clp++;
							}
                                 //          pos1>>C;
				//	   if(C=="M")break;//if the position file reaches the end
                                 //          else pos1>>p>>strand>>context>>mc>>reads>>id;}//go to the next position
                else {//if the position is on the left of the upstream region
                     if(pos1>>C==false)break;
                     pos1>>strand>>p>>context>>mc>>w;
                     c=mystoi(C);
			clp++;
                     }//			if(C=="M")break;//if the position file reaches its end
  //                   else pos1>>p>>strand>>context>>mc>>reads>>id;}//go to the next position
                }

    region1.close();
    pos1.close();			
	cout<<"went through "<<clr<<" r_lines "<<clp << " p_lines for upstream" << endl;

/////////////////////////////////////
//////////////final part !!!  the right-to-interest region(b, 2b-a)//////
////////////////////////////////////

    ifstream pos2;    //the position (mapped)
    pos2.open(pf.c_str());
    if(pos2.fail()){cout<<"position fail";exit(1);}
        
    ifstream region2;   //the interested region 
    region2.open(rf.c_str());
    if(region2.fail()){cout<<"bed fail";exit(1);}

	clr=0;
	clp=0;    
    c=0;   //chromosome
   // int p;     //position in the position file
    //int chr,a,b;    //chromosome,start,end in the region file
    //const int Con=0;//the constant to take into account of the length
/*Above is a very important modification!!!*/

    //int temp1,temp2;
    region2>>chr>>str>>a>>b;
    pos2>>C>>strand>>p>>context>>mc>>w;
	clr++;
	clp++;

    while(true){
//		if(mc+unmc<4){
//                      if(pos2>>C==false)break;
//                      pos2>>strand>>p>>context>>mc>>unmc>>Creads>>Greads>>perc;
//                      continue;
//                      }
                if(c>chr||(c==chr&&(p+Con)>2*b-a+1)){//on the right of the downstream
                                          if(region2>>chr==false)break;
                                          else 		{
							region2>>str>>a>>b;
							clr++;
							}
									}
                else if(c==chr&&b<(p+Con)&&(p+Con)<=2*b-a+1){//fall into the downstream
			unsigned int pDOWNSTREAM = ((p+Con-2*a+b)*Gwidth+10)/(b-a+1);//(p+Con-2*a+b+0.5)*Gwidth/(b-a+1);                     
                     if(context=="CG"){
                     	if(str=='+'){
				mCG[STRAND(strand,str)][pDOWNSTREAM]+=mc;//increase signal(40-59)
				wCG[STRAND(strand,str)][pDOWNSTREAM]+=w;
				    }
			else if(str=='-'){
				mCG[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=mc;
				wCG[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=w;
					 }
			else {cout<<"strand neither + nor -";exit(1);}
				      }
                     else if(context=="CHG"){
                     	if(str=='+'){
				mCHG[STRAND(strand,str)][pDOWNSTREAM]+=mc;//increase signal(40-59)
				wCHG[STRAND(strand,str)][pDOWNSTREAM]+=w;
				    }
			else if(str=='-'){
				mCHG[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=mc;
				wCHG[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=w;
					 }
			else {cout<<"strand neither + nor -";exit(1);}  
					    }                   
		     else if(context=="CHH"){
                     	if(str=='+'){
				mCHH[STRAND(strand,str)][pDOWNSTREAM]+=mc;//increase signal(40-59)
				wCHH[STRAND(strand,str)][pDOWNSTREAM]+=w;
				    }
			else if(str=='-'){
				mCHH[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=mc;
				wCHH[STRAND(strand,str)][Twidth -1-(pDOWNSTREAM)]+=w;
					 }
			else {cout<<"strand neither + nor -";exit(1);}
					    }
                     else {cout<<"error in context"<<context;}  
                                          if(pos2>>C==false)break;
                     pos2>>strand>>p>>context>>mc>>w;
                     c=mystoi(C);
			clp++;
							     }	                                         
                           //                pos2>>C;
				//	   if(C=="M")break;
                                  //         else pos2>>p>>strand>>context>>mc>>reads>>id;}
                else {//if the position is on the left of the upstream region
                     if(pos2>>C==false)break;
                     pos2>>strand>>p>>context>>mc>>w;
			clp++;                     
			c=mystoi(C);
                     }	//		if(C=="M")break;//no more postion to fit
          //           else pos2>>p>>strand>>context>>mc>>reads>>id;}//go to next position
                }
                
    region2.close();
    pos2.close();
	cout<<"went through "<<clr<<" r_lines "<<clp << " p_lines for upstream" << endl;                      
   }
	
int main(){
    cout<<"show me the 24 methylation files' name, please."<<endl;
    string pf[24];
    for(int i = 0;i<24;i++)cin>>pf[i];//get methylation file order

    cout<<"How many region files then?(less than 1000 plz)"<<endl;
    int nr;
    cin >> nr;
    cout<<"And their names?"<<endl;
    string rf[1000];
    for(int j = 0;j<nr;j++)cin>>rf[j];//get region files' list

    for(int j = 0;j<nr;j++){
	double mCG[2][Twidth],mCHG[2][Twidth],mCHH[2][Twidth];
	int wCG[2][Twidth],wCHG[2][Twidth],wCHH[2][Twidth];
	for(int i=0;i<Twidth;i++){
		mCG[0][i]=0;
		wCG[0][i]=0;
		mCHG[0][i]=0;
		wCHG[0][i]=0;
		mCHH[0][i]=0;
		wCHH[0][i]=0;
		mCG[1][i]=0;
		wCG[1][i]=0;
		mCHG[1][i]=0;
		wCHG[1][i]=0;
		mCHH[1][i]=0;
		wCHH[1][i]=0;
			     }
	for(int k = 0;k<24;k++){
		cout << k << '\t' << j << endl;		
		iofile(pf[k],rf[j], mCG, mCHG, mCHH, wCG, wCHG, wCHH);
		
				}

	ofstream outfileCGO, outfileCHGO, outfileCHHO, outfileCGI, outfileCHGI, outfileCHHI;
	string ofCGO = rf[j] + pf[j]+ "CG(Con+).txt";
	string ofCHGO = rf[j] + pf[j]+ "CHG(Con+).txt";
	string ofCHHO = rf[j] + pf[j]+ "CHH(Con+).txt";
	string ofCGI = rf[j] + pf[j]+ "CG(Con-).txt";
	string ofCHGI = rf[j] + pf[j]+ "CHG(Con-).txt";
	string ofCHHI = rf[j] + pf[j]+ "CHH(Con-).txt";
	    
	outfileCGO.open(ofCGO.c_str());//out put file
    	for(int i=0;i<Twidth;i++){ 
		        outfileCGO<<i<<'\t'<<mCG[0][i]<<'\t'<<wCG[0][i]/*+0.00001*/ <<endl;
                             } 
    	outfileCGO.close();
    
    	outfileCHGO.open(ofCHGO.c_str());
    	for(int i=0;i<Twidth;i++){ 
		        outfileCHGO<<i<<'\t'<<mCHG[0][i]<<'\t'<<wCHG[0][i]/*+0.00001*/ <<endl;
                             } 
    	outfileCHGO.close();
    
    	outfileCHHO.open(ofCHHO.c_str());
    	for(int i=0;i<Twidth;i++){ 
		        outfileCHHO<<i<<'\t'<<mCHH[0][i]<<'\t'<<wCHH[0][i]/*+0.00001*/<<endl;
                             } 
    	outfileCHHO.close();

	outfileCGI.open(ofCGI.c_str());//out put file
    	for(int i=0;i<Twidth;i++){ 
		        outfileCGI<<i<<'\t'<<mCG[1][i]<<'\t'<<wCG[1][i]/*+0.00001*/ <<endl;
                             } 
    	outfileCGI.close();
    
    	outfileCHGI.open(ofCHGI.c_str());
    	for(int i=0;i<Twidth;i++){ 
		        outfileCHGI<<i<<'\t'<<mCHG[1][i]<<'\t'<<wCHG[1][i]/*+0.00001*/ <<endl;
                             } 
    	outfileCHGI.close();
    
    	outfileCHHI.open(ofCHHI.c_str());
    	for(int i=0;i<Twidth;i++){ 
		        outfileCHHI<<i<<'\t'<<mCHH[1][i]<<'\t'<<wCHH[1][i]/*+0.00001*/ <<endl;
                             } 
    	outfileCHHI.close();
			      }
    return 0;}
