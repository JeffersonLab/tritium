////////////////////////////////
// Events Generator for SAMC
// -- Zhihong Ye
/////////////////////////////////

/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/

#define CHAR_LEN 1024

using namespace::std;
using namespace::TMath;

string gInput;  //Input filename
string gOutput; //Output filename
int gNEvent;    //Total Event

int getargs(int argc, char**argv);

int main(int argc,char** argv){
  
  int err = getargs(argc,argv);
  if(err<0){
    cerr<<"****** I can not recognize the parameters!"<<endl;
    return -1;
  }
  
  /*Define and Initialize Parameters*/
  double delta_dp=0.14;	//d(dp) dp full width for generator
  double delta_th=0.18;	//d(th) th full width for generator(tan(th))
  double delta_ph=0.09;	//d(ph) th full width for generator(tan(ph))
  double gaus_x_sigma=0;	//if Which_Beam_Profile==1/2 sigma of x beam for generator(cm)
  double gaus_y_sigma=0;	//if Which_Beam_Profile==1/2 sigma of y beam for generator(cm)
  int  Which_Beam_Profile=1; //0=Norm 1=Gaus 2=gaus+norm=triangle
  double raster_x_size=0;	//raster x full size for generator(cm)
  double raster_y_size=0;	//raster y full size for generator(cm)
  double beam_x_center=0;	//beam x center for generator(cm)
  double beam_y_center=0;	//beam y center for generator(cm)
  double z0=0;	//target center for generator(cm)
  double T_L=0;	//target length for generator(cm)
  double T_H=0;	//target height for generator(cm)
  /*}}}*/


  printf("Start to generate seeds ...\n");

  /*Read parameters{{{*/
  vector<string> inputdata;
  FILE* infile;
  infile=fopen(gInput.c_str(),"r");
  char buf[CHAR_LEN];
  char data[CHAR_LEN];
  int i,j,k;
  while ( fgets(buf,CHAR_LEN,infile) )
    {
      i=0;
      while ( buf[i]==' ' )
	{
	  i++;
	}
      if ( buf[i]!='#' )
	{
	  j=0;
	  while ( buf[i]!='#' && buf[i]!='\0' )
	    {
	      data[j]=buf[i];
	      i++; j++;
	    }
	  data[j]='\0';
	  while ( data[--j]==' ' || data[j]=='\t' )
	    {
	      //remove space or tab at the end of data
	      data[j]='\0';
	    }
	  inputdata.push_back(data);
	}
      //else it's comment, skipped
    }
  fclose(infile);
  
  k=0;
  delta_dp=atof(inputdata[k++].c_str());
  delta_th=atof(inputdata[k++].c_str());
  delta_ph=atof(inputdata[k++].c_str());
  Which_Beam_Profile=atoi(inputdata[k++].c_str());
  raster_x_size=atof(inputdata[k++].c_str());
  raster_y_size=atof(inputdata[k++].c_str());
  gaus_x_sigma=atof(inputdata[k++].c_str());
  gaus_y_sigma=atof(inputdata[k++].c_str());
  beam_x_center=atof(inputdata[k++].c_str());
  beam_y_center=atof(inputdata[k++].c_str());
  z0=atof(inputdata[k++].c_str());
  T_L=atof(inputdata[k++].c_str());
  T_H=atof(inputdata[k++].c_str());

  cerr<<Form("--- Delta=%f, Theta=%f, Phi=%f ", delta_dp, delta_th,delta_ph);

  /*}}}*/
  
  /*Generate Seed{{{*/
  printf("Generating and Saving %d seeds in %s...\n", gNEvent, gOutput.c_str());
  
  ofstream outf; outf.open(gOutput.c_str());
  double beam_x,beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen;
  gRandom->SetSeed(0);
  
  for(int iN=0;iN<gNEvent; iN++){
    if ( Which_Beam_Profile==0 ) {
      beam_x=beam_x_center+(gRandom->Rndm()-0.5)*raster_x_size;
      beam_y=beam_y_center+(gRandom->Rndm()-0.5)*raster_y_size;
    }
    else if ( Which_Beam_Profile==1 ) {
      beam_x=gRandom->Gaus(beam_x_center,gaus_x_sigma);
      beam_y=gRandom->Gaus(beam_y_center,gaus_y_sigma);
    }
    else if ( Which_Beam_Profile==2 ) {
      beam_x=beam_x_center+(gRandom->Rndm()-0.5)*raster_x_size;
      beam_y=beam_y_center+(gRandom->Rndm()-0.5)*raster_y_size;
      beam_x=gRandom->Gaus(beam_x,gaus_x_sigma);
      beam_y=gRandom->Gaus(beam_y,gaus_y_sigma);
    }
    th_tg_gen=(gRandom->Rndm()-0.5)*delta_th;
    ph_tg_gen=(gRandom->Rndm()-0.5)*delta_ph;
    dp_gen=(gRandom->Rndm()-0.5)*delta_dp;

    reactz_gen=z0+(gRandom->Rndm()-0.5)*T_L;
    
	outf <<Form("%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f",
		beam_x, beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen)
	 <<endl;

    if(!((iN+1)%100000))
      printf("--- %d events generated ...\n", iN+1);
  }
  /*}}}*/
  
}

int getargs(int argc,char** argv)
{
  char* argptr;
  TString aPar;
  bool noStop;  
  int err=-1;

  for(int i=0;i<argc;i++){
    argptr = argv[i];
    if(*argptr=='-'){
      argptr++;
      switch (*argptr){
      case 'I':
	if(*(++argptr))
	  gInput = argptr;
	err=0;
	break;
      case 'O':
	if(*(++argptr))
	  gOutput = argptr;
	err=0;
	break;
      case 'N':
	if(*(++argptr)){
	  string gNC = argptr;
	  gNEvent = atoi(gNC.c_str());
	}
	err=0;
	break;
      case 'h':
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I Input file to define parameters" <<endl;
	cerr <<"-O Output file to save seeds" <<endl;
	cerr <<"-N Total Number of Event to generate" <<endl;
	cerr <<endl<<"================================================"<<endl;
	noStop=false;
	err=-1;
	goto OUT;
	break;
      default:
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I Input file to define parameters" <<endl;
	cerr <<"-O Output file to save seeds" <<endl;
	cerr <<"-N Total Number of Event to generate" <<endl;
	cerr <<endl<<"================================================"<<endl;
	err=-1;
	break;
      }

    }
    noStop=true;
  OUT:if(!noStop){break;}
    continue;
  }
  return err; 
}
