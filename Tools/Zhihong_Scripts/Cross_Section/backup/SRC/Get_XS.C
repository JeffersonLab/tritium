#include "XGT2_Main.h"
#include "Cross_Section.h"

int Get_XS(){
  TString infile=NULL;
  TString logfile=NULL;
  cerr<<" What is your input file? ...... "; cin >> infile;
  cerr<<" What is your log file?   ...... "; cin >> logfile;
  
  int flag = CrossSection(infile,logfile);
  
  return flag;
}
