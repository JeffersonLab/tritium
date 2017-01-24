  #define NOCASE 1

  /*C/C++ Includes{{{*/
  #include <stdio.h>
  #include <string>
  #include <iostream>
  #include <fstream>
  #include <vector>
  #include <algorithm>
  #include <map>
  #include <cmath>
  /*}}}*/
  /*ROOT Includes{{{*/
  #include <TSystem.h>
  #include <TString.h>
  #include <TStyle.h>
  #include <Riostream.h>
  #include "TObjString.h"
  #include <TNamed.h>
  #include <TPRegexp.h>
  #include <TObjArray.h>
  #include <TChain.h>
  #include <TMath.h>
  #include <TH1.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TFile.h>
  #include <TROOT.h>
  #include <TF1.h>
  #include <TGraph.h>
  #include <TGraphErrors.h>
  #include <TCanvas.h>
  #include <TDatime.h>
  #include <TError.h>
  #include <TVirtualFitter.h>
  #include <TSQLServer.h>
  #include <TSQLResult.h>
  #include <TSQLRow.h>
  #include <TCut.h>
  #include <TMultiGraph.h>
  #include <TCutG.h>
  #include <TLorentzVector.h>
  #include <TMath.h>
  #include <TRandom3.h>
  //#include <TMatrix.h>
  /*}}}*/

  //*ANALYZER Includes{{{*/
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaAnalysisObject.h"
 // #include "THaRasteredBeam.h"
 // #include "THaDecData.h"
 // #include "THaScalerGroup.h"
  //#include <THaApparatus.h>
  //#include <THaHRS.h>
  //#include <THaShower.h>
  //#include <THaCherenkov.h>
  //#include <THaScintillator.h>
  //#include <THaReactionPoint.h>
  //#include <THaElectronKine.h>
  //#include <THaGoldenTrack.h>
  //#include <THaExtTarCor.h>

  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaUnRasteredBeam.h"
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaNormAna.h"
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaAnalyzer.h"
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaRun.h"
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaRunParameters.h"
  #include "/work/halla/e08014/disk1/yez/analyzer/src/THaEvent.h"
  //*}}}*/

  using namespace std;

//All Energy Unit are in MeV, all Length Units are in cm!
  //const Double_t PROTON_MASS=0.938272*1000.00; //MeV Mass of proton

  /*Customization{{{*/
  const TString XGT2_DIR="/work/halla/e08014/disk1/yez/Cross_Section/"; //this is for inputfiles,outfiles,logfiles
  const TString RAWFILES_DIR="/w/cache/mss/halla/e07006/raw/";    //this is for raw
  const TString ROOTFILES_DIR="/work/halla/e08014/disk1/Rootfiles";     //this is for rootfiles
  const TString ROOTFILES_NAME = "e08014"; //such as e08014_4299.root
  //const TString ROOTFILES_DIR="/w/halla-2/e08014/Rootfiles_short";     //this is for rootfiles
  //const TString ROOTFILES_NAME = "e08014_short"; //such as e08014_short_4299.root
  const TString DATAQUALITYFILES_DIR="/work/halla/e08014/disk1/yez/Cross_Section/";//this is for dataqualityfiles
  const TString SAMCFILES_DIR="/work/halla/e08014/disk1/yez/Monte_Calo/SAMC/rootfiles/";    //this is for SAMC Rootfiles
  const TString XEMC_Table_Path = "/work/halla/e08014/disk1/yez/Monte_Calo/Cross_Section_Tables/";
  bool XEMC_Model_New = kTRUE;
  const TString XEMC_ByKin_Sub = "0920/";
  //const TString XEMC_Table_Path = "/w/halla-2/e08014/yez/Monte_Calo/Cross_Section_Tables/0920/";
  ofstream outlog; 
  Double_t Beam_Trip_Cut = 0.0;
  Double_t Target_Length = 0.0;
  Double_t Target_Length_Default = 0.20; // Default length for x>2 cryo-targets are 20cm;
  
  const Double_t GCPMT_For_1Photon=100;
  /*}}}*/
  /*class XGT2_Logger{{{*/
  class XGT2_Logger
  {
  public:
    /*static XGT2_Logger* Instance(){{{*/
    static XGT2_Logger* Instance()
    {
      // Only allow one instance of class to be generated.
      if (!m_pInstance)
        {
          Printf("@@@@@@ Creating XGT2_Logger System.");
          m_pInstance = new XGT2_Logger;
        }
      return m_pInstance;
    }
    /*}}}*/
    /*static XGT2_Logger* Destroy(){{{*/
    static XGT2_Logger* Destroy()
    {
      // Only allow one instance of class to be generated.
      if ( m_pInstance)
        delete m_pInstance;
      Printf("@@@@@@ Deleting XGT2_Logger System.");
      m_pInstance=NULL;
      return m_pInstance;
    }
    /*}}}*/
    /*const TString GetLogFileName(){{{*/
    const TString GetLogFileName()
    {
      return LogFile_Name;
    }
    /*}}}*/

    void SetLogFileName(const char* aLogFile_Name)
    {
      LogFile_Name=aLogFile_Name;
    }


    /*void WriteToLogFile(const char* aSrcFile, const char* aSrcLine,const char* aLog){{{*/
    void WriteToLogFile(const char* aSrcFile, int aSrcLine,const char* aLog)
    {
      ofstream outfile(LogFile_Name.Data(),ios::app);
      if ( !outfile.good() )
        {
          Error(__FILE__,Form("At line %d. %s cannot be appended.",__LINE__,LogFile_Name.Data()));
          return;
        }
      TDatime stamp;
      outfile<<Form("%s File-\"%s\" Line-%d: %s",stamp.AsString(),aSrcFile,aSrcLine,aLog)<<endl;
      outfile.close();
    }
    /*}}}*/


  private:
    XGT2_Logger(){};  // Private so that it can  not be called
    XGT2_Logger(XGT2_Logger const&){};             // copy constructor is private
    XGT2_Logger& operator=(XGT2_Logger const&){};  // assignment operator is private
    static XGT2_Logger* m_pInstance;
    TString LogFile_Name;
  }
  ;/*}}}*/

  XGT2_Logger* XGT2_Logger::m_pInstance = NULL;


  Bool_t IsBackground=kFALSE;

  /*Constants and enum{{{*/
  const Double_t Pi=TMath::Pi();
  const Double_t Na=TMath::Na();             //Avogadro's Number
  const Double_t Qe=TMath::Qe();
  const Int_t    DQ_Name_Length=60;         //for DataQuality Output Name Length
  
  enum E_INPUT_VARS_Def
    {
      k_INPUT_RUNNO,
      k_INPUT_ARM,
      k_INPUT_SPEC_SET_MOM,                   //GeV
      k_INPUT_TARGET_NAME,
      k_INPUT_TARGET_THICKNESS,              //density*Length kg/m^2
      k_INPUT_TARGET_THICKNESS_STAT_ERR,     //density*Length kg/m^2
      k_INPUT_TARGET_A,                      //Atomic mass kg/mol
      k_INPUT_TARGET_D,                      //target density kg/m^3
      k_INPUT_ISEXTGT,                       //1=extended target 0=not
      k_INPUT_MAIN_TRIGGER,
      k_INPUT_S_TREE_NAME,
      k_INPUT_RUNNOCHAIN
    };
  /*}}}*/

  /*Survey Or Constants{{{*/

  //Fix Here, ZYE
  const Double_t Harp_Nom[2][3]=
    {
      //if no offset,beam(0,0),the position of wire is nominal position.
      {23112.8,49637.8,72772.8}, //Harp nominal wire position for A (not encoder position)
      //{39154,81488,118411} encoder nominal position.
      {22319.4,48739.4,71924.4}, //Harp nominal wire position for B
      //{39068,81261,118287} encoder nominal position.
    };
  const Double_t Harp_Encoder_To_Wire_Coef[2][2]=
    {
      {2266,1.596}, //Harp coef to translate encoder to wire pos for A
      {3424,1.597}, //Harp coef to translate encoder to wire pos for B
    };
  const Double_t Harp_Offset[2][3]=
    {
      //mm
      {0.23,-0.38,0},      //Harp survey offset translate wire pos to beam pos for A, See survey DT_A1140.xml, Dy is inversed.
      {-0.64,-0.34,222.7}, //Harp survey offset translate wire pos to beam pos for B, See survey DT_A1140.xml, Dy is inversed.
      //222.7 is additional offset for B which is found from previous experiments. So I add 0 for A.
    };

  const TString VDC_PLANE_NAME[]={"u1","v1","u2","v2",""};

  /*}}}*/

  /*class XGT2_VAR{{{*/
  class XGT2_VAR
  {
  public:
    XGT2_VAR()
      {
        SetValue();
      };
    XGT2_VAR(XGT2_VAR const&){};
    XGT2_VAR& operator=(XGT2_VAR const&){};
    XGT2_VAR& operator*=(const Double_t& aNum)
      {
        this->Value *= aNum;
        this->Stat_Err *= aNum;
        this->Sys_Err *= aNum;
        return *this;
      };

    Double_t Value;//value of this variable
    Double_t Stat_Err;//statistical error of this variable
    Double_t Sys_Err;//systematic error of this variable

    void SetValue
      (
      const Double_t& aValue=0,
      const Double_t& aStat_Err=0,
      const Double_t& aSys_Err=0
      )
    {
      Value=aValue;
      Stat_Err=aStat_Err;
      Sys_Err=aSys_Err;
    };

    void Print()
    {
      Printf("%-*s%e",20,"Value:",Value);
      Printf("%-*s%e",20,"Stat_Err:",Stat_Err);
      Printf("%-*s%e",20,"Sys_Err:",Sys_Err);
    }
  };
  /*}}}*/

  /*inline vector<TString> Tokens(TString aline,TString aDelim=" "){{{*/
  inline vector<TString> Tokens(TString aline,TString aDelim=" ")
  {
    Int_t i;
    TObjArray* InObjArray;
    TObjString* os;
    TString s;
    vector<TString> OutStringVec;
    OutStringVec.clear();

    InObjArray=aline.Tokenize(aDelim);
    for ( i=0; i<InObjArray->GetEntriesFast(); i++ )
      {
        os=(TObjString*)InObjArray->At(i);
        s=os->GetString();
        OutStringVec.push_back(s);
      }
    return OutStringVec;
  }
  /*}}}*/


  /*inline void gGet_Acc_Bin(const Int_t& aRunNo,Double_t aAccBin[4][3]){{{*/
  inline void gGet_Acc_Bin(const TString& aAccBin_File,Double_t aAccBin[4][3])
  {
    //aAccBin[4][3]
    //[4] 0:dp  1:nu  2:Q2 3:xbj
    //[3] 0:min 1:max 2:stepsize
      
    TString infilename=Form("%s/Bin/%s",XGT2_DIR.Data(),aAccBin_File.Data());
    TString aPath = Form("%s/Acc_Bin/",XGT2_DIR.Data());

  /*   gCheckDirectory(aPath);*/
    if ( gSystem->FindFile(aPath.Data(),infilename) )
    {
      ifstream acc_bin_infile(infilename.Data());

    TString varname;
    Double_t min,max,stepsize;

    varname.ReadLine(acc_bin_infile);//read first line
    while ( acc_bin_infile>>varname>>min>>max>>stepsize )
      {
        if ( varname.Contains("dp") )
          {
            aAccBin[0][0]=min;
            aAccBin[0][1]=max;
            aAccBin[0][2]=stepsize;
          }
        else if ( varname.Contains("nu") )
          {
            aAccBin[1][0]=min;
            aAccBin[1][1]=max;
            aAccBin[1][2]=stepsize;
          }      
        else if ( varname.Contains("qsq") )
          {
            aAccBin[2][0]=min; //in GeV/c^2
            aAccBin[2][1]=max; 
            aAccBin[2][2]=stepsize;
          }
        else if ( varname.Contains("xbj") )
          {
            aAccBin[3][0]=min;
            aAccBin[3][1]=max;
            aAccBin[3][2]=stepsize;
          }
        if ( !acc_bin_infile.good() )
          break;
      }
    acc_bin_infile.close();

    }
    else{
      outlog <<"File Cannot be found "<< infilename.Data() << "at "<<aPath <<endl;
      cerr <<"File Cannot be found "<< infilename.Data() << "at "<<aPath <<endl;
    }
  }
  /*}}}*/

 /*inline TChain* gAddTree(const Int_t aRunNo,const TString& aArm,TString aTreeName="T"){{{*/
  inline TChain* gAddTree(const Int_t aRunNo,TString aTreeName="T")
  {
    TChain* aTree=new TChain(aTreeName.Data(),"Main Tree");
    Int_t index=0;
    TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
    if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
      {
	aTree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
      }
    else{  
      outlog<<Form("       *ERROR*, no rootfile founded for Run#%d",aRunNo)<<endl;
      cerr<<Form("       *ERROR*, no rootfile founded for Run#%d",aRunNo)<<endl;
      return 0;
    }
    index++;
    File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
    while ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
      {
	aTree->Add(Form("%s/%s_%d_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo,index));
	index++;
	File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
      }
    aTree->LoadTree(0);
  
    return aTree;
  }
  /*}}}*/

  /*inline TChain* gGetTree(const vector<Int_t>& aRunNoChain,const TString& aArm,TString aTreeName="T"){{{*/
  inline TChain* gGetTree(const vector<Int_t>& aRunNoChain,const TString& aArm,TString aTreeName="T")
  {
    
    outlog<<Form("       Obtaining Tree for %s Arm",aArm.Data())<<endl;
    cerr<<Form("       Obtaining Tree for %s Arm",aArm.Data())<<endl;
    
    TChain* Tree=new TChain(aTreeName.Data(),"Main Tree");
    for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
      {
        Int_t aRunNo=aRunNoChain[i];
        Int_t index=0;
        TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
        if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
          {
            outlog<<Form("      %s->Add(%s/%s_%d.root)",aTreeName.Data(),ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo)<<endl;
            cerr<<Form("      %s->Add(%s/%s_%d.root)",aTreeName.Data(),ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo)<<endl;
            Tree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
          }
        else
          return 0;
        index++;
        File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
        while ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
          {
            outlog<<Form("         Add %s/%s_%d_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo,index)<<endl;
            cerr<<Form("         Add %s/%s_%d_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo,index)<<endl;
            Tree->Add(Form("%s/%s_%d_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo,index));
            index++;
            File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
          }
        //Have to add this line because there is a bug
        //in the ROOT
        //please see
        //http://root.cern.ch/phpBB2/viewtopic.php?t=6145&highlight=ttreeformula+findleafforexpression
        Tree->LoadTree(0);
      }
    return Tree;
  }
  /*}}}*/

  /*inline TChain* gGetRootfiles(const TString& aRootFileName,TString aTreeName="T"){{{*/
  inline TChain* gGetRootFiles(const TString& aRootFilePath,const TString& aRootFileName,TString aTreeName="T")
  {
    outlog<<Form("       Obtaining Rootfiles for %s ",aRootFileName.Data())<<endl;
    cerr<<Form("       Obtaining Rootfiles for %s ",aRootFileName.Data())<<endl;
    TChain* Tree=new TChain(aTreeName.Data(),"Main Tree");
  
    TString File_Form=aRootFileName.Data();
   
    if ( gSystem->FindFile(aRootFilePath.Data(),File_Form) )
      {
        outlog<<("%s->Add(%s/%s)",aTreeName.Data(),aRootFilePath.Data(),aRootFileName.Data())<<endl;
        cerr<<("%s->Add(%s/%s)",aTreeName.Data(),aRootFilePath.Data(),aRootFileName.Data())<<endl;
        Tree->Add(Form("%s/%s",aRootFilePath.Data(),aRootFileName.Data()));
      }
    else
      return 0;
  
   Int_t index=1;
    File_Form=File_Form.ReplaceAll(".root",Form("_%d.root",index));
    while ( gSystem->FindFile(aRootFilePath.Data(),File_Form) )
       {
         outlog<<Form("%s->Add(%s/%s)",aTreeName.Data(),aRootFilePath.Data(),File_Form.Data())<<endl;
         cerr<<Form("%s->Add(%s/%s)",aTreeName.Data(),aRootFilePath.Data(),File_Form.Data())<<endl;
         Tree->Add(Form("%s/%s",aRootFilePath.Data(),File_Form.Data()));      
         index++;
         File_Form=File_Form.ReplaceAll(Form("_%d.root",index-1),Form("_%d.root",index));
       }
        //Have to add this line because there is a bug
        //in the ROOT
        //please see
        //http://root.cern.ch/phpBB2/viewtopic.php?t=6145&highlight=ttreeformula+findleafforexpression
    Tree->LoadTree(0);
    return Tree;
  }
  /*}}}*/


  /*inline Double_t gCal_Simple_Error(map <Double_t,Double_t>& aX,const Double_t& aY){{{*/
  inline Double_t gCal_Simple_Error(map <Double_t,Double_t>& aX,const Double_t& aY)
  {
    //calculate propagation error from simple formula
    //aX(Value,Error) aY=aX1*aX2.../(aXN*aX(N+1)...)
    Double_t error=0;
    for ( map<Double_t,Double_t>::iterator it=aX.begin(); it!=aX.end(); it++ )
      {
        if ( fabs(it->first)>1e-30 )
          error+=it->second*it->second/(it->first*it->first);
      }
    error=aY*sqrt(error);

    return error;
  }
  /*}}}*/

  /*inline vector<Int_t> gGet_RunNoChain(const Int_t& aRunNo,const TString& aString){{{*/
  inline vector<Int_t> gGet_RunNoChain(const Int_t& aRunNo,const TString& aString)
  {
    vector<Int_t> output;
    output.push_back(aRunNo);
    vector<TString> input=Tokens(aString);
    unsigned int i;
    for ( i=0; i<input.size(); i++ )
      {
        if ( isdigit(*input[i].Data()) || input[i].BeginsWith("-") )
          {
            vector<TString> tmp=Tokens(input[i],'-');
            if ( tmp.size()<=2 )
              {
               // Int_t FirstRunNo;
                Int_t LastRunNo;
                Int_t thisrunno;
                size_t index=input[i].First('-');
                if ( index==0 )
                  {
                    thisrunno=aRunNo;
                    LastRunNo=tmp[0].Atoi();
                  }
                else if ( index<input[i].Length() )
                  {
                    thisrunno=tmp[0].Atoi();
                    LastRunNo=tmp[1].Atoi();
                  }
                else
                  {
                    thisrunno=tmp[0].Atoi();
                    LastRunNo=thisrunno;
                  }
                if ( thisrunno>LastRunNo )
                  {
                    thisrunno ^= LastRunNo;
                    LastRunNo ^= thisrunno;
                    thisrunno ^= LastRunNo;
                  }
                while ( thisrunno<=LastRunNo )
                  {
                    if ( thisrunno!=aRunNo )
                      output.push_back(thisrunno);
                    thisrunno++;
                  }
              }
          }
      }
    for ( i=0; i<output.size(); i++ )
      {
        if ( i==0 && output.size()>1 ){
          outlog<<"      RunNoChain is "<<output[i];
          cerr<<"      RunNoChain is "<<output[i];
		}
        else if ( i==0 && i==(output.size()-1) ){
          outlog<<"      RunNoChain is "<<output[i]<<endl;
          cerr<<"      RunNoChain is "<<output[i]<<endl;
		}
        else if ( i==(output.size()-1) ){
          outlog<<" "<<output[i]<<endl;
          cerr<<" "<<output[i]<<endl;
		}
        else{
          outlog<<" "<<output[i];
          cerr<<" "<<output[i];
		}
      }
    return output;
  }
  /*}}}*/

  /*inline const void gGet_Leaf_Mean_and_RMS(TTree* aTree,const TString& aLeaf,const TString& aCut,Double_t& oMean,Double_t& oRMS){{{*/
  inline const void gGet_Leaf_Mean_and_RMS(TTree* aTree,const TString& aLeaf,const TString& aCut,Double_t& oMean,Double_t& oRMS,const TString& aFit="NULL")
  {
    aTree->Draw(aLeaf.Data(),aCut.Data(),"goff");
    TH1* htemp=(TH1*)gROOT->FindObject("htemp");
    if ( aFit!="NULL" )
      {
        htemp->Fit(aFit.Data());
        oMean=htemp->GetFunction(aFit.Data())->GetParameter("Mean");
        oRMS=htemp->GetFunction(aFit.Data())->GetParameter("Sigma");
      }
    else
      {
        oMean=htemp->GetMean();
        oRMS=htemp->GetRMS();
      }
    delete htemp;
  }
  /*}}}*/

  /*inline void gGet_Leaf_Range(TTree* aTree,const TString& aLeaf,const TString& aCut,Double_t & oLow,Double_t & oHigh){{{*/
  inline void gGet_Leaf_Range(TTree* aTree,const TString& aLeaf,const TString& aCut,Double_t & oLow,Double_t & oHigh)
  {
    aTree->Draw(aLeaf.Data(),aCut.Data(),"goff");
    TH1* htemp=(TH1*)gROOT->FindObject("htemp");
    oLow=htemp->GetXaxis()->GetXmin();
    oHigh=htemp->GetXaxis()->GetXmax();
    delete htemp;
  }
  /*}}}*/

  /*inline void gGet_File_Lines(const TString& aFile_Name,Int_t& oLineNo,Int_t& oCommentLineNo){{{*/
  inline void gGet_File_Lines(const TString& aFile_Name,Int_t& oLineNo,Int_t& oCommentLineNo)
  {
    ifstream file(aFile_Name.Data());
    oLineNo=0;
    oCommentLineNo=0;
    TString content;
    while ( content.ReadLine(file) )
      {
        if ( content.BeginsWith("#") )
          oCommentLineNo++;
        else
          oLineNo++;
        if ( !file.good() )
          break;
      }
    file.close();
  }
  /*}}}*/

  /*const inline void gCheckDirectory(const TString& aDir){{{*/
  const inline void gCheckDirectory(const TString& aDir)
  {
    if ( !gSystem->OpenDirectory(aDir.Data()) )
      {
        gSystem->mkdir(aDir);
      }
  }
  /*}}}*/

  /*inline TString gGet_InputFile_Var(const TString& aFile_Name,const Int_t& aLine){{{*/
  inline TString gGet_InputFile_Var(const TString& aFile_Name,const Int_t& aLine)
  {
    ifstream file(aFile_Name.Data());
    TString content;
    Int_t i=0;
    TString outstring;
    while ( content.ReadLine(file) )
      {
        if ( i==aLine )
          {
            outstring=content;
            break;
          }
        i++;
        if ( !file.good() )
          break;
      }
    file.close();
    if ( outstring=="\"\"" )
      outstring="";
    return outstring;
  }
  /*}}}*/

  /*inline Int_t gAskRunNo(Int_t RunNo){{{*/
  inline Int_t gAskRunNo(Int_t RunNo)
  {
    if ( RunNo<0 )
      {
        outlog<<"Please Enter Run number ( -1-quit ): ";
        cerr<<"Please Enter Run number ( -1-quit ): ";
        cin>>RunNo;
      }
    return RunNo;
  }
  /*}}}*/

  /*inline TString gGetScriptname(const char* aFileName){{{*/
  inline TString gGetScriptname(const char* aFileName)
  {
    TString scriptname=gSystem->BaseName(aFileName);
    scriptname.ReplaceAll(".C","");
    scriptname.ReplaceAll("_inputfiles","");
    scriptname.ReplaceAll(".log","");
    scriptname.ReplaceAll(".dat","");
    scriptname.ReplaceAll(".txt","");
    scriptname.ReplaceAll(".root","");
    return scriptname;

  }
  /*}}}*/

  /*const inline void gPrintLine(TString in="-",Int_t length=30,Bool_t NewLine=kTRUE){{{*/
  const inline void gPrintLine(TString in="-",Int_t length=30,Bool_t NewLine=kTRUE)
  {
    for ( Int_t i=0; i<length; i++ )
      printf("%s",in.Data());
    if ( NewLine )
      printf("\n");
  }
  /*}}}*/

