#include "/home/yez/work/Cross_Section/SRC/XGT2_Main.h"

/*int Gen_CrossSection_RootFile(const TString& aFileName,const Int_t& aBaseSize,const TString& aLogFileName){{{*/
int Gen_CrossSection_RootFile(const TString& aTarget_Name,const TString& aFileName,const TString& aLogFileName)
{
  XGT2_Logger* log=XGT2_Logger::Instance();
  log->SetLogFileName(aLogFileName);
	
  /*Read InputFile{{{*/
  ifstream inputfile;
	
  const TString InputFileName=Form("%s/results/%s/%s",XGT2_DIR.Data(),aTarget_Name.Data(),aFileName.Data());
  inputfile.open(InputFileName.Data());
	
  if ( !inputfile.good() )
    {
      log->WriteToLogFile(__FILE__,__LINE__,Form("Error No InputFileName is input."));
      log->Destroy();
      return -2;
    }

  TString F_FileName=Form("%s/results/%s/%s.root",XGT2_DIR.Data(),aTarget_Name.Data(),gGetScriptname(aFileName).Data());
  TFile* F=new TFile(F_FileName.Data(),"recreate");
  TTree* T=new TTree("T",Form("%s Cross Section",gGetScriptname(aFileName).Data()));

  Int_t RunNo;
  Double_t P0,nu,XS,XS_Stat_Err,XS_Sys_Err,Q2,xbj,angle,Ei,Ep,tmp,Livetime,Ntg,Ni,Nf_EX,Nf_SAMC,Nf_SAMC_Gen;
  T->Branch("RunNo",&RunNo,"RunNo/I");
  T->Branch("XS",&XS,"XS/D");
  T->Branch("XS_Stat_Err",&XS_Stat_Err,"XS_Stat_Err/D");
  T->Branch("XS_Sys_Err",&XS_Sys_Err,"XS_Sys_Err/D");
  T->Branch("P0",&P0,"P0/D");
  T->Branch("nu",&nu,"nu/D");
  T->Branch("xbj",&xbj,"xbj/D");
  T->Branch("Q2",&Q2,"Q2/D");
  T->Branch("angle",&angle,"angle/D");
  T->Branch("Ei",&Ei,"Ei/D");
  T->Branch("Ep",&Ep,"Ep/D");
  T->Branch("Livetime",&Livetime,"Livetime/D");
  T->Branch("Ntg",&Ntg,"Ntg/D");
  T->Branch("Ni",&Ni,"Ni/D");
  T->Branch("Nf_EX",&Nf_EX,"Nf_EX/D");
  T->Branch("Nf_SAMC",&Nf_SAMC,"Nf_SAMC/D");
  T->Branch("Nf_SAMC_Gen",&Nf_SAMC_Gen,"Nf_SAMC_Gen/D");

  ifstream run_cs_file;run_cs_file.open(InputFileName.Data());
  TString comment;
  cerr << "Creating Rootfile ..." <<endl;
  comment.ReadLine(run_cs_file);
  Int_t Count =0;
  while ( run_cs_file>>XS>>XS_Stat_Err>>XS_Sys_Err>>xbj>>angle>>Ep>>nu>>Q2>>Ei>>P0>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_SAMC>>Nf_SAMC_Gen )
    {
      if ( angle<180 )
	{
	  Count++;
          cerr << Form("   -- XS = %e, Xbj = %f", XS, xbj ) <<endl;
	  T->Fill();
	}
      if ( !run_cs_file.good() )
	break;
    }
  run_cs_file.close();
  inputfile.close();

  cerr << " -----> Total Data Points = "<< Count <<endl;
  T->Write();

  delete T;
  delete F;

  log->WriteToLogFile(__FILE__,__LINE__,Form("Success. Result is saved at %s",F_FileName.Data()));
  log->Destroy();
  return 0;
}
/*}}}*/
