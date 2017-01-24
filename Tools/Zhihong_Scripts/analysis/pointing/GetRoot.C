void GetRoot(){

	TString Target = ""; cerr<<"--- Target = "; cin >>Target;
	TString Kin = ""   ; cerr<<"--- Kin = ";    cin >> Kin;
	TString Arm = ""   ; cerr<<"--- Arm = ";    cin >> Arm;
	
	TString infile_name = Form("%s_Yield_%s_Kin%s_ep.dat",Target.Data(), Arm.Data(), Kin.Data());
	TString infile_dum = Form("Dummy_Yield_%s_Kin%s_ep.dat", Arm.Data(), Kin.Data());
 
 	TString rootfile_name = infile_name;
	rootfile_name.ReplaceAll(".dat",".root");
	TFile *f1 = new TFile(rootfile_name.Data(),"recreate");
    TTree *T = new TTree("T","A new Tree");
	
    /* Angle {{{*/
    double Theta = 0.0;
    if(Kin.Contains("3.1")||Kin.Contains("3.2"))
		Theta = 21.00;
	else if(Kin.Contains("4.1")||Kin.Contains("4.2"))
		Theta = 23.00;
	else if(Kin.Contains("5.0")||Kin.Contains("5.05")||Kin.Contains("5.1")||Kin.Contains("5.2"))
		Theta = 25.00;
	else if(Kin.Contains("6.5"))
		Theta = 28.00;
    else
		cerr<<"*** ERROR, Wrong Kinematics"<<endl;
	/*}}}*/

	const double E0 = 3.356; //GeV
	double xbj, Ep, VZ, Y, Y_Err;
	double xbj_dum,Ep_dum, VZ_dum, Y_dum, Y_Err_dum;

	T->Branch("Ep",&Ep,"Ep/D");
	T->Branch("xbj",&xbj,"xbj/D");
	T->Branch("VZ",&VZ,"VZ/D");
	T->Branch("Y",&Y,"Y/D");
	T->Branch("Y_Err",&Y_Err,"Y_Err/D");

	T->Branch("Ep_dum",&Ep_dum,"Ep_dum/D");
	T->Branch("xbj_dum",&xbj_dum,"xbj_dum/D");
	T->Branch("VZ_dum",&VZ_dum,"VZ_dum/D");
	T->Branch("Y_dum",&Y_dum,"Y_dum/D");
	T->Branch("Y_Err_dum",&Y_Err_dum,"Y_Err_dum/D");
	
    ifstream infile(infile_name.Data());
	while(!(infile.eof())){
      infile >> Ep >> VZ >> Y >> Y_Err;

	  xbj = 4.0*E0*Ep*pow( sin(Theta*3.14/180.*0.5), 2)/(2.0*(E0-Ep)*0.938272);

	  ifstream inf_dum(infile_dum.Data());
	  while(!(inf_dum.eof())){
		  inf_dum >> Ep_dum >> VZ_dum >> Y_dum >> Y_Err_dum;
		  if(fabs(Ep-Ep_dum)<0.00001&&fabs(VZ-VZ_dum)<0.00001){
		  xbj_dum = 4.0*E0*Ep_dum*pow( sin(Theta*3.14/180.*0.5), 2)/(2.0*(E0-Ep_dum)*0.938272);
		  cerr<<Form("--Dummy: Ep=%e, xbj=%e, VZ=%e, Y=%e", Ep_dum,xbj_dum,VZ_dum,Y_dum)<<endl;
		  cerr<<Form("--%s: Ep=%e, xbj=%e, VZ=%e, Y=%e", Target.Data(), Ep,xbj, VZ, Y)<<endl;
		  T->Fill();
		  break;
		  }
	  }
	  inf_dum.close();

	  //T->Fill();
	}
    infile.close();


	f1->cd(); T->Write(); f1->Close();

}
