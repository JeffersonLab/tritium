void GenRoot(){

	ifstream infile1("XS_C12_T21_0_51.table");
	ifstream infile2("XS_C12_T23_0_54.table");
	ifstream infile3("XS_C12_T25_0_83.table");
	ifstream infile4("XS_C12_T28_0_34.table");
	TFile *file = new TFile("C12_XS.root","recreate");


	TTree *T = new TTree("T","A new Tree");

	//	outfile << Form("%5d %3d %6.3f %10.7f %6d %10.7f %6d %16.6e %16.6e %16.6e", A,Z,Theta,mom,j,theta,i,cs_born,cs_rad,tail)<<endl;

	const double E0 = 3.356;
	const double Mass = 0.938272;
	Double_t Angle, Ep, Theta, xbj, XS_Born, XS_Rad, tail,RC,Q2;
	Int_t A, Z, i,j;

	T->Branch("E0", &E0, "E0/D");
	T->Branch("Ep", &Ep, "Ep/D");
	T->Branch("xbj", &xbj, "xbj/D");
	T->Branch("Q2", &Q2, "Q2/D");
	T->Branch("Angle", &Angle, "Angle/D");
	T->Branch("Theta", &Theta, "Theta/D");
	T->Branch("XS_Born", &XS_Born, "XS_Born/D");
	T->Branch("XS_Rad", &XS_Rad, "XS_Rad/D");
	T->Branch("tail", &tail, "tail/D");
	T->Branch("RC", &RC, "RC/D");

	T->Branch("A", &A, "A/I");
	T->Branch("Z", &Z, "Z/I");
	T->Branch("i", &i, "i/I");
	T->Branch("j", &j, "j/I");

	TString Com;
	infile1 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile1.eof())){
		infile1 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born >> XS_Rad >> tail;
		RC = XS_Rad / XS_Born;
		Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile1.close();

	infile2 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile2.eof())){
		infile2 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born >> XS_Rad >> tail;
		RC = XS_Rad / XS_Born;
		double Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile2.close();

	infile3 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile3.eof())){
		infile3 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born >> XS_Rad >> tail;
		RC = XS_Rad / XS_Born;
		double Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile3.close();

	infile4 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile4.eof())){
		infile4 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born >> XS_Rad >> tail;
		RC = XS_Rad / XS_Born;
		double Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile4.close();

	file->cd(); T->Write(); file->Close();

}
