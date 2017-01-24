void GetSys_Rad_He4(){
	ifstream infile11("stat/He4/XS_He4_T21_0_51.table");
	ifstream infile12("stat/He4/XS_He4_T23_0_54.table");
	ifstream infile13("stat/He4/XS_He4_T25_0_83.table");

	ifstream infile21("stat_q2/He4/XS_He4_T21_0_51.table");
	ifstream infile22("stat_q2/He4/XS_He4_T23_0_54.table");
	ifstream infile23("stat_q2/He4/XS_He4_T25_0_83.table");

	ifstream infile31("stat_xbj/He4/XS_He4_T21_0_51.table");
	ifstream infile32("stat_xbj/He4/XS_He4_T23_0_54.table");
	ifstream infile33("stat_xbj/He4/XS_He4_T25_0_83.table");

	TFile *file = new TFile("He4_XS.root","recreate");

	TTree *T = new TTree("T","A new Tree");

	//	outfile << Form("%5d %3d %6.3f %10.7f %6d %10.7f %6d %16.6e %16.6e %16.6e", A,Z,Theta,mom,j,theta,i,cs_born,cs_rad,tail)<<endl;

	const double E0 = 3.356;
	const double Mass = 0.938272;
	Double_t Angle, Ep, Theta, xbj,  tail,RC,Q2;
	Int_t A, Z, i,j;

	Double_t XS_Born1, XS_Rad1,RadC1;
	Double_t XS_Born2, XS_Rad2,RadC2;
	Double_t XS_Born3, XS_Rad3,RadC3;

	T->Branch("E0", &E0, "E0/D");
	T->Branch("Ep", &Ep, "Ep/D");
	T->Branch("xbj", &xbj, "xbj/D");
	T->Branch("Q2", &Q2, "Q2/D");
	T->Branch("Angle", &Angle, "Angle/D");
	T->Branch("Theta", &Theta, "Theta/D");
	T->Branch("tail", &tail, "tail/D");
	T->Branch("RC", &RC, "RC/D");

	T->Branch("XS_Born1", &XS_Born1,"XS_Born1/D");
	T->Branch("XS_Rad1",  &XS_Rad1, "XS_Rad1/D");
	T->Branch("RadC1",    &RadC1,   "RadC1/D");

	T->Branch("XS_Born2", &XS_Born2,"XS_Born2/D");
	T->Branch("XS_Rad2",  &XS_Rad2, "XS_Rad2/D");
	T->Branch("RadC2",    &RadC2,   "RadC2/D");

	T->Branch("XS_Born3", &XS_Born3,"XS_Born3/D");
	T->Branch("XS_Rad3",  &XS_Rad3, "XS_Rad3/D");
	T->Branch("RadC3",    &RadC3,   "RadC3/D");

	T->Branch("A", &A, "A/I");
	T->Branch("Z", &Z, "Z/I");
	T->Branch("i", &i, "i/I");
	T->Branch("j", &j, "j/I");

	TString Com;
	infile11 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile21 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile31 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile11.eof())){
		infile11 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born1 >> XS_Rad1 >> tail;
		RadC1 = XS_Rad1 / XS_Born1;
		infile21 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born2 >> XS_Rad2 >> tail;
		RadC2 = XS_Rad2 / XS_Born2;
		infile31 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born3 >> XS_Rad3 >> tail;
		RadC3 = XS_Rad3 / XS_Born3;

		Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile11.close(); infile21.close();	infile31.close();

	infile12 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile22 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile32 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile12.eof())){
		infile12 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born1 >> XS_Rad1 >> tail;
		RadC1 = XS_Rad1 / XS_Born1;
		infile22 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born2 >> XS_Rad2 >> tail;
		RadC2 = XS_Rad2 / XS_Born2;
		infile32 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born3 >> XS_Rad3 >> tail;
		RadC3 = XS_Rad3 / XS_Born3;

		Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile12.close(); infile22.close();	infile32.close();

	infile13 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile23 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	infile33 >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com >> Com;
	while(!(infile13.eof())){
		infile13 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born1 >> XS_Rad1 >> tail;
		RadC1 = XS_Rad1 / XS_Born1;
		infile23 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born2 >> XS_Rad2 >> tail;
		RadC2 = XS_Rad2 / XS_Born2;
		infile33 >> A >> Z >> Angle >> Ep >> j>> Theta >> i >> XS_Born3 >> XS_Rad3 >> tail;
		RadC3 = XS_Rad3 / XS_Born3;

		Q2 = 4.0 * E0 *Ep*pow( sin(0.5*Theta*3.1415926/180) ,2);
		xbj = Q2 / (2.0 * Mass * (E0-Ep ));

		T->Fill();	
	}
	infile13.close(); infile23.close();	infile33.close();

	file->cd(); T->Write(); file->Close();

}
