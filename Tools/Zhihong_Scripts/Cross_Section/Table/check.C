{
	int runno;
	double mom_epics, mom_set, mom_calc, err;
	TString com;

	ifstream inf1("E08014_Mom_Right.table");
	ifstream inf2("E08014_Mom_Left.table");

	TH2F *h1 = new TH2F("h1","Central Momentum Deviation on HRS-R", 500, 3550, 4310,500, -0.1,0.1);
	TH2F *h2 = new TH2F("h2","Central Momentum Deviation on HRS-L", 500, 3550, 4310,500, -0.1,0.1);

	inf1 >> com >> com >> com >> com >> com;
	while(!(inf1.eof())){
		inf1 >> runno >> mom_calc >> mom_set >> com >> com;
		err = (mom_set - mom_calc) / mom_set * 100;
		h1->Fill(runno, err);
	}
	inf1.close();
	inf2 >> com >> com >> com >> com >> com;
	while(!(inf2.eof())){
		inf2 >> runno >> mom_calc >> mom_set >> com >> com;
		err = (mom_set - mom_calc) / mom_set * 100;
		h2->Fill(runno, err);
	}
	inf2.close();

	TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
    h1->SetMarkerColor(2);
	h1->SetMarkerStyle(21);
	h1->Draw();
    h2->SetMarkerColor(4);
	h2->SetMarkerStyle(23);
	h2->Draw("same");
}
