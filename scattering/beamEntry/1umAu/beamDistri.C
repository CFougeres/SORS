#include <iostream>
void beamDistri(){
	int Nion = 9999;
	TH1F* theta_distri = new TH1F("theta_distri_entrance","Au 1#mum window; #theta_{beam} (deg);counts/0.1degree",500,0,5);
	double theta,x,y,z, i; double Pi = TMath::Pi();
	Double_t X[Nion], Y[Nion], Z[Nion];
	std::ifstream in;
	in.open("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/scattering/beamEntry/1umAu/range.dat");
	Int_t nlines = 0;
	while (1) {
		in >> i >> z >> x >> y;
		Z[nlines] = z / 10000000.; X[nlines] = x / 10000000.;  Y[nlines] = y / 10000000.;
		if (!in.good()) break;
		theta = acos(Z[nlines] / sqrt(X[nlines] * X[nlines] + Y[nlines] * Y[nlines] + Z[nlines] * Z[nlines])) * 180. / Pi;
		theta_distri->Fill(theta);
		nlines++;
	}
	printf(" found %d points\n", nlines);
	in.close();
	TGraph* g = new TGraph(Nion, X, Y);
	g->SetTitle("Au 1#mum window;X (mm); Y (mm))");
	TCanvas* c = new TCanvas("c", "", 800, 800); c->Divide(2, 1);
	c->cd(1);	g->Draw("A*");
	c->cd(2);	theta_distri->Draw();
}
