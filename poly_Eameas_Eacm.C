//***************************************
//Simulation resonant scattering reactions in gaseous target
//+ Si detection at 0 degree
//***************************************

#include <iostream>
#include <stdlib.h>

//***************************************
//Extraction function for SRIM Stop. Pow.
//***************************************
int extraction_data_SRIM(Char_t* file, Double_t EnergyPart[], Double_t RangePart_in_Au[], Double_t SPelPart_in_Au[], Double_t SPnuPart_in_Au[]) {
    std::ifstream in;
    in.open(file);
    Int_t nlines = 0;
    while (1) {
        in >> EnergyPart[nlines] >> SPelPart_in_Au[nlines] >> SPnuPart_in_Au[nlines] >> RangePart_in_Au[nlines]; // SP in keV/microns, E in keV, Range in Angstrom
        if (!in.good()) break;
        nlines++;
    }
    in.close();
    return 1;
}


//***************************************
//Interpolation on SRIM Stop. Pow.
//***************************************
double interpol(double x[], double y[], int const n, double p)
{
    int i = 1;
    double val, b1, c1, d1, c2, d2;
    // b1 false if something is wrong
    bool bb1 = true;
    // loop breaking
    bool bb2 = true;
    if (p <= x[1])
    {
        bb1 = false;
    }
    else
    {
        do
        {
            i++;
            if (p <= x[i])bb2 = false;
        } while (bb2 && (i < n - 2));
    }
    if (p > x[n - 2])bb1 = false;
    if (bb1 == true)
    {
        b1 = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);
        b1 = b1 * p - b1 * x[i - 2] + y[i - 2];
        c1 = (y[i] - y[i - 2]) / (x[i] - x[i - 2]);
        c1 = c1 * p - c1 * x[i - 2] + y[i - 2];
        d1 = (y[i + 1] - y[i - 2]) / (x[i + 1] - x[i - 2]);
        d1 = d1 * p - d1 * x[i - 2] + y[i - 2];
        c2 = (c1 - b1) / (x[i] - x[i - 1]);
        c2 = c2 * p - c2 * x[i - 1] + b1;
        d2 = (d1 - b1) / (x[i + 1] - x[i - 1]);
        d2 = d2 * p - d2 * x[i - 1] + b1;
        val = (d2 - c2) / (x[i] - x[i - 1]);
        val = val * p - val * x[i] + c2;
    }
    else
    {
        val = 0.;
    }
    return val;
}


//***************************************
// Local derivation of energy losses
// energy in MeV, distance in mm,
// Sto. Pow. in keV/microns
//***************************************
double loss_E_srim(double Ei, double distance, Double_t EnergyPart[], Double_t SPel_Part[], Double_t SPnu_Part[], int N_energy) {
    //distance in mm
    double Etemp = Ei;
    double dpart = 0;
    double ddx = 5.0; //micron
    while (dpart <= distance*1000) {
        Etemp = Etemp - 0.001 * interpol(EnergyPart, SPel_Part, N_energy, Etemp * 1000.0) * ddx - 0.001 * interpol(EnergyPart, SPnu_Part, N_energy, Etemp * 1000.0) * ddx;
        dpart = dpart + ddx;
    }
    return Ei - Etemp;
}
double loss_E_srim_wind(double Ei, double distance, Double_t EnergyPart[], Double_t SPel_Part[], Double_t SPnu_Part[], int N_energy) {
    //distance in mm
    double Etemp = Ei;
    double dpart = 0;
    double ddx = 0.1; //micron
    while (dpart <= distance*1000) {
        Etemp = Etemp - 0.001 * interpol(EnergyPart, SPel_Part, N_energy, Etemp * 1000.0) * ddx - 0.001 * interpol(EnergyPart, SPnu_Part, N_energy, Etemp * 1000.0) * ddx;
        dpart = dpart + ddx;
    }
    return Ei - Etemp;
}

//***************************************
//Derivation scattered beam after entrance window
//***************************************
double scattered_theta() {
    TFile* file1 = new TFile("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/scattering/beamEntry/theta_distri_entrance.root", "READ");
    TH1F* srim_scattered = (TH1F*)file1->Get("theta_distri_entrance");
    gRandom->SetSeed(0);
    double theta = srim_scattered->GetRandom();
    file1->Close();
    return theta;
}

//***************************************
// Main Function
//***************************************
void poly_Eameas_Eacm(){
    
    //***************************************
    //Experiment conditions
    //***************************************
    int choice_pressure_target = 3;
    double EbeamEntry = 1.8*15.001138;// MeV
    int aivalable_pressure_target[4]={92,183,274,350};
    double target_width = 90;//mm Al==85  Au105.
    /*Al windows*/
  /*  double exitWindow_width = 0.020; //mm
    double entranceWindow_width=0.004;//mm*/
    /*Au windows*/
    double exitWindow_width = 0.005; //mm
    double entranceWindow_width=0.0048;//mm
    double angle_max_Si = 2;//degree
    //Si detector response
    TRandom* rEloss = new TRandom();//gaussian reponse of Si detector
    double sigma_e = 0.001*10.0;//MeV
    srand(time(NULL));
    gRandom->SetSeed(0);
    
    //***************************************
    //Graph Ea_meas = f(Ea_cm) to fit
    //***************************************
    TCanvas* c = new TCanvas("c","c",800,800);
    /*Al*/
 //   double Ea_meas[1*1215-1];  double Ea_cm[1*1215-1];
    /*Au*/
  //  double Ea_meas[1*1500-1];  double Ea_cm[1*1500-1];
    /*Ti*/
    double Ea_meas[1*1285-1];  double Ea_cm[1*1285-1];
    
    
    //***************************************
    //SRIM Sto. Pow.
    //***************************************
    int N_energy = 125;
    int extraction_data;
    Double_t EnergyRecoil[2][125];     Double_t RangeRecoil_in_4He[2][125];    Double_t SPelRecoil_in_4He[2][125];    Double_t SPnuRecoil_in_4He[2][125];
    Double_t EnergyRecoil_window[2][125];     Double_t RangeRecoil_in_4He_window[2][125];    Double_t SPelRecoil_in_4He_window[2][125];    Double_t SPnuRecoil_in_4He_window[2][125];
    //15O
    /*  Char_t* fileSRIM[2] ={Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/4He_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target]),Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15O_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target])};
   std::cout<<Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15O_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target])<<std::endl;
   
 //  Char_t* fileSRIM_window[2] ={"/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/4He_27Al.dat","/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15O_27Al.dat"};
   Char_t* fileSRIM_window[2] ={"/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/4He_Ti.dat","/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15O_Ti.dat"};
   */
    
    //15N
    Char_t* fileSRIM[2] ={Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/4He_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target]),Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15N_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target])};
    std::cout<<Form("/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15N_4He_%iTorr.dat", aivalable_pressure_target[choice_pressure_target])<<std::endl;
    Char_t* fileSRIM_window[2] ={"/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/4He_Ti.dat","/Users/chloefougeres/Documents/experiment/experiments_proposed/4He_15O_4He_15O/analysis/sp_srim/15N_Ti.dat"};

 for(int i=0;i<2;i++){
        extraction_data = extraction_data_SRIM(fileSRIM[i], EnergyRecoil[i], RangeRecoil_in_4He[i],SPelRecoil_in_4He[i], SPnuRecoil_in_4He[i]);
        extraction_data = extraction_data_SRIM(fileSRIM_window[i], EnergyRecoil_window[i], RangeRecoil_in_4He_window[i],SPelRecoil_in_4He_window[i], SPnuRecoil_in_4He_window[i]);
    }
    
    //***************************************
    //Reaction initialization
    //***************************************
    double mu = 931.5;
//    double m15O = 15 * mu +2.856 ;
    double m15N = 15 * mu +0.1014;
    double m4He = 4 * mu + 2.425;
    double mu4He=m4He /(m4He + m15N);
    double mu15N=m15N /(m4He + m15N);
    gSystem->Load("libPhysics");
    double Pfais = sqrt(EbeamEntry * EbeamEntry + 2. * EbeamEntry * m15N);
    TLorentzVector target(0.0, 0.0, 0.0, m4He * 0.001);
    TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (EbeamEntry + m15N) * 0.001);
    Double_t masses[2] ={0.001 * (m15N), 0.001 * m4He};
    TLorentzVector Tr = beam + target;
    TLorentzVector * pRecoil;
    TLorentzVector * pEjectil;
    TGenPhaseSpace event;
    event.SetDecay(Tr, 2, masses);
    //***************************************
    //Variables MC simu
    //***************************************
    int Nreac =500000;//500000
    int Nbeam =5000;//5000
    double dx=0.07;//mm
    int Nentry = (target_width/dx)*double(Nbeam)-1;
    double Ebeam, thetaBeam, beam_target_width, vbeamlab, ELossBeam, ElossRecoil;
    Ebeam = EbeamEntry;
    double thetaEjectil_lab, vEjectil_lab, eEjectil_lab, ElossEjectil, distEjectil, remainingDist, vEjectil_lab_meas, Ecm_reaction;
    double Emeas; int compteur=0; int compteur_reaction=0;
    double positionBeam=0;

    
    //********************************
    //Derivation reactions
    //********************************
    double etemplossbeam=0;
    for(int b=0;b<Nbeam;b++){
        positionBeam= 0 ;    Ebeam = EbeamEntry;
        thetaBeam = scattered_theta();
        beam_target_width =  target_width/cos(3.14*(thetaBeam/180.));
        ELossBeam = loss_E_srim_wind(Ebeam, entranceWindow_width/cos(3.14*(thetaBeam/180.)), EnergyRecoil_window[1], SPelRecoil_in_4He_window[1], SPnuRecoil_in_4He_window[1], N_energy);
        Ebeam = Ebeam - ELossBeam;
       std::cout<< Ebeam/15.001138<<std::endl;
    while(positionBeam<(beam_target_width-dx)){
        compteur_reaction=0;
        Ecm_reaction=mu4He*Ebeam;
       // dx = derivation_eff_target_width(Ecm_reaction, bin_to_E, fileSRIM[1]);
        positionBeam = positionBeam + dx;
        ELossBeam = loss_E_srim(Ebeam, dx, EnergyRecoil[1], SPelRecoil_in_4He[1], SPnuRecoil_in_4He[1], N_energy);
        Ebeam = Ebeam - ELossBeam;
        vbeamlab = sqrt(1 - 1 / pow((Ebeam / m15N + 1), 2));
        Pfais = sqrt(Ebeam * Ebeam + 2. * Ebeam * m15N);
        beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Ebeam+ m15N) * 0.001);
        Tr = beam + target;
        event.SetDecay(Tr, 2, masses);
        Ecm_reaction = (Ecm_reaction+mu4He*Ebeam)/2.;
        if(etemplossbeam-ELossBeam>0){
            std::cout<< positionBeam<<" Bragg peak "<<std::endl;
        }
        etemplossbeam=ELossBeam;
        for(int j=0;j<Nreac;j++){
            event.Generate();
            pRecoil = event.GetDecay(0); pEjectil = event.GetDecay(1);
            thetaEjectil_lab = 180./3.14*pEjectil->Theta();//Rad
            vEjectil_lab =pEjectil->Beta();
            eEjectil_lab =(1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * m4He;
            remainingDist = beam_target_width-positionBeam;
            if(eEjectil_lab>0 && thetaEjectil_lab<angle_max_Si && compteur_reaction==0){
                compteur_reaction=1;
                distEjectil = (remainingDist)/cos(3.14*thetaEjectil_lab/180.);
                ElossRecoil = loss_E_srim(eEjectil_lab, distEjectil , EnergyRecoil[0],SPelRecoil_in_4He[0], SPnuRecoil_in_4He[0], N_energy);
                eEjectil_lab = eEjectil_lab-ElossRecoil;
                ElossRecoil = loss_E_srim_wind(eEjectil_lab, exitWindow_width/cos(3.14*thetaEjectil_lab/180.), EnergyRecoil_window[0],SPelRecoil_in_4He_window[0], SPnuRecoil_in_4He_window[0], N_energy);
                Emeas= rEloss->Gaus(eEjectil_lab-ElossRecoil,sigma_e);
               if(compteur<Nentry){
                    Ea_meas[compteur]=Emeas; Ea_cm[compteur]=Ecm_reaction;
                }
                j=Nreac;
            }
        }
        compteur=compteur+1;
    }
        std::cout<<"beam after entrance "<<  beam_target_width <<std::endl;
    }
    std::cout<<compteur<< " "<<Nentry<<std::endl;
    std::cout<<"beam at the end "<<  Ebeam <<std::endl;
    
    //***************************************
    //Plotting
    //***************************************
    TGraph* function_EaMeas_EaCM = new TGraph(Nentry, Ea_meas, Ea_cm);
    std::cout<<Ea_cm[Nentry-1]<< " "<<Ea_meas[Nentry-1]<<std::endl;
    function_EaMeas_EaCM->SetTitle(";E_{#alpha}^{lab} at Si detector (MeV); E_{#alpha}^{cm} at elastic reaction (MeV)");
    c->cd();function_EaMeas_EaCM->Draw("A*");
}
