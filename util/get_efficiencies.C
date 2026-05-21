void get_efficiencies(std::vector<double> peaks){
    int Npeaks = peaks.size();
    int TotEvents = 8000000;

    std::vector<double> eff;
    for (int i=0; i < Npeaks; i++){
        TFile *ffile = new TFile(Form("hist%d.root",(int) peaks[i]),"READ");
        TH1 *hh = (TH1*) ffile->Get("gretsim/CoreEnergy_FEP");
        eff.push_back(1.0*hh->GetEntries()/TotEvents);
        ffile->Close();
    }

    TF1 *f1 = new TF1("f1","[0]*43.0/32.0 * 4.532 * pow(100 + x,-0.621)",500,1100);

    TGraph *gg = new TGraph(Npeaks, &peaks[0], &eff[0]);

    gg->Draw("A*");
    gg->Fit(f1,"R");
    std::cout<<1.0/f1->GetParameter(0)<<std::endl;
    // f1->Draw("same");
}

void get_efficiencies(std::string listfile){
    int TotEvents = 8000000;

    std::ifstream file(listfile);

    std::string line;

    std::vector<double> eff;
    std::vector<double> peaks;
    while (std::getline(file, line)) {
        //get peak energy
        std::stringstream sline = std::stringstream(line);
        int pk;
        sline >> pk;

        //open hist
        TFile *ffile = new TFile(Form("hist%d.root",pk),"READ");
        TH1 *hh = (TH1*) ffile->Get("gretsim/dopEn_fep");

        double raw_eff = 1.0*hh->GetEntries()/TotEvents;

        printf("%d %f\n",pk,raw_eff*1.2794);

        eff.push_back(raw_eff*1.2794);
        peaks.push_back(pk);

        ffile->Close();
    }

    TF1 *f1 = new TF1("f1","43.0/32.0 * 4.532 * pow(100 + x,-0.621)",200,7000);

    TGraph *gr = new TGraph((int) eff.size(), &peaks[0], &eff[0]);
    gr->Draw("A*");
    f1->Draw("same");
}