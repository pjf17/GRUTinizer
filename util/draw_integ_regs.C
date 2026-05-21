void readFile(std::string filename, std::map<int,std::pair<double,double>> &gmap){
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int first;
        double glo, ghi;
        std::string temp;

        // Read first number
        ss >> first;

        // Skip 7 fields
        for (int i = 0; i < 7; ++i) {
            ss >> temp;
        }

        // Read second-to-last and last number
        ss >> glo >> ghi;

        gmap[first] = std::make_pair(glo,ghi);
    }
}

std::vector<int> get_map_keys(const std::map<int,std::pair<double,double>> & my_map) {
    std::vector<int> keys;
    for (const auto& pair : my_map) {
        keys.push_back(pair.first);
    }
    return keys;
}

void draw_integ_regs(std::string fpara, std::string fperp, std::vector<int> keys = {}){
    //get the scale of the histogram
    TH1D *h;
    double ymax = 0;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            h = (TH1D*)obj;
            double tempy = h->GetMaximum();
            if (tempy > ymax) ymax = tempy;
        }
    }
    //read the file
    std::map<int,std::pair<double,double>> para, perp;
    readFile(fpara,para);
    readFile(fperp,perp);

    //draw all peaks if none specified
    if (keys.size() == 0) keys = get_map_keys(para);

    //draw the gates
    int Np = keys.size();
    for (int i = 0; i < Np; i++){
        int peak = keys[i];
        std::vector<TLine *> lpara = {new TLine(para[peak].first,0,para[peak].first,ymax), new TLine(para[peak].second,0,para[peak].second,ymax)};
        std::vector<TLine *> lperp = {new TLine(perp[peak].first,0,perp[peak].first,ymax), new TLine(perp[peak].second,0,perp[peak].second,ymax)};
        
        for (int j=0; j < 2; j++){
            lpara[j]->SetLineColor(kRed);
            lperp[j]->SetLineColor(kBlue);
            lpara[j]->SetLineStyle(kDashed);
            lperp[j]->SetLineStyle(kDashed);
            lperp[j]->SetLineWidth(3);
            lpara[j]->SetLineWidth(3);
            lpara[j]->Draw("same");
            lperp[j]->Draw("same");
        }
    }

    return;
}

