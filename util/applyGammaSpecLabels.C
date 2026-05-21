void applyGammaSpecLabels() {
    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }

    std::string xtitle = "Energy [keV]";
    std::string ytitle = "Counts / keV";
    
    int bw = hh[0]->GetBinWidth(1);
    if (bw > 1) 
        ytitle = Form("Counts / %d keV",bw);

    int nHists = hh.size();
    for (int i=0; i < nHists; i++){
        hh[i]->GetXaxis()->SetTitle(xtitle.c_str());
        hh[i]->GetYaxis()->SetTitle(ytitle.c_str());

        hh[i]->GetXaxis()->SetLabelSize(0.04);
        hh[i]->GetYaxis()->SetLabelSize(0.04);

        hh[i]->GetXaxis()->SetTitleOffset(1.0);
        hh[i]->GetYaxis()->SetTitleOffset(1.1);

        hh[i]->GetXaxis()->SetTitleSize(0.045);
        hh[i]->GetYaxis()->SetTitleSize(0.045);
    }

    gPad->Update();
}

void rebin(int rbw){
    //get the histograms on the canvas
    TH1 *hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh = (TH1*)obj;
            hh->Rebin(rbw);
        }
    }
    gPad->Modified();
    gPad->Update();
    return;
}

void reSizeCanvas(double ww=800, double hh=500){
    gPad->GetCanvas()->SetWindowSize(ww,hh);
    gPad->Modified();
    gPad->Update();
}

void SetAxisrange(double xlo, double xhi, bool doX){
    if (xlo > xhi) std::swap(xlo,xhi);
    //get the histograms on the canvas
    TH1 *hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh = (TH1*)obj;
            if (doX) hh->GetXaxis()->SetRangeUser(xlo,xhi);
            else hh->GetYaxis()->SetRangeUser(xlo,xhi);
        }
    }
    gPad->Modified();
    gPad->Update();
    return;
}

void SetXrange(double xlo, double xhi){
    SetAxisrange(xlo,xhi,true);
}

void SetYrange(double xlo, double xhi){
    SetAxisrange(xlo,xhi,false);
}

void SetTitleOpt(bool isX, double size, double offset){
    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }

    int nHists = hh.size();
    for (int i=0; i < nHists; i++){
        if (size > 0.0){
            if (isX) hh[i]->GetXaxis()->SetTitleSize(size);
            else hh[i]->GetYaxis()->SetTitleSize(size);
        }

        if (offset > 0.0){
            if (isX) hh[i]->GetXaxis()->SetTitleOffset(offset);
            else hh[i]->GetYaxis()->SetTitleOffset(offset);
        }
    }

    gPad->Modified();
    gPad->Update();
}

void SetLabelSize(bool isX, double size){
    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }

    int nHists = hh.size();
    for (int i=0; i < nHists; i++){
        if (size > 0.0){
            if (isX) hh[i]->GetXaxis()->SetLabelSize(size);
            else hh[i]->GetYaxis()->SetLabelSize(size);
        }
    }

    gPad->Modified();
    gPad->Update();
}

void SetTitleXsize(double size=0.045) {SetTitleOpt(true,size,0.0);} 
void SetTitleXoffset(double offset=1.0) {SetTitleOpt(true,0.0,offset);} 
void SetLabelXsize(double size=0.04) {SetLabelSize(true,size);}

void SetTitleYsize(double size=0.045) {SetTitleOpt(false,size,0.0);} 
void SetTitleYoffset(double offset=1.1) {SetTitleOpt(false,0.0,offset);}
void SetLabelYsize(double size=0.04) {SetLabelSize(false,size);} 

void quickFormat(int wx=800, int wy=500, std::string xtitle ="", std::string ytitle ="") {
    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }
    gPad->GetCanvas()->SetWindowSize(wx,wy);

    //montecarlo res vs beta settings
    // gPad->SetTopMargin(0.08);
    // gPad->SetBottomMargin(0.11);
    // gPad->SetLeftMargin(0.11);
    // gPad->SetRightMargin(0.08);

    //normal settings
    double max = hh[0]->GetMaximum();

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.11);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.025);
    

    //if titles are blank do default spec labels
    if (xtitle == "") xtitle = "Energy [keV]";
    if (ytitle == "") {
        ytitle = "Counts / keV";
        if ( hh[0]->GetBinWidth(1) > 1)
            ytitle = std::string(Form("Counts %d / keV",(int) hh[0]->GetBinWidth(1)));
    }

    int nHists = hh.size();
    for (int i=0; i < nHists; i++){
        hh[i]->GetXaxis()->SetTitle(xtitle.c_str());
        hh[i]->GetXaxis()->SetTitleFont(43);
        hh[i]->GetXaxis()->SetTitleSize(25);
        hh[i]->GetXaxis()->SetTitleOffset(0.9);
        hh[i]->GetXaxis()->SetTitleColor(kBlack);
        hh[i]->GetXaxis()->SetLabelFont(43);
        hh[i]->GetXaxis()->SetLabelSize(22);
        hh[i]->GetXaxis()->SetLabelOffset(0.01);

        hh[i]->GetYaxis()->SetTitle(ytitle.c_str());
        hh[i]->GetYaxis()->SetTitleFont(43);
        hh[i]->GetYaxis()->SetTitleSize(25);
        hh[i]->GetYaxis()->SetTitleOffset(0.95);
        hh[i]->GetYaxis()->SetLabelFont(43);
        hh[i]->GetYaxis()->SetLabelSize(22);
        hh[i]->GetYaxis()->SetLabelOffset(0.01);
    }

    gPad->Modified();
    gPad->Update();
}

void ggQuickFormat() {
    // int wx=800;
    // int wy=400;
    int wx=600;
    int wy=400;
    
    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }
    gPad->GetCanvas()->SetWindowSize(wx,wy);
    gPad->SetTopMargin(0.0175);
    gPad->SetBottomMargin(0.145);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.04);
    // gPad->SetLeftMargin(0.1);
    // gPad->SetRightMargin(0.035);
    double max = hh[0]->GetMaximum();

    //if titles are blank do default spec labels
    std::string xtitle = "Energy [keV]";
    std::string ytitle = "Counts / keV";
    if ( hh[0]->GetBinWidth(1) > 1)
        ytitle = std::string(Form("Counts %d / keV",(int) hh[0]->GetBinWidth(1)));

    int nHists = hh.size();
    for (int i=0; i < nHists; i++){
        hh[i]->GetXaxis()->SetTitle(xtitle.c_str());
        hh[i]->GetXaxis()->SetTitleFont(43);
        hh[i]->GetXaxis()->SetTitleSize(25);
        hh[i]->GetXaxis()->SetTitleOffset(0.9);
        hh[i]->GetXaxis()->SetLabelFont(43);
        hh[i]->GetXaxis()->SetLabelSize(22);
        hh[i]->GetXaxis()->SetLabelOffset(0.01);

        hh[i]->GetYaxis()->SetTitle(ytitle.c_str());
        hh[i]->GetYaxis()->SetTitleFont(43);
        hh[i]->GetYaxis()->SetTitleSize(25);
        hh[i]->GetYaxis()->SetTitleOffset(0.7);
        hh[i]->GetYaxis()->SetLabelFont(43);
        hh[i]->GetYaxis()->SetLabelSize(22);
        hh[i]->GetYaxis()->SetLabelOffset(0.01);
    }

    gPad->Modified();
    gPad->Update();
}

void readEnergies(std::string filename, std::vector<std::pair<int,bool>> &energies){
    std::ifstream inFile(filename);
    std::string line;
    while (std::getline(inFile,line)){
        std::stringstream ss = std::stringstream(line);
        int eng;
        bool scheme;
        ss >> eng >> scheme;
        energies.push_back(std::make_pair(eng,scheme));
    }
    inFile.close();
    return;
}

void labelPeaks(int includeMode = 0, std::vector<int> highlight = {}, int includeHighlight = true, std::string filename="") {
    if (filename == "") filename = "peak_labels.txt";

    //get the histograms on the canvas
    std::vector<TH1 *> hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh.push_back( (TH1*)obj );
        }
    }
    //if multiple, get the histogram with the highest value
    TH1 *h = hh[0];
    if (hh.size() > 1){
        for (int i=1; i < hh.size(); i++){
            if (hh[i]->GetMaximum() > h->GetMaximum())
                h = hh[i];
        }
    }

    std::vector<std::pair<int,bool>> peaks;
    readEnergies(filename,peaks);
    std::sort(peaks.begin(), peaks.end());

    //x range
    double xmin = h->GetBinCenter(h->GetXaxis()->GetFirst());     
    double xmax = h->GetBinCenter(h->GetXaxis()->GetLast());
    double xRange = (xmax-xmin);
    //y range
    double ymax = h->GetMaximum();
    double rescaleYaxis = ymax*1.25;
    double labelHeight = ymax*0.85;
    // h->GetYaxis()->SetRangeUser(0,rescaleYaxis);

    //find peak index range and the gated peak
    int nPeaks = peaks.size();
    int idxLo = 0;
    int idxHi = nPeaks-1;
    for (int i=0; i < nPeaks-1; i++){
        if ( (peaks[i].first*1.0 - xmax)*(peaks[i+1].first*1.0 - xmax) < 0 )
            idxHi = i;
        if ( (peaks[i].first*1.0 - xmin)*(peaks[i+1].first*1.0 - xmin) < 0 )
            idxLo = i+1;
    }

    //figure out which peak is being gated on and don't include it if this is gamma-gamma
    std::string hname(h->GetName());
    int gatedPeak = -1;
    bool is_gammagamma = hname.find("gamma_gamma") != std::string::npos;  
    if (is_gammagamma) {
        std::vector <std::string> tokens;
        stringstream check1(hname);
        string intermediate;

        while(getline(check1, intermediate, '_')) tokens.push_back(intermediate);
        int binLo = std::stoi(tokens[tokens.size()-2]);
        int binHi = std::stoi(tokens[tokens.size()-1]);
        gatedPeak = (h->GetBinCenter(binLo/2) + h->GetBinCenter(binHi/2))*1.0/2;
        int gpi = -1;
        double smallestDiff = 1000;

        for (int i=0; i < nPeaks-1; i++) {
            if (std::abs(peaks[i].first - gatedPeak) < smallestDiff ){
                smallestDiff = std::abs(peaks[i].first - gatedPeak);
                gpi = i;
            }
        }
        gatedPeak = peaks[gpi].first;
    }

    std::vector<TLine *> peakLines;
    for (int i=idxLo; i < idxHi+1; i++) {
        bool isInHighlight = false;
        if (highlight.size() > 0) {
            for (auto e: highlight){
                if (e == peaks[i].first){
                    isInHighlight = true;
                    break;
                } 
            }
        } 
        if (highlight.size() > 0 && (isInHighlight != includeHighlight)) continue;
        
        int drawColor = kBlack;

        switch (includeMode){
            case -1: //draw only unplaced
                if (peaks[i].second == 1) continue; 
                break;
            case 0: //draw only placed
                if (peaks[i].second == 0) continue; 
                break;
            case 1: //draw everything, draw unplaced with blue
                if (peaks[i].second == 0) drawColor = kBlue;
                break;
            case 2: //draw everything the same
                break;
            default:
                std::cout<<"BAD DRAW MODE\n";
                break;
        }

        //get the total x dimensions
        double lfm = gPad->GetLeftMargin();
        double rtm = gPad->GetRightMargin();
        double xAxis = 1.0 - lfm - rtm;
        double canvWidthKeV = xRange*(1.0 + rtm/xAxis + lfm/xAxis);
        bool gate_clip = peaks[i].first/xmax > 0.9 && is_gammagamma;
        double clip_factor = gate_clip ? 0.8 : 1.0;
        
        int peakBin = h->FindBin(peaks[i].first);
        double peakHeight = std::max(h->GetBinContent(peakBin-1),h->GetBinContent(peakBin));
        peakHeight = std::max(peakHeight,h->GetBinContent(peakBin+1));
        double peakLength = 0.05*rescaleYaxis;
        double peakYmax = peakHeight+0.02*rescaleYaxis + peakLength;

        bool isShortPeak = true;
        double Xshift = 0.0;
        if (std::abs(peakHeight - ymax)/ymax < 0.05) {
            isShortPeak = false;
            Xshift = -0.02*canvWidthKeV;
        }

        // peakLines.push_back( new TLine(peaks[i].first,peakHeight+0.02*rescaleYaxis,peaks[i].first,peakYmax) );
        peakLines.push_back( new TLine(peaks[i].first+Xshift,peakHeight+0.02*ymax,peaks[i].first+Xshift,labelHeight*clip_factor) );
        peakLines.back()->SetLineStyle(kDashed);
        peakLines.back()->SetLineColor(drawColor);
        if (isShortPeak) peakLines.back()->Draw("same");

        // TText *tt = new TText(peaks[i].first,peakYmax,Form("%d",peaks[i].first));
        TText *tt = new TText(peaks[i].first,labelHeight*clip_factor,Form("%d",peaks[i].first));
        // double ttSize = 0.04; //normal size
        
        // check neighbors for close peaks
        bool narrowToLeft = i > idxLo && abs(peaks[i].first - peaks[i-1].first) < 0.028*canvWidthKeV;
        bool narrowToRight = i < nPeaks-2 && abs(peaks[i].first - peaks[i+1].first) < 0.038*canvWidthKeV;
        if (narrowToLeft && !narrowToRight) tt->SetX(peaks[i].first + 0.016*canvWidthKeV);
        
        // if (narrowToLeft && !narrowToRight){
        //     if (narrowToRight) ttSize = abs(peaks[i].first - peaks[i+1])/canvWidthKeV + 0.05;
        //     else tt->SetX(peaks[i].first + 0.02*canvWidthKeV);
        // }
        // if (ttSize > 0.035) ttSize = 0.035;
        
        tt->SetTextAngle(90);
        tt->SetTextFont(43);
        tt->SetTextSize(20);
        tt->SetTextColor(drawColor);
        tt->Draw();
    }

    if (is_gammagamma) {
        TLatex tgated;
        tgated.SetTextFont(43);
        tgated.SetTextSize(22);
        tgated.SetTextAlign(33);
        tgated.DrawLatexNDC(0.93,0.96,Form("#gamma_{gate} %d",gatedPeak));
        gPad->Update();
    }
}

void drawSubFig(std::string letter, double xx=0.15, double yy=0.88){
    TText text;
    text.SetTextAlign(31);
    text.SetTextFont(43);
    text.SetTextSize(25);
    text.DrawTextNDC(xx,yy,letter.c_str());
}