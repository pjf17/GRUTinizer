#include "TH2D.h"

void norm2dhist(GH2D *h2, std::string xy = "y") {
    int nbinsY = h2->GetNbinsY();
    int nbinsX = h2->GetNbinsX();
    TH2D *hclone = new TH2D("normalized","normalized",nbinsX,h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),nbinsY,h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax());

    bool isX = xy.find("x") != std::string::npos;
    bool isY = xy.find("y") != std::string::npos;
    
    if (isY && !isX){
        for (int by=1; by <= nbinsY; by++){
            double norm = h2->Integral(1,nbinsX,by,by);
            if (norm < 100) continue;
            for (int bx=1; bx <= nbinsX; bx++){
                hclone->SetBinContent(hclone->GetBin(bx,by),h2->GetBinContent(h2->GetBin(bx,by))*1.0/norm);
                hclone->SetBinError(hclone->GetBin(bx,by),h2->GetBinError(h2->GetBin(bx,by))*1.0/norm);
            }
        }
    } 
    else if (!isY && isX) { 
        for (int bx=1; bx <= nbinsX; bx++){
            double norm = h2->Integral(bx,bx,1,nbinsY);
            if (norm < 100) continue;
            for (int by=1; by <= nbinsY; by++){
                hclone->SetBinContent(hclone->GetBin(bx,by),h2->GetBinContent(h2->GetBin(bx,by))*1.0/norm);
                hclone->SetBinError(hclone->GetBin(bx,by),h2->GetBinError(h2->GetBin(bx,by))*1.0/norm);
            }
        }
    }
    else if (isY && isX) {
        for (int bx=1; bx <= nbinsX; bx++){
            double norm = h2->Integral(bx,bx,1,nbinsY);
            if (norm < 100) continue;
            for (int by=1; by <= nbinsY; by++){
                hclone->SetBinContent(hclone->GetBin(bx,by),h2->GetBinContent(h2->GetBin(bx,by))*1.0/norm);
                hclone->SetBinError(hclone->GetBin(bx,by),h2->GetBinError(h2->GetBin(bx,by))*1.0/norm);
            }
        }

        for (int by=1; by <= nbinsY; by++){
            double norm = h2->Integral(1,nbinsX,by,by);
            if (norm < 100) continue;
            for (int bx=1; bx <= nbinsX; bx++){
                hclone->SetBinContent(hclone->GetBin(bx,by),hclone->GetBinContent(hclone->GetBin(bx,by))*1.0/norm);
                hclone->SetBinError(hclone->GetBin(bx,by),hclone->GetBinError(hclone->GetBin(bx,by))*1.0/norm);
            }
        }
    }
    GH2D *hgrut = (GH2D*) hclone;
    hgrut->Draw("colz");
}