#include "TH2D.h"

void norm2dhist(GH2D *h2) {
    int nbinsY = h2->GetNbinsY();
    int nbinsX = h2->GetNbinsX();
    TH2D *hclone = new TH2D("normalized","normalized",nbinsX,h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),nbinsY,h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax());
    for (int by=1; by <= nbinsY; by++){
        double norm = h2->Integral(1,nbinsX,by,by);
        // std::cout<<norm<<std::endl;
        if (norm < 100) continue;
        for (int bx=1; bx <= nbinsX; bx++){
            hclone->SetBinContent(hclone->GetBin(bx,by),h2->GetBinContent(h2->GetBin(bx,by))*1.0/norm);
            hclone->SetBinError(hclone->GetBin(bx,by),h2->GetBinError(h2->GetBin(bx,by))*1.0/norm);
        }
    }
    GH2D *hgrut = (GH2D*) hclone;
    hgrut->Draw("colz");
}