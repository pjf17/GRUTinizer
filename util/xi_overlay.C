#include "MultiPlotter.h"

// MultiPlotter *mp = new MultiPlotter();
// 
// void xi_overlay(TFile *f1, TFile *f2){
    // 
    // mp->Add(f1,"Bkg Subtracted");
    // mp->Add(f2,"Bkg Subtracted");
// 
    // mp->SetLegendEntry("Bkg Subtracted","#tilde{#xi} < 30#circ");
    // mp->SetLegendEntry("Bkg Subtracted_perp","#tilde{#xi} > 60#circ");
// 
    // mp->SetLineColor("Bkg Subtracted_perp",kBlue);
    // mp->SetLineColor("Bkg Subtracted",kRed);
    // mp->SetFill(0.4);
    // mp->SetLegendCoordinates(0.78,0.84,0.98,0.98);
    // mp->Draw();
// }

void do_format(TH1* hh){
    hh->GetXaxis()->SetTitle("Energy [keV]");
    if (hh->GetBinWidth(1) > 1)
        hh->GetYaxis()->SetTitle(Form("Counts / %d keV",int(hh->GetBinWidth(1))));
    else 
        hh->GetYaxis()->SetTitle("Counts / keV");

    hh->GetXaxis()->SetLabelSize(0.055);
    hh->GetYaxis()->SetLabelSize(0.055);

    hh->GetXaxis()->SetTitleSize(0.06);
    hh->GetYaxis()->SetTitleSize(0.06);

    hh->GetXaxis()->SetTitleOffset(0.9);
    hh->GetYaxis()->SetTitleOffset(1.55);

    hh->GetXaxis()->SetNdivisions(705);
}

void xi_overlay(TFile *f1, TFile *f2){
    TH1* para = (TH1*) f1->Get("Bkg Subtracted");
    TH1* perp = (TH1*) f2->Get("Bkg Subtracted");
    TLegend *leg = new TLegend(0.785,0.77,0.97,0.98);
    leg->SetTextAlign(22);
    leg->AddEntry(perp,"#tilde{#xi} > 60#circ","l");
    leg->AddEntry(para,"#tilde{#xi} < 30#circ","l");

    para->SetLineColor(kRed);
    para->SetFillColorAlpha(kRed,0.4);

    perp->SetLineColor(kBlue);
    perp->SetFillColorAlpha(kBlue,0.4);

    do_format(perp);
    do_format(para);

    TCanvas *canv = new TCanvas("polarization","polarization",500,450);

    perp->Draw("hist");
    para->Draw("hist same");

    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.02);

    canv->Update();

    leg->Draw("same");
}

void drawSubFig(std::string str){
    TText text;
    text.SetTextAlign(12);
    text.SetTextFont(43);
    text.SetTextSize(26);
    text.DrawTextNDC(0.22,0.92,str.c_str());

    // TLatex eng;
    // eng.SetTextAlign(12);
    // eng.SetTextFont(43);
    // eng.SetTextSize(26);
    // eng.DrawLatexNDC(0.22,0.82,Form("E_{#gamma} = %d",energy));
}