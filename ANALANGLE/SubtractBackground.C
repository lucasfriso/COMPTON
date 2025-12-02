#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include "TString.h"
#include <iostream>

void SubtractBackground() {

    TString fileDataName = "DataM_Aluminum100.root";
    TString fileBkgName  = "DataM_Background100.root";
    TString treeName     = "Data_B";


    Double_t calib_Tag_p0 = -13.5147;
    Double_t calib_Tag_p1 = 2.351;


    Double_t calib_Det_p0 = -11.0774;
    Double_t calib_Det_p1 = 2.412935;


    Double_t calib_Sct_p0 = -14.25;
    Double_t calib_Sct_p1 = 2.85775;


    Double_t timeData = 2000;
    Double_t timeBkg  = 1000;


    TCut fold2Cut = "Channel_0.Energy > 10 && Channel_2.Energy > 10";


    TString varFormula;
    varFormula.Form("(%f + %f * Channel_2.Energy)", calib_Det_p0, calib_Det_p1);

    std::cout << "Plotting Variable: " << varFormula << std::endl;


    Double_t scaleFactor = timeData / timeBkg;


    TFile *fData = new TFile(fileDataName);
    TTree *tData = (TTree*)fData->Get(treeName);
    TFile *fBkg = new TFile(fileBkgName);
    TTree *tBkg = (TTree*)fBkg->Get(treeName);

    if (!tData || !tBkg) {
        std::cout << "Error opening files or trees." << std::endl;
        return;
    }


    TH1F *hData = new TH1F("hData", "Calibrated Energy (Al);Energy (keV);Counts", 100, 0, 500);
    TH1F *hBkg  = new TH1F("hBkg",  "Background Run;Energy (keV);Counts", 100, 0, 500);

    hData->Sumw2();
    hBkg->Sumw2();


    tData->Draw(varFormula + ">>hData", fold2Cut);
    tBkg->Draw(varFormula + ">>hBkg", fold2Cut);


    TH1F *hSignal = (TH1F*)hData->Clone("hSignal");
    hSignal->SetTitle("Background Subtracted Signal (Calibrated)");

    hBkg->Scale(scaleFactor);
    hSignal->Add(hBkg, -1);


    TCanvas *c1 = new TCanvas("c1", "Calibrated Spectra", 900, 600);
    hData->SetLineColor(kBlue);
    hBkg->SetLineColor(kRed); hBkg->SetLineStyle(2);
    hSignal->SetLineColor(kGreen+2); hSignal->SetFillColorAlpha(kGreen+2, 0.3);
    Double_t binWidth = hSignal->GetBinWidth(1);

    hData->Draw("HIST E");
    hBkg->Draw("HIST SAME");
    hSignal->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg->AddEntry(hData, "Data (Al)", "l");
    leg->AddEntry(hBkg, "Scaled Bkg", "l");
    leg->AddEntry(hSignal, "Signal", "f");
    leg->Draw();
    std::cout << "Bin Width: " << binWidth << " keV" << std::endl; //
}