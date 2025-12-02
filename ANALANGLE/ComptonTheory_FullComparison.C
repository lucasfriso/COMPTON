#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TLatex.h>


Double_t calib_Tag_p0 = -13.5147;   Double_t calib_Tag_p1 = 2.351;
Double_t calib_Det_p0 = -11.0774;   Double_t calib_Det_p1 = 2.412935;
Double_t calib_Sct_p0 = -14.25;     Double_t calib_Sct_p1 = 2.8577;

Double_t tag_min = 470.0; Double_t tag_max = 540.0; Double_t scat_th = 20.0;
const double E0 = 511.0;

// --- THEORY ---
Double_t TheoryPhoton(Double_t *x, Double_t *par) {
    Double_t theta_rad = x[0] * TMath::DegToRad();
    return E0 / (1 + (E0/511.0) * (1 - TMath::Cos(theta_rad)));
}
Double_t TheoryElectron(Double_t *x, Double_t *par) { return E0 - TheoryPhoton(x, par); }

struct Measurement {
    double angle;
    double e_gam, e_gam_err;
    double e_ele, e_ele_err;
    double e_sum, e_sum_err;
    bool valid;
};

Measurement AnalyzeRun(const char* filename, double angle) {
    Measurement m;
    m.angle = angle;
    m.e_gam = -100; m.e_ele = -100; m.e_sum = -100;
    m.valid = (angle >= 50);

    TFile *f = new TFile(filename);
    if (!f || f->IsZombie()) return m;

    TTree *t = (TTree*)f->Get("Data_B");
    TLeaf *l_tag = t->GetLeaf("Channel_0.Energy");
    TLeaf *l_sct = t->GetLeaf("Channel_1.Energy");
    TLeaf *l_det = t->GetLeaf("Channel_2.Energy");
    if (!l_tag) { l_tag = t->GetLeaf("ch0.Energy"); l_sct = t->GetLeaf("ch1.Energy"); l_det = t->GetLeaf("ch2.Energy"); }

    TH1F *h_det = new TH1F("temp_det", "", 100, 0, 1000);
    TH1F *h_sct = new TH1F("temp_sct", "", 100, 0, 1000);
    TH1F *h_sum = new TH1F("temp_sum", "", 100, 0, 1200);

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);
        Double_t v_tag = calib_Tag_p0 + calib_Tag_p1 * l_tag->GetValue();
        Double_t v_sct = calib_Sct_p0 + calib_Sct_p1 * l_sct->GetValue();
        Double_t v_det = calib_Det_p0 + calib_Det_p1 * l_det->GetValue();

        if (v_tag > tag_min && v_tag < tag_max) {
            if (angle < 50 || (v_sct > scat_th && v_det > 10)) {
                h_det->Fill(v_det); h_sct->Fill(v_sct); h_sum->Fill(v_det + v_sct);
            }
        }
    }

    h_det->GetXaxis()->SetRangeUser(50, 1000);
    TF1 *f1 = new TF1("f1", "gaus");
    double peak = h_det->GetBinCenter(h_det->GetMaximumBin());
    h_det->Fit(f1, "Q N", "", peak-60, peak+60);
    m.e_gam = f1->GetParameter(1); m.e_gam_err = f1->GetParError(1);

    if (m.valid) {
        h_sct->GetXaxis()->SetRangeUser(50, 1000);
        peak = h_sct->GetBinCenter(h_sct->GetMaximumBin());
        h_sct->Fit(f1, "Q N", "", peak-80, peak+80);
        m.e_ele = f1->GetParameter(1); m.e_ele_err = f1->GetParError(1);

        h_sum->GetXaxis()->SetRangeUser(400, 600);
        peak = h_sum->GetBinCenter(h_sum->GetMaximumBin());
        h_sum->Fit(f1, "Q N", "", peak-50, peak+50);
        m.e_sum = f1->GetParameter(1); m.e_sum_err = f1->GetParError(1);
    }
    f->Close(); return m;
}

void ComptonTheory_FullComparison() {
    const int nTotal = 5;
    double angles[nTotal] = {10, 30, 60, 90, 110};
    const char* files[nTotal] = {
        "DataM_10t1.root", "DataM_30t1.root",
        "DataM_60t1.root", "DataM_90t1.root", "DataM_110t1.root"
    };

    TGraphErrors *gGamV = new TGraphErrors();
    TGraphErrors *gGamG = new TGraphErrors();
    TGraphErrors *gEleV = new TGraphErrors();
    TGraphErrors *gSumV = new TGraphErrors();

    int cV=0, cG=0;
    for(int i=0; i<nTotal; i++) {
        Measurement m = AnalyzeRun(files[i], angles[i]);
        if (m.valid) {
            gGamV->SetPoint(cV, m.angle, m.e_gam); gGamV->SetPointError(cV, 0, m.e_gam_err);
            gEleV->SetPoint(cV, m.angle, m.e_ele); gEleV->SetPointError(cV, 0, m.e_ele_err);
            gSumV->SetPoint(cV, m.angle, m.e_sum); gSumV->SetPointError(cV, 0, m.e_sum_err);
            cV++;
        } else {
            gGamG->SetPoint(cG, m.angle, m.e_gam); gGamG->SetPointError(cG, 0, m.e_gam_err);
            cG++;
        }
    }


    TCanvas *c1 = new TCanvas("c_final", "Physics Results", 1000, 1000);


    TPad *pad1 = new TPad("pad1", "Energies", 0.0, 0.35, 1.0, 0.93);
    TPad *pad2 = new TPad("pad2", "Sum",      0.0, 0.00, 1.0, 0.35);

    pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.02);

    pad2->SetLeftMargin(0.12); pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);

    pad1->Draw(); pad2->Draw();


    pad1->cd(); pad1->SetGrid();
    TH1F *hr1 = pad1->DrawFrame(0, 0, 140, 600);
    hr1->GetYaxis()->SetTitle("Energy [keV]");
    hr1->GetYaxis()->SetTitleSize(0.06); hr1->GetYaxis()->SetLabelSize(0.05);
    hr1->GetXaxis()->SetLabelSize(0);

    // Theory
    TF1 *fGam = new TF1("fGam", TheoryPhoton, 0, 140, 0);
    fGam->SetLineColor(kBlue+2); fGam->SetLineWidth(2); fGam->Draw("SAME"); fGam->SetLineStyle(2);

    TF1 *fEle = new TF1("fEle", TheoryElectron, 0, 140, 0);
    fEle->SetLineColor(kRed+2); fEle->SetLineWidth(2); fEle->Draw("SAME"); fEle->SetLineStyle(2);


    gGamV->SetMarkerStyle(kFullCircle); gGamV->SetMarkerSize(1); gGamV->SetMarkerColor(kBlue); gGamV->SetLineColor(kBlue);
    gGamV->Draw("P E1 SAME");

    gEleV->SetMarkerStyle(kFullCircle); gEleV->SetMarkerSize(1); gEleV->SetMarkerColor(kRed); gEleV->SetLineColor(kRed);
    gEleV->Draw("P E1 SAME");


    gGamG->SetMarkerStyle(kOpenCircle); gGamG->SetMarkerSize(1); gGamG->SetMarkerColor(kBlue); gGamG->SetLineColor(kBlue);
    gGamG->Draw("P E1 SAME");

    TLegend *leg = new TLegend(0.55, 0.65, 0.93, 0.88);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(1);
    leg->AddEntry(fGam, "Photon Theory", "l");
    leg->AddEntry(fEle, "Electron Theory", "l");
    leg->AddEntry(gGamV, "Photon", "p");
    leg->AddEntry(gEleV, "Electron", "p");
    leg->AddEntry(gGamG, "Photon (Artifacts)", "p");
    leg->Draw();


    pad2->cd(); pad2->SetGrid();
    TH1F *hr2 = pad2->DrawFrame(0, 450, 140, 550);
    hr2->GetYaxis()->SetTitle("Sum [keV]");
    hr2->GetXaxis()->SetTitle("Scattering Angle [deg]");

    hr2->GetYaxis()->SetTitleSize(0.08); hr2->GetYaxis()->SetLabelSize(0.08);
    hr2->GetXaxis()->SetTitleSize(0.08); hr2->GetXaxis()->SetLabelSize(0.08);
    hr2->GetYaxis()->SetTitleOffset(0.6);
    hr2->GetYaxis()->SetNdivisions(505);

    TLine *lRef = new TLine(0, 511, 140, 511);
    lRef->SetLineColor(kGreen+2); lRef->SetLineStyle(7); lRef->SetLineWidth(2);
    lRef->Draw();

    gSumV->SetMarkerStyle(kFullCircle); gSumV->SetMarkerSize(1); gSumV->SetMarkerColor(kBlack); gSumV->SetLineColor(kBlack);
    gSumV->Draw("P E1 SAME");


    c1->cd();
    TLatex *title = new TLatex();
    title->SetNDC();
    title->SetTextFont(62);
    title->SetTextSize(0.03);
    title->SetTextAlign(22);

    title->DrawLatex(0.5, 0.965, "");

    c1->SaveAs("ComptonFinalTheory_Fixed.png");
    c1->SaveAs("ComptonFinalTheory_Fixed.pdf");
}