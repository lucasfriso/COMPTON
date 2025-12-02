

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>

Double_t Calibrate(Double_t val, Double_t p0, Double_t p1) {
    return p0 + (p1 * val);
}

void AnalyzeFinal(const char* filename) {




    Double_t calib_Tag_p0 = -13.5147;    Double_t calib_Tag_p1 = 2.351;

    Double_t calib_Det_p0 = -11.0774;    Double_t calib_Det_p1 = 2.412935;

    Double_t calib_Sct_p0 = -14.25;   Double_t calib_Sct_p1 =2.85775;

    // =========================
    // 2. CUTS
    // =========================
    Double_t tagger_min = 470.0;
    Double_t tagger_max = 530.0;
    Double_t scat_threshold = 20;

    // =========================
    // 3. FILE READING
    // =========================
    TFile *f = new TFile(filename);
    if (!f || f->IsZombie()) { std::cout << "Error: File not found." << std::endl; return; }

    TTree *t = (TTree*)f->Get("Data_B");
    if (!t) { std::cout << "Error: Tree 'Data_B' not found." << std::endl; return; }

    TLeaf *l_tag = t->GetLeaf("Channel_0.Energy");
    TLeaf *l_sct = t->GetLeaf("Channel_1.Energy");
    TLeaf *l_det = t->GetLeaf("Channel_2.Energy");

    if (!l_tag || !l_sct || !l_det) {
        l_tag = t->GetLeaf("ch0.Energy");
        l_sct = t->GetLeaf("ch1.Energy");
        l_det = t->GetLeaf("ch2.Energy");
        if (!l_tag) return;
    }

    // =========================
    // 4. HISTOGRAMS
    // =========================

    TH1F *h_tag = new TH1F("h_tag", "Tagger Spectrum (Check Gate);Energy [keV]", 100, 0, 1000);
    TH1F *h_det = new TH1F("h_det", "Scattered Photon;Energy [keV]", 100, 0, 1000);
    TH1F *h_sct = new TH1F("h_sct", "Recoil Electron;Energy [keV]", 100, 0, 1000);
    TH1F *h_sum = new TH1F("h_sum", "Sum Energy;Sum [keV]", 100, 0, 1200);
    TH2F *h_2D  = new TH2F("h_2D",  "Det vs Scat;Det [keV];Scat [keV]", 100, 0, 800, 100, 0, 800);

    Long64_t nEntries = t->GetEntries();
    Long64_t accepted = 0;

    std::cout << "Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        t->GetEntry(i);

        Double_t E_tag = Calibrate(l_tag->GetValue(), calib_Tag_p0, calib_Tag_p1);
        Double_t E_sct = Calibrate(l_sct->GetValue(), calib_Sct_p0, calib_Sct_p1);
        Double_t E_det = Calibrate(l_det->GetValue(), calib_Det_p0, calib_Det_p1);

        h_tag->Fill(E_tag);

        if (E_tag > tagger_min && E_tag < tagger_max) {
            if (E_sct > scat_threshold && E_det > 10) {
                accepted++;
                h_det->Fill(E_det);
                h_sct->Fill(E_sct);
                h_sum->Fill(E_det + E_sct);
                h_2D->Fill(E_det, E_sct);
            }
        }
    }

    // =========================
    // 5. PLOTTING & FITS
    // =========================
    TCanvas *c_check = new TCanvas("c_check", "Tagger Check", 800, 600);
    h_tag->Draw();


    TLine *l_min = new TLine(tagger_min, 0, tagger_min, h_tag->GetMaximum());
    l_min->SetLineColor(kRed); l_min->SetLineWidth(2);
    l_min->Draw();

    TLine *l_max = new TLine(tagger_max, 0, tagger_max, h_tag->GetMaximum());
    l_max->SetLineColor(kRed); l_max->SetLineWidth(2);
    l_max->Draw();

    TCanvas *c1 = new TCanvas("c1", "Compton Analysis", 1200, 800);
    c1->Divide(2, 2);
    gStyle->SetOptFit(1111);

    // --- Pad 1: Detector Fit ---
    c1->cd(1);
    if (h_det->GetEntries() > 0) {

        h_det->GetXaxis()->SetRangeUser(100, 400);
        Double_t peak = h_det->GetBinCenter(h_det->GetMaximumBin());
        h_det->GetXaxis()->UnZoom();

        TF1 *f1 = new TF1("f1", "gaus", peak - 80, peak + 80);
        f1->SetLineColor(kRed);
        h_det->Fit(f1, "R");
        h_det->Draw();
        std::cout << "Photon Energy: " << f1->GetParameter(1) << " keV" << std::endl;
    }

    // --- Pad 2: Electron Fit
    c1->cd(2);
    if (h_sct->GetEntries() > 0) {

        h_sct->GetXaxis()->SetRangeUser(100, 400);
        Double_t peak_sct = h_sct->GetBinCenter(h_sct->GetMaximumBin());
        h_sct->GetXaxis()->UnZoom();

        TF1 *f2 = new TF1("f2", "gaus", peak_sct - 80, peak_sct + 80);
        f2->SetLineColor(kRed);
        h_sct->Fit(f2, "R");
        h_sct->Draw();
        std::cout << "Electron Energy: " << f2->GetParameter(1) << " keV" << std::endl;
    }

    // --- Pad 3: Sum Fit ---
    c1->cd(3);
    if (h_sum->GetEntries() > 0) {
        h_sum->GetXaxis()->SetRangeUser(400, 500); // Look near 511
        Double_t peak_sum = h_sum->GetBinCenter(h_sum->GetMaximumBin());
        h_sum->GetXaxis()->UnZoom();

        TF1 *f3 = new TF1("f3", "gaus", peak_sum - 50, peak_sum + 70);
        f3->SetLineColor(kRed);
        h_sum->Fit(f3, "R");
        h_sum->Draw();
        std::cout << "Sum Energy: " << f3->GetParameter(1) << " keV" << std::endl;
    }

    // --- Pad 4: 2D ---
    c1->cd(4);
    h_2D->Draw("COLZ");

    std::cout << "Total Accepted: " << accepted << std::endl;
}