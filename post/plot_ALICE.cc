// plot_ALICE.cc
// Makes v2sq_vs_pt, v2_vs_pt, and V2D_pt from pairings.

#include <iostream>
#include <string>
#include <algorithm>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

static const int NSP = 5;
static const char* SP_NAME[NSP]  = {"pi", "K", "p", "K0S", "Lam"};
static const char* SP_LABEL[NSP] = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}",
                                    "K^{0}_{S}", "#Lambda+#bar{#Lambda}"};
static const int   SP_COLOR[NSP] = {kBlue, kAzure, kRed, kCyan, kOrange};
static const int   SP_MARKER[NSP] = {20, 34, 21, 22, 33};

static int frameCounter = 0;

// Draws all species of one kind (prefix = v2, v2sq, V2D_TPC_FMD12, V2D_TPC_FMD3)
// into the current pad. clamp > 0 limits the y axis to +-clamp.
static void drawPad(TFile* f, const char* prefix, const char* ytitle,
                    double clamp, const char* refText) {
    TGraphErrors* g[NSP];
    double lo = 0, hi = 0, xlo = 0, xhi = 0;
    bool any = false;

    for (int s = 0; s < NSP; ++s) {
        g[s] = dynamic_cast<TGraphErrors*>(f->Get(Form("%s_%s", prefix, SP_NAME[s])));
        if (!g[s] || g[s]->GetN() == 0) { g[s] = nullptr; continue; }
        any = true;
        for (int i = 0; i < g[s]->GetN(); ++i) {
            const double y = g[s]->GetY()[i], ey = g[s]->GetErrorY(i);
            const double x = g[s]->GetX()[i], ex = g[s]->GetErrorX(i);
            lo  = std::min(lo,  y - ey);
            hi  = std::max(hi,  y + ey);
            xhi = std::max(xhi, x + ex);
        }
    }
    if (!any) { std::cerr << "No graphs found for " << prefix << std::endl; return; }

    if (clamp > 0) { lo = std::max(lo, -clamp); hi = std::min(hi, clamp); }
    if (hi <= lo) { hi = lo + 1e-3; }
    const double margin = 0.1 * (hi - lo);

    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.06);
    gPad->SetRightMargin(0.04);

    TH1F* fr = new TH1F(Form("frame%d", frameCounter++), "", 100, xlo, xhi);
    fr->SetStats(0);
    fr->SetMinimum(lo - margin);
    fr->SetMaximum(hi + margin);
    fr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fr->GetYaxis()->SetTitle(ytitle);
    fr->GetYaxis()->SetTitleOffset(1.9);
    fr->Draw("AXIS");

    TLine* zero = new TLine(xlo, 0, xhi, 0);
    zero->SetLineStyle(2);
    zero->Draw();

    TLegend* leg = new TLegend(0.60, 0.68, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (int s = 0; s < NSP; ++s) {
        if (!g[s]) continue;
        g[s]->SetMarkerStyle(SP_MARKER[s]);
        g[s]->SetMarkerColor(SP_COLOR[s]);
        g[s]->SetLineColor(SP_COLOR[s]);
        g[s]->SetMarkerSize(1.1);
        g[s]->Draw("PZ SAME");
        leg->AddEntry(g[s], SP_LABEL[s], "p");
    }
    leg->Draw();

    if (refText && refText[0]) {
        TLatex* t = new TLatex();
        t->SetNDC();
        t->SetTextSize(0.035);
        t->DrawLatex(0.19, 0.89, refText);
    }
}

static void save(TCanvas* c, const std::string& outdir, const char* name) {
    c->SaveAs((outdir + "/" + name + ".pdf").c_str());
    c->SaveAs((outdir + "/" + name + ".png").c_str());
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " results.root outdir" << std::endl;
        return 1;
    }
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    TFile* f = TFile::Open(argv[1], "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << argv[1] << std::endl;
        return 1;
    }
    const std::string outdir = argv[2];

    // reference V2D(FMD12-FMD3), shown as text on every plot
    std::string ref;
    TGraphErrors* gr = dynamic_cast<TGraphErrors*>(f->Get("V2D_FMD12_FMD3"));
    if (gr && gr->GetN() > 0)
        ref = Form("V_{2#Delta}(FMD1,2-FMD3) = %.2e #pm %.2e", gr->GetY()[0], gr->GetErrorY(0));

    // v2 vs pT
    TCanvas* c1 = new TCanvas("c_v2", "v2", 800, 600);
    c1->cd();
    drawPad(f, "v2", "v_{2}", 1.0, ref.c_str());
    save(c1, outdir, "v2_vs_pt");

    // v2^2 vs pT
    TCanvas* c2 = new TCanvas("c_v2sq", "v2sq", 800, 600);
    c2->cd();
    drawPad(f, "v2sq", "v_{2}^{2}", 1.0, ref.c_str());
    save(c2, outdir, "v2sq_vs_pt");

    // V2D of the two identified pairings side by side
    TCanvas* c3 = new TCanvas("c_pair", "V2D pairings", 1400, 600);
    c3->Divide(2, 1);
    c3->cd(1);
    drawPad(f, "V2D_TPC_FMD12", "V_{2#Delta} (TPC-FMD1,2)", 0.05, ref.c_str());
    c3->cd(2);
    drawPad(f, "V2D_TPC_FMD3", "V_{2#Delta} (TPC-FMD3)", 0.05, ref.c_str());
    save(c3, outdir, "V2D_pairings");

    f->Close();
    std::cout << "Plots written to " << outdir << std::endl;
    return 0;
}
