// Plotter for the post_ALICE output (v2_results.root).
//   build:  g++ -O2 -o post/plot_ALICE post/plot_ALICE.cc `root-config --cflags --libs`
//   run:    ./post/plot_ALICE [data/data_ALICE/v2_results.root] [plots/plot_ALICE]
// Writes: v2_vs_pt.pdf, v2sq_vs_pt.pdf, V2D_pairings.pdf (+ .png of each)
#include <cstdio>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TROOT.h"

static const int NSP = 5;
static const char* SP_NAME[NSP]  = {"pi", "K", "p", "K0S", "Lam"};
static const char* SP_LABEL[NSP] = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}", "K^{0}_{S}", "#Lambda+#bar{#Lambda}"};
// mesons: blue-ish, baryons: red-ish
static const int SP_COLOR[NSP]  = {kBlue + 1, kAzure + 7, kRed + 1, kCyan + 2, kOrange + 7};
static const int SP_MARKER[NSP] = {20, 21, 24, 25, 22};

static TFile* gF = nullptr;

static TGraphErrors* styled(const std::string& name, int s) {
    TGraphErrors* g = (TGraphErrors*)gF->Get(name.c_str());
    if (!g || g->GetN() == 0) return nullptr;
    g->SetMarkerStyle(SP_MARKER[s]);
    g->SetMarkerColor(SP_COLOR[s]);
    g->SetLineColor(SP_COLOR[s]);
    g->SetMarkerSize(1.1);
    return g;
}

// one canvas, all species of one quantity vs pT
static void drawAll(const char* prefix, const char* ytitle, const std::string& out) {
    TCanvas c("c", "", 800, 600);
    c.SetLeftMargin(0.15); c.SetBottomMargin(0.13);
    TMultiGraph mg;
    TLegend leg(0.62, 0.62, 0.89, 0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    int n = 0;
    for (int s = 0; s < NSP; ++s) {
        TGraphErrors* g = styled(Form("%s_%s", prefix, SP_NAME[s]), s);
        if (!g) { fprintf(stderr, "  [no points: %s_%s]\n", prefix, SP_NAME[s]); continue; }
        mg.Add(g, "PZ");
        leg.AddEntry(g, SP_LABEL[s], "p");
        ++n;
    }
    if (!n) { fprintf(stderr, "Nothing to draw for %s\n", prefix); return; }
    mg.Draw("A");
    mg.GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mg.GetYaxis()->SetTitle(ytitle);
    mg.GetXaxis()->SetLimits(0.0, 10.0);
    mg.GetXaxis()->SetTitleSize(0.05); mg.GetYaxis()->SetTitleSize(0.05);
    mg.GetYaxis()->SetTitleOffset(1.4);
    c.Modified();
    TLine zero(0.0, 0.0, 10.0, 0.0);
    zero.SetLineStyle(2); zero.SetLineColor(kGray + 2);
    zero.Draw();
    leg.Draw();
    c.SaveAs((out + ".pdf").c_str());
    c.SaveAs((out + ".png").c_str());
}

// V2Delta of the two TPC pairings, side by side
static void drawPairings(const std::string& out) {
    TCanvas c("c2", "", 1400, 600);
    c.Divide(2, 1);
    const char* pair[2]  = {"V2D_TPC_FMD12", "V2D_TPC_FMD3"};
    const char* title[2] = {"V_{2#Delta}  TPC-FMD1,2", "V_{2#Delta}  TPC-FMD3"};
    TMultiGraph* mg[2];
    TLegend* leg[2];
    for (int p = 0; p < 2; ++p) {
        c.cd(p + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.13);
        mg[p] = new TMultiGraph();
        leg[p] = new TLegend(0.60, 0.62, 0.89, 0.89);
        leg[p]->SetBorderSize(0); leg[p]->SetFillStyle(0);
        int n = 0;
        for (int s = 0; s < NSP; ++s) {
            TGraphErrors* g = styled(Form("%s_%s", pair[p], SP_NAME[s]), s);
            if (!g) continue;
            mg[p]->Add(g, "PZ");
            leg[p]->AddEntry(g, SP_LABEL[s], "p");
            ++n;
        }
        if (!n) continue;
        mg[p]->Draw("A");
        mg[p]->SetTitle(title[p]);
        mg[p]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        mg[p]->GetYaxis()->SetTitle("V_{2#Delta}");
        mg[p]->GetXaxis()->SetLimits(0.0, 10.0);
        mg[p]->GetXaxis()->SetTitleSize(0.05); mg[p]->GetYaxis()->SetTitleSize(0.05);
        mg[p]->GetYaxis()->SetTitleOffset(1.7);
        TLine* zero = new TLine(0.0, 0.0, 10.0, 0.0);
        zero->SetLineStyle(2); zero->SetLineColor(kGray + 2);
        zero->Draw();
        leg[p]->Draw();
    }
    c.SaveAs((out + ".pdf").c_str());
    c.SaveAs((out + ".png").c_str());
}

int main(int argc, char* argv[]) {
    const std::string inFile = (argc > 1) ? argv[1] : "data/data_ALICE/v2_results.root";
    const std::string outDir = (argc > 2) ? argv[2] : "plots/plot_ALICE";

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    gF = TFile::Open(inFile.c_str());
    if (!gF || gF->IsZombie()) { fprintf(stderr, "Cannot open %s\n", inFile.c_str()); return 1; }

    drawAll("v2",   "v_{2}  (signed sqrt)", outDir + "/v2_vs_pt");
    drawAll("v2sq", "v_{2}^{2}  (signed)",  outDir + "/v2sq_vs_pt");
    drawPairings(outDir + "/V2D_pairings");

    gF->Close();
    printf("Plots written to %s/\n", outDir.c_str());
    return 0;
}
