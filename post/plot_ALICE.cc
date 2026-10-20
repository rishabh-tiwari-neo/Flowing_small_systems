// plot_ALICE.cc
// Makes v2sq_vs_pt, V2D_TPC_FMD12, V2D_TPC_FMD3 and V2D_FMD12_FMD3.

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TParameter.h"
#include "TLatex.h"

// ---- input / output locations and the classes to plot ----
static const char* DATA_DIR      = "data/data_ALICE";
static const char* PLOT_DIR_BASE = "plots/plot_ALICE";
static const std::vector<std::string> CLASSES = {"HM", "HM35"};

// ---- species: names (must match post_ALICE), legend labels, colours, markers ----
static const int NSP = 5;
static const char* SP_NAME[NSP]  = {"pi", "K", "p", "K0S", "Lam"};
static const char* SP_LABEL[NSP] = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}",
                                    "K^{0}_{S}", "#Lambda+#bar{#Lambda}"};
static const int   SP_COLOR[NSP] = {kBlue, kAzure, kRed, kCyan, kOrange};
static const int   SP_MARKER[NSP] = {20, 34, 21, 22, 33};

// reads the shoving strength g stored in merged.root by gen_ALICE; returns -1 if absent (older data)
static double readG() {
    TFile* f = TFile::Open((std::string(DATA_DIR) + "/merged.root").c_str(), "READ");
    if (!f || f->IsZombie()) return -1;
    TParameter<double>* p = dynamic_cast<TParameter<double>*>(f->Get("g"));
    const double g = p ? p->GetVal() : -1;
    f->Close();
    return g;
}

// gives every axis frame a unique name
static int frameCounter = 0;

// pad margins shared by all plots
static void setMargins() {
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.06);
    gPad->SetRightMargin(0.04);
}

// Draws all species of one kind (prefix = v2sq, V2D_TPC_FMD12, V2D_TPC_FMD3)

static void drawPad(TFile* f, const char* prefix, const char* ytitle, double ylo, double yhi, double gShove)
{
    TGraphErrors* g[NSP];
    double xlo = 0, xhi = 0;
    bool any = false;

    // ---- fetch the graphs and find the x range ----
    for (int s = 0; s < NSP; ++s) {
        g[s] = dynamic_cast<TGraphErrors*>(f->Get(Form("%s_%s", prefix, SP_NAME[s])));
        if (!g[s] || g[s]->GetN() == 0) { g[s] = nullptr; continue; }   // missing or empty: skip this species
        any = true;
        for (int i = 0; i < g[s]->GetN(); ++i)
            xhi = std::max(xhi, g[s]->GetX()[i] + g[s]->GetErrorX(i));
    }
    if (!any) { std::cerr << "No graphs found for " << prefix << std::endl; return; }

    setMargins();

    // ---- empty frame that carries the axes, titles and y range (set by hand: ylo, yhi) ----
    TH1F* fr = new TH1F(Form("frame%d", frameCounter++), "", 100, xlo, xhi);
    fr->SetStats(0);
    fr->SetMinimum(ylo);
    fr->SetMaximum(yhi);
    fr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fr->GetYaxis()->SetTitle(ytitle);
    fr->GetYaxis()->CenterTitle(true);      // y title in the middle of the axis
    fr->GetYaxis()->SetTitleOffset(1.9);
    fr->Draw("AXIS");

    // ---- dashed zero line ----
    TLine* zero = new TLine(xlo, 0, xhi, 0);
    zero->SetLineStyle(2);
    zero->Draw();

    // ---- draw each species with its own style and add it to the legend ----
    TLegend* leg = new TLegend(0.84, gShove >= 0 ? 0.64 : 0.70, 0.97, 0.92);
    if (gShove >= 0) leg->SetHeader(Form("g = %g", gShove));
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);       
    leg->SetEntrySeparation(0.09);
    leg->SetMargin(0.30);
    for (int s = 0; s < NSP; ++s) {
        if (!g[s]) continue;
        g[s]->SetMarkerStyle(SP_MARKER[s]);
        g[s]->SetMarkerColor(SP_COLOR[s]);
        g[s]->SetLineColor(SP_COLOR[s]);
        g[s]->SetMarkerSize(1.6);
        g[s]->Draw("PZ SAME");
        leg->AddEntry(g[s], SP_LABEL[s], "p");
    }
    leg->Draw();
}

// Draws the single FMD1,2-FMD3 point into the current pad.
static void drawRef(TFile* f, double g) {
    TGraphErrors* gr = dynamic_cast<TGraphErrors*>(f->Get("V2D_FMD12_FMD3"));
    if (!gr || gr->GetN() == 0) { std::cerr << "No graph found for V2D_FMD12_FMD3" << std::endl; return; }

    // ---- y range: include zero, the point and its error bar, plus 10% margin ----
    const double y = gr->GetY()[0], ey = gr->GetErrorY(0);
    const double lo = std::min(0.0, y - ey), hi = std::max(0.0, y + ey);
    const double margin = 0.1 * (hi - lo > 0 ? hi - lo : 1e-3);

    setMargins();

    TH1F* fr = new TH1F(Form("frame%d", frameCounter++), "", 1, -1, 1);
    fr->SetStats(0);
    fr->SetMinimum(lo - margin);
    fr->SetMaximum(hi + margin);
    fr->GetXaxis()->SetTitle("FMD1,2 - FMD3");
    fr->GetXaxis()->SetNdivisions(0);
    fr->GetYaxis()->SetTitle("V_{2#Delta} (FMD1,2-FMD3)");
    fr->GetYaxis()->CenterTitle(true);
    fr->GetYaxis()->SetTitleOffset(1.9);
    fr->Draw("AXIS");

    // ---- dashed zero line ----
    TLine* zero = new TLine(-1, 0, 1, 0);
    zero->SetLineStyle(2);
    zero->Draw();

    // ---- the point itself ----
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->SetLineColor(kBlack);
    gr->SetMarkerSize(1.3);
    gr->Draw("PZ SAME");

    if (g >= 0) {
        TLatex* t = new TLatex(0.20, 0.88, Form("g = %g", g));
        t->SetNDC();
        t->SetTextSize(0.04);
        t->Draw();
    }
}

// saves one canvas as PDF and PNG, then deletes it (canvas names stay unique per class)
static void save(TCanvas* c, const std::string& outdir, const char* name) {
    c->SaveAs((outdir + "/" + name + ".pdf").c_str());
    c->SaveAs((outdir + "/" + name + ".png").c_str());
    delete c;
}

// Opens the results file of one class and writes its five plots. Returns false on failure.
static bool plotClass(const std::string& cls) {
    const std::string inFile = std::string(DATA_DIR) + "/v2_results_" + cls + ".root";
    const std::string outdir = std::string(PLOT_DIR_BASE) + "_" + cls;

    TFile* f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << inFile << std::endl;
        return false;
    }
    const double g = readG();
    gSystem->mkdir(outdir.c_str(), true);       // creates the folder (and parents) if missing

    // v2^2 vs pT
    TCanvas* c2 = new TCanvas(Form("c_v2sq_%s", cls.c_str()), "v2sq", 800, 600);
    c2->cd();
    drawPad(f, "v2sq", "v_{2}^{2}", -1.0, 1.0, g);
    save(c2, outdir, "v2sq_vs_pt");

    // V2D TPC-FMD1,2 vs pT
    TCanvas* c3 = new TCanvas(Form("c_tpc12_%s", cls.c_str()), "V2D TPC-FMD12", 800, 600);
    c3->cd();
    drawPad(f, "V2D_TPC_FMD12", "V_{2#Delta} (TPC-FMD1,2)", -0.05, 0.05, g);
    save(c3, outdir, "V2D_TPC_FMD12");

    // V2D TPC-FMD3 vs pT
    TCanvas* c4 = new TCanvas(Form("c_tpc3_%s", cls.c_str()), "V2D TPC-FMD3", 800, 600);
    c4->cd();
    drawPad(f, "V2D_TPC_FMD3", "V_{2#Delta} (TPC-FMD3)", -0.05, 0.05, g);
    save(c4, outdir, "V2D_TPC_FMD3");

    // V2D FMD1,2-FMD3 (single number)
    TCanvas* c5 = new TCanvas(Form("c_fmdfmd_%s", cls.c_str()), "V2D FMD12-FMD3", 800, 600);
    c5->cd();
    drawRef(f, g);
    save(c5, outdir, "V2D_FMD12_FMD3");

    f->Close();
    std::cout << "Class " << cls << ": plots written to " << outdir << std::endl;
    return true;
}

int main(int argc, char* argv[]) {
    // ---- arguments and ROOT setup ----
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [class]   (no argument = all classes)" << std::endl;
        return 1;
    }
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    // ---- which classes to plot: all of them, or the one given on the command line ----
    std::vector<std::string> todo = CLASSES;
    if (argc == 2) todo = {argv[1]};

    // ---- load and plot the classes one by one ----
    int failed = 0;
    for (const std::string& cls : todo)
        if (!plotClass(cls)) ++failed;

    return failed ? 1 : 0;
}