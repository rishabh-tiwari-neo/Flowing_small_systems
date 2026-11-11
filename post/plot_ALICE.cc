// plot_ALICE.cc
// Makes v2sq_vs_pt, v2_vs_pt, V2D_TPC_FMD12, V2D_TPC_FMD3, V2D_FMD12_FMD3,
// and the NCQ-scaling plots NCQ_pT (v2/nq vs pT/nq) and NCQ_KET (v2/nq vs KET/nq).
// Everything is read from v2_results_<class>.root (the "v2_<sp>" graphs post_ALICE
// writes); no separate tool or CSV parsing needed for the NCQ plots.
//
// v2 = sqrt(v2sq) is always positive, so it and the NCQ plots below lose the sign of
// the underlying V2Delta -- check v2sq_vs_pt or the per-pairing V2D plots for that.

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
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

// ---- for NCQ scaling only: PDG mass (GeV/c^2) and number of constituent quarks ----
static const double SP_MASS[NSP] = {0.13957, 0.49368, 0.93827, 0.49761, 1.11568};
static const int    SP_NQ[NSP]   = {2, 2, 3, 2, 3};

// significance cut for the NCQ plots: keep a point only if v2/err(v2) exceeds this.
// (equivalent to a cut of MIN_SIG_V2/2 on v2sq/err(v2sq), since err(v2) = err(v2sq)/(2 v2))
static const double MIN_SIG_V2 = 4.0;

static inline double ket(double pt, double m) { return std::sqrt(pt * pt + m * m) - m; }

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

// Draws the NCQ-scaling plot: v2/nq vs pT/nq (kind 0) or v2/nq vs KET/nq (kind 1).
// Reads the "v2_<sp>" graphs post_ALICE writes (x = pT bin centre, y = v2, with their
// errors) and rescales each point by the species' nq; for kind 1 the x axis becomes
// KET = sqrt(pT^2+m^2) - m, evaluated at the bin centre, with the x error propagated
// from the bin edges (x -+ errorX) through the same KET transform.
static void drawNCQ(TFile* f, int kind, const std::string& cls, double gShove) {
    TGraphErrors* g[NSP];
    double xhi = 0, yhi = 0;
    bool any = false;

    for (int s = 0; s < NSP; ++s) {
        TGraphErrors* src = dynamic_cast<TGraphErrors*>(f->Get(Form("v2_%s", SP_NAME[s])));
        g[s] = new TGraphErrors();
        if (!src || src->GetN() == 0) continue;
        const double m = SP_MASS[s], nq = SP_NQ[s];
        for (int i = 0; i < src->GetN(); ++i) {
            const double pt = src->GetX()[i], ept = src->GetErrorX(i);
            const double v2 = src->GetY()[i], ev2 = src->GetErrorY(i);
            if (ev2 <= 0 || v2 / ev2 < MIN_SIG_V2) continue;
            const double lo = pt - ept, hi = pt + ept;
            double x, ex;
            if (kind == 0) { x = pt / nq;        ex = ept / nq; }
            else           { x = ket(pt, m) / nq; ex = 0.5 * (ket(hi, m) - ket(lo, m)) / nq; }
            const int n = g[s]->GetN();
            g[s]->SetPoint(n, x, v2 / nq);
            g[s]->SetPointError(n, ex, ev2 / nq);
            xhi = std::max(xhi, x + ex);
            yhi = std::max(yhi, (v2 + ev2) / nq);
        }
        if (g[s]->GetN() > 0) any = true;
    }
    if (!any) { std::cerr << "No points pass the NCQ cuts for " << cls << " (kind " << kind << ")" << std::endl; return; }

    setMargins();

    TH1F* fr = new TH1F(Form("frameN%d", frameCounter++), "", 100, 0, xhi * 1.05);
    fr->SetStats(0);
    fr->SetMinimum(0);
    fr->SetMaximum(yhi * 1.15);
    fr->GetXaxis()->SetTitle(kind == 0 ? "p_{T}/n_{q} (GeV/c)" : "KE_{T}/n_{q} (GeV)");
    fr->GetYaxis()->SetTitle("v_{2}/n_{q}");
    fr->GetYaxis()->CenterTitle(true);
    fr->GetYaxis()->SetTitleOffset(1.9);
    fr->Draw("AXIS");

    TLegend* leg = new TLegend(0.18, gShove >= 0 ? 0.62 : 0.68, 0.45, 0.92);
    if (gShove >= 0) leg->SetHeader(Form("%s, g = %g", cls.c_str(), gShove));
    else             leg->SetHeader(cls.c_str());
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04); leg->SetMargin(0.30);
    for (int s = 0; s < NSP; ++s) {
        if (g[s]->GetN() == 0) continue;
        g[s]->SetMarkerStyle(SP_MARKER[s]); g[s]->SetMarkerColor(SP_COLOR[s]);
        g[s]->SetLineColor(SP_COLOR[s]);    g[s]->SetMarkerSize(1.6);
        g[s]->Draw("PZ SAME");
        leg->AddEntry(g[s], SP_LABEL[s], "p");
    }
    leg->Draw();
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

    // v2 = sqrt(v2sq) vs pT (always positive -- see the note at the top of this file)
    TCanvas* c2b = new TCanvas(Form("c_v2_%s", cls.c_str()), "v2", 800, 600);
    c2b->cd();
    drawPad(f, "v2", "v_{2}", 0.0, 1, g);
    save(c2b, outdir, "v2_vs_pt");

    // NCQ scaling: v2/nq vs pT/nq and v2/nq vs KET/nq
    TCanvas* c2c = new TCanvas(Form("c_ncqpt_%s", cls.c_str()), "NCQ pT", 800, 600);
    c2c->cd();
    drawNCQ(f, 0, cls, g);
    save(c2c, outdir, "NCQ_pT");

    TCanvas* c2d = new TCanvas(Form("c_ncqket_%s", cls.c_str()), "NCQ KET", 800, 600);
    c2d->cd();
    drawNCQ(f, 1, cls, g);
    save(c2d, outdir, "NCQ_KET");

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