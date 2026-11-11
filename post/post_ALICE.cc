// post_ALICE.cc
// Builds per-trigger yields Y(dphi) = S/B projected over deta, then fits them with the template.
// S and B are 2D (dphi x deta) histograms from gen_ALICE. B is normalised with one constant:
// its dphi average at the deta bin where B is largest. S is divided by the number of triggers
// (for FMD1,2-FMD3: the number of FMD1,2 particles).
// Weighted linear least-squares fitting with analytical covariance-based error estimation.
// The LM template is treated as exact, so its error is not in the V2Δ uncertainty.
// Any species can be analysed in coarser pT bins: the 2D S and B histograms of
// several generator pT bins are summed before dividing (see the REBIN table and ptGroups).
// The three V2Δ (TPC-FMD1,2, TPC-FMD3, FMD1,2-FMD3) are combined into v2^2 (3x2PC).
// Plain v2 = sqrt(v2^2) is also written (graph "v2_<sp>", CSV columns v2/err), defined
// only where v2^2 > 0. v2 is always positive and so loses the sign of the underlying
// V2Delta -- check v2^2 or the per-pairing V2Delta before drawing conclusions from v2.
//
// Per-pairing deta window: gen_ALICE.cc fills all three pairings into the same 1.1-6.4
// deta axis, but each pairing can only physically produce pairs in a narrower range set
// by its two eta acceptances (trigger |eta|<0.8, FMD1,2 1.7-5.1, FMD3 -3.1 to -1.7).
// Outside that kinematic range a deta slice is exactly empty and is skipped anyway
// (B<=0). Just inside the edge, though, the number of (eta_trig, eta_ref) combinations
// that can produce a deta that close to the kinematic maximum tapers toward zero, so a
// handful of slices there can have only a few raw pairs even with huge total statistics.
// Those thin slices get heavy weight in the variance sum (it scales like 1/B^4) and can
// dominate a point's error, or even make the template fit unable to tell the template
// scale F apart from the flow amplitude a2 (seen as |rho| close to 1 in diag_ALICE.cc).
// WIN_* below trims DETA_EDGE_MARGIN off each kinematic edge to stay out of that taper.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <utility>
#include <memory>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TString.h"

static const int NSP = 5;
static const char* SP_NAME[NSP] = {"pi", "K", "p", "K0S", "Lam"};
static const char* DATA_DIR = "data/data_ALICE";
static const std::vector<std::string> CLASSES = {"HM", "HM35"};

// kinematic deta = eta_ref - eta_trig (or eta_a - eta_b) ranges, from the acceptance
// cuts in gen_ALICE.cc: trigger |eta|<0.8, FMD1,2 1.7<eta<5.1, FMD3 -3.1<eta<-1.7.
//   TPC-FMD1,2:  deta = eta_FMD12 - eta_trig  -> (1.7-0.8, 5.1+0.8)   = (0.9, 5.9)
//   TPC-FMD3:    deta = eta_trig - eta_FMD3   -> (-0.8-(-1.7), 0.8-(-3.1)) = (0.9, 3.9)
//   FMD1,2-FMD3: deta = eta_FMD12 - eta_FMD3  -> (1.7-(-1.7), 5.1-(-3.1))  = (3.4, 8.2)
// margin trimmed off each edge to stay clear of the low-count taper near the kinematic
// boundary (see the comment at the top of this file). Tune DETA_EDGE_MARGIN if needed.
static const double DETA_EDGE_MARGIN = 0.4;
struct DetaWindow { double lo, hi; };
static const DetaWindow WIN_TPC_FMD12  = {0.9 + DETA_EDGE_MARGIN, 5.9 - DETA_EDGE_MARGIN};
static const DetaWindow WIN_TPC_FMD3   = {0.9 + DETA_EDGE_MARGIN, 3.9 - DETA_EDGE_MARGIN};
static const DetaWindow WIN_FMD12_FMD3 = {3.4 + DETA_EDGE_MARGIN, 8.2 - DETA_EDGE_MARGIN};

// initialising the fit results
struct FitResult {
    bool   ok      = false;
    double V2      = 0;
    double err     = 0;
    double chi2ndf = 0;
};

// 4x4 inversion (Gauss-Jordan with row swaps)
static bool invert4(const double a[4][4], double out[4][4]) {
    double m[4][8];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            m[i][j]     = a[i][j];
            m[i][j + 4] = (i == j) ? 1.0 : 0.0;
        }
    for (int c = 0; c < 4; ++c) {
        int piv = c;
        for (int r = c + 1; r < 4; ++r)
            if (std::fabs(m[r][c]) > std::fabs(m[piv][c])) piv = r;
        if (std::fabs(m[piv][c]) < 1e-300) return false;
        if (piv != c) for (int k = 0; k < 8; ++k) std::swap(m[c][k], m[piv][k]);
        const double d = m[c][c];
        for (int k = 0; k < 8; ++k) m[c][k] /= d;
        for (int r = 0; r < 4; ++r) {
            if (r == c) continue;
            const double f = m[r][c];
            if (f == 0.0) continue;
            for (int k = 0; k < 8; ++k) m[r][k] -= f * m[c][k];
        }
    }
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) out[i][j] = m[i][j + 4];
    return true;
}

static FitResult templateFit(const TH1D* hs, const TH1D* ht) {
    FitResult res;
    if (!hs || !ht) return res;
    if (hs->Integral() <= 0 || ht->Integral() <= 0) return res;

    // inputs are already per-trigger yields, so no further normalisation
    std::vector<double> y, e, t, ph;
    for (int i = 1; i <= hs->GetNbinsX(); ++i) {
        const double err = hs->GetBinError(i);
        if (err <= 0) continue;
        y.push_back(hs->GetBinContent(i));
        e.push_back(err);
        t.push_back(ht->GetBinContent(i));
        ph.push_back(hs->GetBinCenter(i));
    }
    const int npts = (int)y.size();
    if (npts < 5) return res;

    double A[4][4] = {{0}}, B[4] = {0};
    for (int i = 0; i < npts; ++i) {
        const double f[4] = {t[i], 1.0, std::cos(2.0 * ph[i]), std::cos(3.0 * ph[i])};
        const double w = 1.0 / (e[i] * e[i]);
        for (int a = 0; a < 4; ++a) {
            B[a] += w * f[a] * y[i];
            for (int b = 0; b < 4; ++b) A[a][b] += w * f[a] * f[b];
        }
    }
    double C[4][4];
    if (!invert4(A, C)) return res;

    double p[4];
    for (int a = 0; a < 4; ++a) {
        p[a] = 0;
        for (int b = 0; b < 4; ++b) p[a] += C[a][b] * B[b];
    }
    const double G = p[1], a2 = p[2];
    if (G <= 0) return res;

    // V2_delta = a2 / (2G)
    const double dA = 1.0 / (2.0 * G);
    const double dG = -a2 / (2.0 * G * G);
    const double var = dA * dA * C[2][2] + dG * dG * C[1][1] + 2.0 * dA * dG * C[2][1];
    if (var < 0) return res;

    double chi2 = 0;
    for (int i = 0; i < npts; ++i) {
        const double f[4] = {t[i], 1.0, std::cos(2.0 * ph[i]), std::cos(3.0 * ph[i])};
        double m = 0;
        for (int a = 0; a < 4; ++a) m += p[a] * f[a];
        const double d = (y[i] - m) / e[i];
        chi2 += d * d;
    }

    res.ok      = true;
    res.V2      = a2 / (2.0 * G);
    res.err     = std::sqrt(var);
    res.chi2ndf = chi2 / (npts - 4);
    return res;
}

static TH1D* getH(TFile* f, const TString& name) {
    TH1D* h = dynamic_cast<TH1D*>(f->Get(name));
    if (!h) std::cerr << "WARNING: missing histogram " << name << std::endl;
    return h;
}

static TH2D* getH2(TFile* f, const TString& name) {
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name));
    if (!h) std::cerr << "WARNING: missing histogram " << name << std::endl;
    return h;
}

// pT rebinning table, one entry per species in SP_NAME order: pi, K, p, K0S, Lam.
// Each {first, last} is an inclusive range of generator pT bins (0..19) merged into one output bin.
// An empty list keeps the 20 generator bins unchanged.
static const std::vector<std::vector<std::pair<int, int>>> REBIN = {
    /* pi  */ {},
    /* K   */ {},
    /* p   */ {},
    /* K0S */ {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 14}, {15, 16}, {17, 18}, {19, 19}},
    /* Lam */ {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 14}, {15, 16}, {17, 18}, {19, 19}}
};

static std::vector<std::pair<int, int>> ptGroups(int species, int nbGen) {
    if (nbGen == 20 && species >= 0 && species < NSP && !REBIN[species].empty())
        return REBIN[species];
    std::vector<std::pair<int, int>> g;
    for (int b = 0; b < nbGen; ++b) g.push_back({b, b});
    return g;
}

// names of the 2D histograms of generator pT bins b0..b1
static std::vector<TString> binNames(const char* kind, const std::string& cls, const char* pair,
                                     const char* sp, int b0, int b1) {
    std::vector<TString> v;
    for (int b = b0; b <= b1; ++b) v.push_back(Form("%s_%s_%s_%s_pt%d", kind, cls.c_str(), pair, sp, b));
    return v;
}

// Per-trigger yield: (S / (B/c)) summed over deta, divided by ntrig.
// S and B are the sums of the listed histograms (merged pT bins).
// c = dphi average of B at the deta bin where B is largest.
// detaLo/detaHi restrict the deta slices that are summed, to the pairing's trimmed
// kinematic window (see WIN_TPC_FMD12 etc. above); a slice whose bin centre falls
// outside [detaLo, detaHi] is skipped entirely, same as an empty (B<=0) slice.
// c is computed from, and val/var are summed over, only the slices inside the window,
// so a thin slice just beyond the cut does not leak into either the normalisation or
// the fit input. Defaults keep the full axis, for any caller that does not pass a window.
static std::unique_ptr<TH1D> yieldDphi(TFile* f, const std::vector<TString>& namesS,
                                       const std::vector<TString>& namesB, double ntrig,
                                       double detaLo = -1e9, double detaHi = 1e9) {
    if (namesS.empty() || namesS.size() != namesB.size() || ntrig <= 0) return nullptr;
    TH2D* h0 = getH2(f, namesS[0]);
    if (!h0) return nullptr;
    const int nx = h0->GetNbinsX(), ny = h0->GetNbinsY();
    const double xmin = h0->GetXaxis()->GetXmin(), xmax = h0->GetXaxis()->GetXmax();

    // which deta slices (y bins) fall inside the requested window
    std::vector<bool> inWin(ny, true);
    for (int iy = 1; iy <= ny; ++iy) {
        const double center = h0->GetYaxis()->GetBinCenter(iy);
        inWin[iy - 1] = (center >= detaLo && center <= detaHi);
    }

    std::vector<double> S(nx * ny, 0.0), eS2(nx * ny, 0.0), B(nx * ny, 0.0), eB2(nx * ny, 0.0);
    for (size_t k = 0; k < namesS.size(); ++k) {
        TH2D* hS = getH2(f, namesS[k]);
        TH2D* hB = getH2(f, namesB[k]);
        if (!hS || !hB) return nullptr;
        for (int ix = 1; ix <= nx; ++ix)
            for (int iy = 1; iy <= ny; ++iy) {
                if (!inWin[iy - 1]) continue;
                const int i = (ix - 1) * ny + (iy - 1);
                S[i] += hS->GetBinContent(ix, iy);  eS2[i] += std::pow(hS->GetBinError(ix, iy), 2);
                B[i] += hB->GetBinContent(ix, iy);  eB2[i] += std::pow(hB->GetBinError(ix, iy), 2);
            }
    }

    double c = 0;
    for (int iy = 1; iy <= ny; ++iy) {
        if (!inWin[iy - 1]) continue;
        double s = 0;
        for (int ix = 1; ix <= nx; ++ix) s += B[(ix - 1) * ny + (iy - 1)];
        s /= nx;
        if (s > c) c = s;
    }
    if (c <= 0) return nullptr;

    std::unique_ptr<TH1D> h(new TH1D(Form("y_%s", namesS[0].Data()), "", nx, xmin, xmax));
    for (int ix = 1; ix <= nx; ++ix) {
        double val = 0, var = 0;
        for (int iy = 1; iy <= ny; ++iy) {
            if (!inWin[iy - 1]) continue;
            const int i = (ix - 1) * ny + (iy - 1);
            const double b = B[i];
            if (b <= 0) continue;
            const double s = S[i];
            const double k = c / b;
            val += s * k;
            var += eS2[i] * k * k + (s * c) * (s * c) * eB2[i] / (b * b * b * b);
        }
        h->SetBinContent(ix, val / ntrig);
        h->SetBinError(ix, std::sqrt(var) / ntrig);
    }
    return h;
}

static double count1(TFile* f, const TString& name) {
    TH1D* h = getH(f, name);
    return h ? h->GetBinContent(1) : 0.0;
}

static std::string num(double v, bool ok) {
    if (!ok) return "nan";
    char b[32];
    std::snprintf(b, sizeof b, "%.5e", v);
    return b;
}

static int runClass(const std::string& inFile, const std::string& outBase, const std::string& cls) {
    TFile* fin = TFile::Open(inFile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open " << inFile << std::endl;
        return 1;
    }

    TFile* fout = new TFile((outBase + ".root").c_str(), "RECREATE");
    std::ofstream csv(outBase + ".csv");
    // generator settings (PYTHIA version, g) written by gen_ALICE next to the data, kept as '#' header lines
    {
        std::ifstream info(std::string(DATA_DIR) + "/generator_info.txt");
        std::string line;
        bool any = false;
        while (std::getline(info, line)) { csv << "# " << line << "\n"; any = true; }
        if (!any) {
            csv << "# generator_info.txt not found: PYTHIA version and g unknown\n";
            std::cerr << "WARNING: " << DATA_DIR << "/generator_info.txt not found" << std::endl;
        }
        csv << "# class: " << cls << "\n";
        csv << "# deta window (TPC-FMD12): " << WIN_TPC_FMD12.lo << "-" << WIN_TPC_FMD12.hi << "\n";
        csv << "# deta window (TPC-FMD3):  " << WIN_TPC_FMD3.lo  << "-" << WIN_TPC_FMD3.hi  << "\n";
        csv << "# deta window (FMD12-FMD3):" << WIN_FMD12_FMD3.lo << "-" << WIN_FMD12_FMD3.hi << "\n";
    }
    csv << "species,pt_low,pt_high,V2D_TPC_FMD12,err,chi2ndf,V2D_TPC_FMD3,err,chi2ndf,"
           "V2D_FMD12_FMD3,err,chi2ndf,v2sq,err,v2,err\n";

    // ---- reference: FMD12 - FMD3 (one number, no pT dependence) -----
    const double nFF   = count1(fin, Form("hNassoc_%s_FMD12", cls.c_str()));   // FMD1,2 particles = triggers
    const double nFFlm = count1(fin, "hNassoc_LM_FMD12");
    const FitResult fc = templateFit(
        yieldDphi(fin, {Form("hCorr_%s_FMD12_FMD3", cls.c_str())}, {Form("hMix_%s_FMD12_FMD3", cls.c_str())}, nFF,
                  WIN_FMD12_FMD3.lo, WIN_FMD12_FMD3.hi).get(),
        yieldDphi(fin, {"hCorr_LM_FMD12_FMD3"}, {"hMix_LM_FMD12_FMD3"}, nFFlm,
                  WIN_FMD12_FMD3.lo, WIN_FMD12_FMD3.hi).get());
    std::cout << "Class " << cls << ", template LM\n";
    if (fc.ok)
        std::printf("V2D(FMD12-FMD3) = %.4e +- %.4e   chi2/ndf = %.2f\n", fc.V2, fc.err, fc.chi2ndf);
    else
        std::printf("V2D(FMD12-FMD3): fit failed\n");

    TGraphErrors* gRef = new TGraphErrors();
    gRef->SetName("V2D_FMD12_FMD3");
    if (fc.ok) { gRef->SetPoint(0, 0.0, fc.V2); gRef->SetPointError(0, 0.0, fc.err); }
    gRef->Write();

    // ---- identified species -----
    for (int s = 0; s < NSP; ++s) {
        const char* sp = SP_NAME[s];
        TH1D* hN = getH(fin, Form("hNtrig_%s_%s", cls.c_str(), sp));
        if (!hN) continue;
        TH1D* hNlm = getH(fin, Form("hNtrig_LM_%s", sp));
        if (!hNlm) continue;
        const std::vector<std::pair<int, int>> groups = ptGroups(s, hN->GetNbinsX());

        TGraphErrors* gA = new TGraphErrors(); gA->SetName(Form("V2D_TPC_FMD12_%s", sp));
        TGraphErrors* gB = new TGraphErrors(); gB->SetName(Form("V2D_TPC_FMD3_%s",  sp));
        TGraphErrors* gS = new TGraphErrors(); gS->SetName(Form("v2sq_%s", sp));
        // plain v2 = sqrt(v2sq), filled only where v2sq > 0 (v2sq <= 0 has no real sqrt,
        // and is itself a sign that the 3x2PC combination did not work for that point --
        // see the "Sign of V2Delta" note in Readme.md). This is a positive quantity: it
        // carries no information about the sign of the underlying V2D(TPC-FMD12/FMD3).
        TGraphErrors* gV2 = new TGraphErrors(); gV2->SetName(Form("v2_%s", sp));

        std::printf("\n== %s ==\n", sp);
        std::printf("%-11s %-24s %-24s %-24s\n", "pT", "V2D TPC-FMD12", "V2D TPC-FMD3", "v2sq");

        for (const auto& g : groups) {
            const int b0 = g.first, b1 = g.second;
            const double lo = hN->GetBinLowEdge(b0 + 1);
            const double hi = hN->GetBinLowEdge(b1 + 1) + hN->GetBinWidth(b1 + 1);
            const double x  = 0.5 * (lo + hi), ex = 0.5 * (hi - lo);

            double nT = 0, nTlm = 0;
            for (int b = b0; b <= b1; ++b) { nT += hN->GetBinContent(b + 1); nTlm += hNlm->GetBinContent(b + 1); }

            const FitResult fa = templateFit(
                yieldDphi(fin, binNames("hCorr", cls, "TPC_FMD12", sp, b0, b1),
                               binNames("hMix",  cls, "TPC_FMD12", sp, b0, b1), nT,
                               WIN_TPC_FMD12.lo, WIN_TPC_FMD12.hi).get(),
                yieldDphi(fin, binNames("hCorr", "LM", "TPC_FMD12", sp, b0, b1),
                               binNames("hMix",  "LM", "TPC_FMD12", sp, b0, b1), nTlm,
                               WIN_TPC_FMD12.lo, WIN_TPC_FMD12.hi).get());
            const FitResult fb = templateFit(
                yieldDphi(fin, binNames("hCorr", cls, "TPC_FMD3", sp, b0, b1),
                               binNames("hMix",  cls, "TPC_FMD3", sp, b0, b1), nT,
                               WIN_TPC_FMD3.lo, WIN_TPC_FMD3.hi).get(),
                yieldDphi(fin, binNames("hCorr", "LM", "TPC_FMD3", sp, b0, b1),
                               binNames("hMix",  "LM", "TPC_FMD3", sp, b0, b1), nTlm,
                               WIN_TPC_FMD3.lo, WIN_TPC_FMD3.hi).get());

            if (fa.ok) { int n = gA->GetN(); gA->SetPoint(n, x, fa.V2); gA->SetPointError(n, ex, fa.err); }
            if (fb.ok) { int n = gB->GetN(); gB->SetPoint(n, x, fb.V2); gB->SetPointError(n, ex, fb.err); }

            // signed v2^2 = Va*Vb/Vc (no square root taken)
            bool   okC = fa.ok && fb.ok && fc.ok && fa.V2 != 0 && fb.V2 != 0 && fc.V2 != 0;
            double v2sq = 0, ev2sq = 0;
            if (okC) {
                v2sq = fa.V2 * fb.V2 / fc.V2;
                const double ra = fa.err / fa.V2, rb = fb.err / fb.V2, rc = fc.err / fc.V2;
                const double rel = std::sqrt(ra * ra + rb * rb + rc * rc);
                ev2sq = std::fabs(v2sq) * rel;
                int n = gS->GetN(); gS->SetPoint(n, x, v2sq); gS->SetPointError(n, ex, ev2sq);
            }

            // plain v2 = sqrt(v2sq); only defined where v2sq > 0
            const bool   okV2 = okC && v2sq > 0;
            double v2 = 0, ev2 = 0;
            if (okV2) {
                v2  = std::sqrt(v2sq);
                ev2 = ev2sq / (2.0 * v2);
                int n = gV2->GetN(); gV2->SetPoint(n, x, v2); gV2->SetPointError(n, ex, ev2);
            }

            std::printf("%4.2f-%-5.2f  %s  %s  %s\n", lo, hi,
                        fa.ok ? Form("%+.3e +- %.2e", fa.V2, fa.err) : "fit failed             ",
                        fb.ok ? Form("%+.3e +- %.2e", fb.V2, fb.err) : "fit failed             ",
                        okC   ? Form("%+.3e +- %.2e", v2sq, ev2sq)   : "undefined              ");

            csv << sp << "," << lo << "," << hi << ","
                << num(fa.V2, fa.ok) << "," << num(fa.err, fa.ok) << "," << num(fa.chi2ndf, fa.ok) << ","
                << num(fb.V2, fb.ok) << "," << num(fb.err, fb.ok) << "," << num(fb.chi2ndf, fb.ok) << ","
                << num(fc.V2, fc.ok) << "," << num(fc.err, fc.ok) << "," << num(fc.chi2ndf, fc.ok) << ","
                << num(v2sq, okC) << "," << num(ev2sq, okC) << ","
                << num(v2, okV2) << "," << num(ev2, okV2) << "\n";
        }

        gA->Write(); gB->Write(); gS->Write(); gV2->Write();
    }

    csv.close();
    fout->Close();
    fin->Close();
    std::cout << "\nWrote " << outBase << ".root and " << outBase << ".csv\n" << std::endl;
    return 0;
}

int main() {
    TH1::AddDirectory(kFALSE);   // temporary yield histograms stay out of the ROOT directory
    int failed = 0;
    for (const std::string& cls : CLASSES)
        if (runClass(std::string(DATA_DIR) + "/merged.root",
                     std::string(DATA_DIR) + "/v2_results_" + cls, cls) != 0) ++failed;
    return failed ? 1 : 0;
}