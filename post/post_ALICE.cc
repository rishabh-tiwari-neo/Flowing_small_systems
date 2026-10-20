// post_ALICE.cc
// Builds per-trigger yields Y(dphi) = S/B projected over deta, then fits them with the template.
// S and B are 2D (dphi x deta) histograms from gen_ALICE. B is normalised with one constant:
// its dphi average at the deta bin where B is largest. S is divided by the number of triggers
// (for FMD1,2-FMD3: the number of FMD1,2 particles).
// Weighted linear least-squares fitting with analytical covariance-based error estimation.
// The LM template is treated as exact, so its error is not in the V2Δ uncertainty.
// Sparse species (K0S, Lambda) are analysed in coarser pT bins: the 2D S and B histograms of
// several generator pT bins are summed before dividing (see ptGroups).
// The three V2Δ (TPC-FMD1,2, TPC-FMD3, FMD1,2-FMD3) are combined into v2^2 (3x2PC). 

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

// pT bins of the generator that are merged for each species. {first, last} bin, inclusive.
// pi, K, p keep the 20 generator bins; K0S and Lambda use coarser bins because they are sparse.
static const int K0S_IDX = 3, LAM_IDX = 4;
static std::vector<std::pair<int, int>> ptGroups(int species, int nbGen) {
    static const int v0[][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 14}, {15, 16}, {17, 18}, {19, 19}};
    std::vector<std::pair<int, int>> g;
    if ((species == K0S_IDX || species == LAM_IDX) && nbGen == 20)
        for (const auto& r : v0) g.push_back({r[0], r[1]});
    else
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
static std::unique_ptr<TH1D> yieldDphi(TFile* f, const std::vector<TString>& namesS,
                                       const std::vector<TString>& namesB, double ntrig) {
    if (namesS.empty() || namesS.size() != namesB.size() || ntrig <= 0) return nullptr;
    TH2D* h0 = getH2(f, namesS[0]);
    if (!h0) return nullptr;
    const int nx = h0->GetNbinsX(), ny = h0->GetNbinsY();
    const double xmin = h0->GetXaxis()->GetXmin(), xmax = h0->GetXaxis()->GetXmax();

    std::vector<double> S(nx * ny, 0.0), eS2(nx * ny, 0.0), B(nx * ny, 0.0), eB2(nx * ny, 0.0);
    for (size_t k = 0; k < namesS.size(); ++k) {
        TH2D* hS = getH2(f, namesS[k]);
        TH2D* hB = getH2(f, namesB[k]);
        if (!hS || !hB) return nullptr;
        for (int ix = 1; ix <= nx; ++ix)
            for (int iy = 1; iy <= ny; ++iy) {
                const int i = (ix - 1) * ny + (iy - 1);
                S[i] += hS->GetBinContent(ix, iy);  eS2[i] += std::pow(hS->GetBinError(ix, iy), 2);
                B[i] += hB->GetBinContent(ix, iy);  eB2[i] += std::pow(hB->GetBinError(ix, iy), 2);
            }
    }

    double c = 0;
    for (int iy = 1; iy <= ny; ++iy) {
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
    }
    csv << "species,pt_low,pt_high,V2D_TPC_FMD12,err,chi2ndf,V2D_TPC_FMD3,err,chi2ndf,"
           "V2D_FMD12_FMD3,err,chi2ndf,v2sq,err\n";

    // ---- reference: FMD12 - FMD3 (one number, no pT dependence) -----
    const double nFF   = count1(fin, Form("hNassoc_%s_FMD12", cls.c_str()));   // FMD1,2 particles = triggers
    const double nFFlm = count1(fin, "hNassoc_LM_FMD12");
    const FitResult fc = templateFit(
        yieldDphi(fin, {Form("hCorr_%s_FMD12_FMD3", cls.c_str())}, {Form("hMix_%s_FMD12_FMD3", cls.c_str())}, nFF).get(),
        yieldDphi(fin, {"hCorr_LM_FMD12_FMD3"}, {"hMix_LM_FMD12_FMD3"}, nFFlm).get());
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
                               binNames("hMix",  cls, "TPC_FMD12", sp, b0, b1), nT).get(),
                yieldDphi(fin, binNames("hCorr", "LM", "TPC_FMD12", sp, b0, b1),
                               binNames("hMix",  "LM", "TPC_FMD12", sp, b0, b1), nTlm).get());
            const FitResult fb = templateFit(
                yieldDphi(fin, binNames("hCorr", cls, "TPC_FMD3", sp, b0, b1),
                               binNames("hMix",  cls, "TPC_FMD3", sp, b0, b1), nT).get(),
                yieldDphi(fin, binNames("hCorr", "LM", "TPC_FMD3", sp, b0, b1),
                               binNames("hMix",  "LM", "TPC_FMD3", sp, b0, b1), nTlm).get());

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

            std::printf("%4.2f-%-5.2f  %s  %s  %s\n", lo, hi,
                        fa.ok ? Form("%+.3e +- %.2e", fa.V2, fa.err) : "fit failed             ",
                        fb.ok ? Form("%+.3e +- %.2e", fb.V2, fb.err) : "fit failed             ",
                        okC   ? Form("%+.3e +- %.2e", v2sq, ev2sq)   : "undefined              ");

            csv << sp << "," << lo << "," << hi << ","
                << num(fa.V2, fa.ok) << "," << num(fa.err, fa.ok) << "," << num(fa.chi2ndf, fa.ok) << ","
                << num(fb.V2, fb.ok) << "," << num(fb.err, fb.ok) << "," << num(fb.chi2ndf, fb.ok) << ","
                << num(fc.V2, fc.ok) << "," << num(fc.err, fc.ok) << "," << num(fc.chi2ndf, fc.ok) << ","
                << num(v2sq, okC) << "," << num(ev2sq, okC) << "\n";
        }

        gA->Write(); gB->Write(); gS->Write();
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