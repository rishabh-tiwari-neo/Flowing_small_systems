// post_ALICE.cc
// weighted linear least-squares fitting with analytical covariance-based error estimation.
// The LM template is treated as exact, so its error is not in the V2Δ uncertainty.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"

static const int NSP = 5;
static const char* SP_NAME[NSP] = {"pi", "K", "p", "K0S", "Lam"};

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
    const double is = hs->Integral(), it = ht->Integral();
    if (is <= 0 || it <= 0) return res;

    // both histograms normalised to unit area 
    std::vector<double> y, e, t, ph;
    for (int i = 1; i <= hs->GetNbinsX(); ++i) {
        const double err = hs->GetBinError(i) / is;
        if (err <= 0) continue;
        y.push_back(hs->GetBinContent(i) / is);
        e.push_back(err);
        t.push_back(ht->GetBinContent(i) / it);
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

static std::string num(double v, bool ok) {
    if (!ok) return "nan";
    char b[32];
    std::snprintf(b, sizeof b, "%.5e", v);
    return b;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " merged.root out_base [HM|HM35]" << std::endl;
        return 1;
    }
    const std::string inFile  = argv[1];
    const std::string outBase = argv[2];
    const std::string cls     = (argc > 3) ? argv[3] : "HM";

    TFile* fin = TFile::Open(inFile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open " << inFile << std::endl;
        return 1;
    }

    TFile* fout = new TFile((outBase + ".root").c_str(), "RECREATE");
    std::ofstream csv(outBase + ".csv");
    csv << "species,pt_low,pt_high,V2D_TPC_FMD12,err,chi2ndf,V2D_TPC_FMD3,err,chi2ndf,"
           "V2D_FMD12_FMD3,err,chi2ndf,v2sq,err,v2,err\n";

    // ---- reference: FMD12 - FMD3 (one number, no pT dependence) -----
    const FitResult fc = templateFit(getH(fin, Form("hCorr_%s_FMD12_FMD3", cls.c_str())),
                                     getH(fin, "hCorr_LM_FMD12_FMD3"));
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
        const int nb = hN->GetNbinsX();

        TGraphErrors* gA = new TGraphErrors(); gA->SetName(Form("V2D_TPC_FMD12_%s", sp));
        TGraphErrors* gB = new TGraphErrors(); gB->SetName(Form("V2D_TPC_FMD3_%s",  sp));
        TGraphErrors* gS = new TGraphErrors(); gS->SetName(Form("v2sq_%s", sp));
        TGraphErrors* gV = new TGraphErrors(); gV->SetName(Form("v2_%s",   sp));

        std::printf("\n== %s ==\n", sp);
        std::printf("%-11s %-24s %-24s %-24s\n", "pT", "V2D TPC-FMD12", "V2D TPC-FMD3", "v2");

        for (int b = 0; b < nb; ++b) {
            const double lo = hN->GetBinLowEdge(b + 1);
            const double hi = lo + hN->GetBinWidth(b + 1);
            const double x  = 0.5 * (lo + hi), ex = 0.5 * (hi - lo);

            const FitResult fa = templateFit(
                getH(fin, Form("hCorr_%s_TPC_FMD12_%s_pt%d", cls.c_str(), sp, b)),
                getH(fin, Form("hCorr_LM_TPC_FMD12_%s_pt%d", sp, b)));
            const FitResult fb = templateFit(
                getH(fin, Form("hCorr_%s_TPC_FMD3_%s_pt%d", cls.c_str(), sp, b)),
                getH(fin, Form("hCorr_LM_TPC_FMD3_%s_pt%d", sp, b)));

            if (fa.ok) { int n = gA->GetN(); gA->SetPoint(n, x, fa.V2); gA->SetPointError(n, ex, fa.err); }
            if (fb.ok) { int n = gB->GetN(); gB->SetPoint(n, x, fb.V2); gB->SetPointError(n, ex, fb.err); }

            // signed v2^2 = Va*Vb/Vc ; v2 = sign*sqrt(|v2^2|)
            bool   okC = fa.ok && fb.ok && fc.ok && fa.V2 != 0 && fb.V2 != 0 && fc.V2 != 0;
            double v2sq = 0, ev2sq = 0, v2 = 0, ev2 = 0;
            if (okC) {
                v2sq = fa.V2 * fb.V2 / fc.V2;
                const double ra = fa.err / fa.V2, rb = fb.err / fb.V2, rc = fc.err / fc.V2;
                const double rel = std::sqrt(ra * ra + rb * rb + rc * rc);
                ev2sq = std::fabs(v2sq) * rel;
                v2    = (v2sq >= 0 ? 1.0 : -1.0) * std::sqrt(std::fabs(v2sq));
                ev2   = 0.5 * std::fabs(v2) * rel;
                int n = gS->GetN(); gS->SetPoint(n, x, v2sq); gS->SetPointError(n, ex, ev2sq);
                n = gV->GetN();     gV->SetPoint(n, x, v2);   gV->SetPointError(n, ex, ev2);
            }

            std::printf("%4.2f-%-5.2f  %s  %s  %s\n", lo, hi,
                        fa.ok ? Form("%+.3e +- %.2e", fa.V2, fa.err) : "fit failed             ",
                        fb.ok ? Form("%+.3e +- %.2e", fb.V2, fb.err) : "fit failed             ",
                        okC   ? Form("%+.3e +- %.2e", v2, ev2)       : "undefined              ");

            csv << sp << "," << lo << "," << hi << ","
                << num(fa.V2, fa.ok) << "," << num(fa.err, fa.ok) << "," << num(fa.chi2ndf, fa.ok) << ","
                << num(fb.V2, fb.ok) << "," << num(fb.err, fb.ok) << "," << num(fb.chi2ndf, fb.ok) << ","
                << num(fc.V2, fc.ok) << "," << num(fc.err, fc.ok) << "," << num(fc.chi2ndf, fc.ok) << ","
                << num(v2sq, okC) << "," << num(ev2sq, okC) << ","
                << num(v2, okC)   << "," << num(ev2, okC) << "\n";
        }

        gA->Write(); gB->Write(); gS->Write(); gV->Write();
    }

    csv.close();
    fout->Close();
    fin->Close();
    std::cout << "\nWrote " << outBase << ".root and " << outBase << ".csv" << std::endl;
    return 0;
}
