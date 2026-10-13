// Post-processing for the gen_ALICE output: template fit -> V2Delta -> v2(pT).
//   build:  g++ -O2 -o post/post_ALICE post/post_ALICE.cc `root-config --cflags --libs`
//   run:    ./post/post_ALICE [data/data_ALICE/merged.root] [data/data_ALICE/v2_results]
//
// Per (species, pT bin):
//   V2D(TPC_FMD12), V2D(TPC_FMD3) from the identified-trigger correlations,
//   V2D(FMD12_FMD3) once from the reference-reference correlation,
//   each from the template fit  Y_HM = F*Y_LM + G*[1 + 2 V2 cos(2dphi) + 2 V3 cos(3dphi)]
//   v2 = sqrt( V2D(TPC_FMD12) * V2D(TPC_FMD3) / V2D(FMD12_FMD3) )   (signed, see below)
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"

static const int NSP = 5, NPT = 11;
static const char* SP_NAME[NSP]       = {"pi", "K", "p", "K0S", "Lam"};
static const char* PAIR_NAME[2]       = {"TPC_FMD12", "TPC_FMD3"};
static const double PT_EDGES[NPT + 1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5,
                                         1.75, 2.0, 3.0, 4.0, 5.0, 10.0};

static TFile* gF = nullptr;
static TH1D* get(const std::string& name) {
    TH1D* h = (TH1D*)gF->Get(name.c_str());
    if (!h) fprintf(stderr, "  [missing histogram: %s]\n", name.c_str());
    return h;
}

// ---- 4x4 matrix inverse (Gauss-Jordan, partial pivoting) -------------
static bool invert4(double M[4][4], double C[4][4]) {
    double a[4][8];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) { a[i][j] = M[i][j]; a[i][j + 4] = (i == j) ? 1.0 : 0.0; }
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int r = col + 1; r < 4; ++r) if (std::fabs(a[r][col]) > std::fabs(a[piv][col])) piv = r;
        if (std::fabs(a[piv][col]) < 1e-300) return false;
        if (piv != col) for (int j = 0; j < 8; ++j) std::swap(a[piv][j], a[col][j]);
        const double d = a[col][col];
        for (int j = 0; j < 8; ++j) a[col][j] /= d;
        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            const double f = a[r][col];
            for (int j = 0; j < 8; ++j) a[r][j] -= f * a[col][j];
        }
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) C[i][j] = a[i][j + 4];
    return true;
}

struct Fit { double v2 = 0, err = 0, chi2ndf = 0; bool ok = false; };

// Template fit. The model is linear in (F, G, G*V2, G*V3):
//   Y_HM(x) = F*Y_LM(x) + G + (2 G V2) cos(2x) + (2 G V3) cos(3x)
// so a weighted linear least squares gives the exact minimum and the covariance matrix.
// Both histograms are normalised to unit area first (F and G are free, V2 does not depend on it).
// The statistical error of the LM template is neglected (LM has >100x the statistics).
static Fit templateFit(TH1D* hHM, TH1D* hLM) {
    Fit r;
    if (!hHM || !hLM || hHM->GetNbinsX() != hLM->GetNbinsX()) return r;
    const double sHM = hHM->Integral(), sLM = hLM->Integral();
    if (sHM <= 0 || sLM <= 0) return r;

    double M[4][4] = {{0}}, b[4] = {0};
    int npts = 0;
    for (int i = 1; i <= hHM->GetNbinsX(); ++i) {
        const double y = hHM->GetBinContent(i) / sHM;
        const double e = hHM->GetBinError(i) / sHM;
        if (e <= 0) continue;
        const double x = hHM->GetBinCenter(i);
        const double f[4] = {hLM->GetBinContent(i) / sLM, 1.0, std::cos(2 * x), std::cos(3 * x)};
        const double w = 1.0 / (e * e);
        for (int j = 0; j < 4; ++j) {
            b[j] += w * f[j] * y;
            for (int k = 0; k < 4; ++k) M[j][k] += w * f[j] * f[k];
        }
        ++npts;
    }
    if (npts < 5) return r;

    double C[4][4];
    if (!invert4(M, C)) return r;
    double th[4] = {0, 0, 0, 0};
    for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k) th[j] += C[j][k] * b[k];

    const double G = th[1], a = th[2];
    if (G <= 0) return r;
    r.v2 = a / (2.0 * G);
    const double dA = 1.0 / (2.0 * G), dG = -a / (2.0 * G * G);
    const double var = dA * dA * C[2][2] + dG * dG * C[1][1] + 2.0 * dA * dG * C[2][1];
    if (var < 0) return r;
    r.err = std::sqrt(var);

    double chi2 = 0;
    for (int i = 1; i <= hHM->GetNbinsX(); ++i) {
        const double y = hHM->GetBinContent(i) / sHM;
        const double e = hHM->GetBinError(i) / sHM;
        if (e <= 0) continue;
        const double x = hHM->GetBinCenter(i);
        const double m = th[0] * hLM->GetBinContent(i) / sLM + th[1] + th[2] * std::cos(2 * x) + th[3] * std::cos(3 * x);
        chi2 += (y - m) * (y - m) / (e * e);
    }
    r.chi2ndf = chi2 / (npts - 4);
    r.ok = true;
    return r;
}

int main(int argc, char* argv[]) {
    const std::string inFile = (argc > 1) ? argv[1] : "data/data_ALICE/merged.root";
    const std::string outBase = (argc > 2) ? argv[2] : "data/data_ALICE/v2_results";

    gF = TFile::Open(inFile.c_str());
    if (!gF || gF->IsZombie()) { fprintf(stderr, "Cannot open %s\n", inFile.c_str()); return 1; }

    // ---- reference-reference V2Delta (species and pT independent) ----
    Fit ff = templateFit(get("hCorr_HM_FMD12_FMD3"), get("hCorr_LM_FMD12_FMD3"));
    printf("File: %s\n", inFile.c_str());
    if (!ff.ok) { fprintf(stderr, "FMD12_FMD3 template fit failed\n"); return 1; }
    printf("V2D(FMD12_FMD3) = %.6f +- %.6f   chi2/ndf = %.2f\n\n", ff.v2, ff.err, ff.chi2ndf);

    TFile fout((outBase + ".root").c_str(), "RECREATE");
    FILE* csv = fopen((outBase + ".csv").c_str(), "w");
    fprintf(csv, "species,pt_lo,pt_hi,V2D_tpcfmd12,err,chi2ndf,V2D_tpcfmd3,err,chi2ndf,"
                 "V2D_fmd12fmd3,err,v2sq,v2sq_err,v2,v2_err,valid\n");

    for (int s = 0; s < NSP; ++s) {
        printf("=== %s ===\n", SP_NAME[s]);
        printf("  %-11s %20s %20s %22s %20s\n", "pT (GeV/c)", "V2D TPC-FMD12", "V2D TPC-FMD3",
               "v2^2 (signed)", "v2 (signed sqrt)");

        TGraphErrors gV2, gSq, gA, gB;
        for (int b = 0; b < NPT; ++b) {
            Fit f[2];
            for (int p = 0; p < 2; ++p)
                f[p] = templateFit(get(Form("hCorr_HM_%s_%s_pt%d", PAIR_NAME[p], SP_NAME[s], b)),
                                   get(Form("hCorr_LM_%s_%s_pt%d", PAIR_NAME[p], SP_NAME[s], b)));

            const double x = 0.5 * (PT_EDGES[b] + PT_EDGES[b + 1]);
            const double xe = 0.5 * (PT_EDGES[b + 1] - PT_EDGES[b]);

            double v2sq = 0, v2sqE = 0, v2 = 0, v2E = 0;
            bool valid = f[0].ok && f[1].ok && ff.ok &&
                         f[0].v2 != 0 && f[1].v2 != 0 && ff.v2 != 0;
            if (valid) {
                // errors treated as independent (the two TPC pairings and the FMD pairing share arrays,
                // so this is an approximation)
                v2sq = f[0].v2 * f[1].v2 / ff.v2;
                const double rel = std::sqrt(std::pow(f[0].err / f[0].v2, 2) +
                                             std::pow(f[1].err / f[1].v2, 2) +
                                             std::pow(ff.err / ff.v2, 2));
                v2sqE = std::fabs(v2sq) * rel;
                v2 = (v2sq >= 0 ? 1.0 : -1.0) * std::sqrt(std::fabs(v2sq));   // signed sqrt
                v2E = 0.5 * std::fabs(v2) * rel;
            }

            printf("  %4.2f-%-5.2f  %9.5f+-%-9.5f %9.5f+-%-9.5f %10.3e+-%-10.3e %8.5f+-%-8.5f%s\n",
                   PT_EDGES[b], PT_EDGES[b + 1],
                   f[0].v2, f[0].err, f[1].v2, f[1].err, v2sq, v2sqE, v2, v2E, valid ? "" : "  (invalid)");

            fprintf(csv, "%s,%.2f,%.2f,%.8f,%.8f,%.3f,%.8f,%.8f,%.3f,%.8f,%.8f,%.8e,%.8e,%.8f,%.8f,%d\n",
                    SP_NAME[s], PT_EDGES[b], PT_EDGES[b + 1],
                    f[0].v2, f[0].err, f[0].chi2ndf, f[1].v2, f[1].err, f[1].chi2ndf,
                    ff.v2, ff.err, v2sq, v2sqE, v2, v2E, valid ? 1 : 0);

            if (f[0].ok) { int n = gA.GetN(); gA.SetPoint(n, x, f[0].v2); gA.SetPointError(n, xe, f[0].err); }
            if (f[1].ok) { int n = gB.GetN(); gB.SetPoint(n, x, f[1].v2); gB.SetPointError(n, xe, f[1].err); }
            if (valid) {
                int n = gV2.GetN(); gV2.SetPoint(n, x, v2);   gV2.SetPointError(n, xe, v2E);
                n = gSq.GetN();     gSq.SetPoint(n, x, v2sq); gSq.SetPointError(n, xe, v2sqE);
            }
        }
        printf("\n");

        gA.Write(Form("V2D_TPC_FMD12_%s", SP_NAME[s]));
        gB.Write(Form("V2D_TPC_FMD3_%s",  SP_NAME[s]));
        gV2.Write(Form("v2_%s",   SP_NAME[s]));
        gSq.Write(Form("v2sq_%s", SP_NAME[s]));
    }

    fclose(csv);
    fout.Close();
    gF->Close();
    printf("Wrote %s.root and %s.csv\n", outBase.c_str(), outBase.c_str());
    return 0;
}
