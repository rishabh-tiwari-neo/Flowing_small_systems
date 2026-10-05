// Analyser for the stat_check output (tree "events").
//   build:  make post/analyser          (from ~/project/V2delta/main)
//   run:    ./post/analyser [data/test_run/merged.root]
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

static const int NCLS = 4, NWIN = 5, NSP = 5;
static const char* kClsName[NCLS] = {"ALL", "HM", "LM", "V0M-cut only"};
static const char* kWinName[NWIN] = {"nchCent  (|eta|<0.8, 0.2<pT<3)",
                                     "nchV0A   (2.8<eta<5.1)",
                                     "nchV0C   (-3.7<eta<-1.7)",
                                     "nchFMD12 (1.7<eta<5.1)",
                                     "nchFMD3  (-3.1<eta<-1.7)"};
static const char* kSpName[NSP] = {"pi+-", "K+-", "p+pbar", "K0S", "Lambda+Lbar"};
static const char* kSpCut[NSP]  = {"|eta|<0.8, 0.2<pT<10",
                                   "|eta|<0.8, 0.3<pT<10",
                                   "|eta|<0.8, 0.5<pT<10",
                                   "|y|<0.5, pT<10",
                                   "|y|<0.5, pT<10"};

struct Acc {
    long   nEv = 0;
    double win[NWIN];
    long   sp[NSP];
    Acc() { for (int i = 0; i < NWIN; ++i) win[i] = 0; for (int i = 0; i < NSP; ++i) sp[i] = 0; }
};

int main(int argc, char* argv[]) {
    std::string inFile = (argc > 1) ? argv[1] : "data/test_run/merged.root";

    TFile* f = TFile::Open(inFile.c_str());
    if (!f || f->IsZombie()) { fprintf(stderr, "Cannot open %s\n", inFile.c_str()); return 1; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { fprintf(stderr, "No TTree 'events' in %s\n", inFile.c_str()); return 1; }

    int nchCent = 0, nchV0A = 0, nchV0C = 0, nchFMD12 = 0, nchFMD3 = 0;
    std::vector<float> *pt = nullptr, *eta = nullptr, *rap = nullptr;
    std::vector<int>   *pid = nullptr;

    int bad = 0;
    bad += tree->SetBranchAddress("nchCent",  &nchCent)  < 0;
    bad += tree->SetBranchAddress("nchV0A",   &nchV0A)   < 0;
    bad += tree->SetBranchAddress("nchV0C",   &nchV0C)   < 0;
    bad += tree->SetBranchAddress("nchFMD12", &nchFMD12) < 0;
    bad += tree->SetBranchAddress("nchFMD3",  &nchFMD3)  < 0;
    bad += tree->SetBranchAddress("pt",  &pt)  < 0;
    bad += tree->SetBranchAddress("eta", &eta) < 0;
    bad += tree->SetBranchAddress("y",   &rap) < 0;
    bad += tree->SetBranchAddress("pid", &pid) < 0;
    if (bad) {
        fprintf(stderr, "%d branch(es) missing: file is not from the current stat_check.\n", bad);
        return 1;
    }

    const long nEntries = tree->GetEntries();

    // ---- pass 1: V0M threshold (top 0.07% of events with a V0A and V0C hit) ----
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("nchV0A", 1);
    tree->SetBranchStatus("nchV0C", 1);

    TH1D hV("hV", "V0M proxy", 400, 0, 400);
    for (long i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (nchV0A > 0 && nchV0C > 0) hV.Fill(nchV0A + nchV0C);
    }
    double q = 0.9993, thr = 0;
    hV.GetQuantiles(1, &thr, &q);
    thr = std::round(thr * 10.0) / 10.0;

    // ---- pass 2: everything else ----
    tree->SetBranchStatus("*", 0);
    const char* active[] = {"nchCent", "nchV0A", "nchV0C", "nchFMD12", "nchFMD3",
                            "pt", "eta", "y", "pid"};
    for (const char* b : active) tree->SetBranchStatus(b, 1);

    Acc acc[NCLS];

    for (long i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        const bool coinc = (nchV0A > 0 && nchV0C > 0);
        const int  v0m   = nchV0A + nchV0C;
        const bool v0cut = (v0m > thr);
        const bool hm    = v0cut && nchCent > 25;   // no coincidence requirement, as in the ROOT-prompt numbers
        const bool lm    = coinc && nchCent < 20;

        long sp[NSP] = {0, 0, 0, 0, 0};
        const size_t n = pid->size();
        for (size_t j = 0; j < n; ++j) {
            const int    a  = std::abs((*pid)[j]);
            const double e  = std::fabs((*eta)[j]);
            const double yy = std::fabs((*rap)[j]);
            const double p  = (*pt)[j];
            if (a == 211  && e < 0.8 && p > 0.2 && p < 10.0) ++sp[0];
            if (a == 321  && e < 0.8 && p > 0.3 && p < 10.0) ++sp[1];
            if (a == 2212 && e < 0.8 && p > 0.5 && p < 10.0) ++sp[2];
            if (a == 310  && yy < 0.5 && p < 10.0)           ++sp[3];
            if (a == 3122 && yy < 0.5 && p < 10.0)           ++sp[4];
        }

        const double w[NWIN] = {(double)nchCent, (double)nchV0A, (double)nchV0C,
                                (double)nchFMD12, (double)nchFMD3};
        const bool in[NCLS] = {true, hm, lm, v0cut};
        for (int c = 0; c < NCLS; ++c) {
            if (!in[c]) continue;
            ++acc[c].nEv;
            for (int k = 0; k < NWIN; ++k) acc[c].win[k] += w[k];
            for (int k = 0; k < NSP;  ++k) acc[c].sp[k]  += sp[k];
        }
    }

    // ---- report ----
    printf("File: %s\n", inFile.c_str());
    printf("V0M proxy threshold (top 0.07%% of events with V0A and V0C hit): > %.1f\n\n", thr);

    printf("=== Events ===\n");
    for (int c = 0; c < NCLS; ++c)
        printf("  %-14s %10ld   (%.4f%% of all)\n", kClsName[c], acc[c].nEv,
               100.0 * acc[c].nEv / (double)nEntries);
    printf("  %-14s %10ld\n", "Neither HM nor LM", acc[0].nEv - acc[1].nEv - acc[2].nEv);
    printf("  HM  = (nchV0A+nchV0C) > thr && nchCent > 25\n");
    printf("  LM  = nchV0A>0 && nchV0C>0 && nchCent < 20\n\n");

    printf("=== Charged particles per window: mean per event / total ===\n");
    for (int k = 0; k < NWIN; ++k) {
        printf("  %s\n", kWinName[k]);
        for (int c = 0; c < 3; ++c) {
            const double mean = acc[c].nEv ? acc[c].win[k] / acc[c].nEv : 0.0;
            printf("     %-4s mean = %8.2f   total = %.0f\n", kClsName[c], mean, acc[c].win[k]);
        }
    }
    printf("  (V0A and FMD1,2 overlap in 2.8<eta<5.1: do not add them)\n\n");

    printf("=== Identified particles: total / per event ===\n");
    for (int k = 0; k < NSP; ++k) {
        printf("  %-12s [%s]\n", kSpName[k], kSpCut[k]);
        for (int c = 0; c < 3; ++c) {
            const double pe = acc[c].nEv ? (double)acc[c].sp[k] / acc[c].nEv : 0.0;
            printf("     %-4s total = %10ld   per event = %7.3f\n", kClsName[c], acc[c].sp[k], pe);
        }
    }
    printf("\n");

    printf("=== Ratios (from the totals above; cuts differ between rows) ===\n");
    printf("  %-5s %8s %8s %8s %10s %8s\n", "", "K/pi", "p/pi", "K0S/K", "Lam/K0S", "Lam/p");
    for (int c = 0; c < 3; ++c) {
        auto r = [&](int a, int b) { return acc[c].sp[b] ? (double)acc[c].sp[a] / acc[c].sp[b] : 0.0; };
        printf("  %-5s %8.4f %8.4f %8.4f %10.4f %8.4f\n", kClsName[c],
               r(1, 0), r(2, 0), r(3, 1), r(4, 3), r(4, 2));
    }

    f->Close();
    return 0;
}
