// Analyser for the gen_ALICE output (histograms only, no tree).
//   build:  g++ -O2 -o post/analyser_ALICE post/analyser_ALICE.cc `root-config --cflags --libs`
//   run:    ./post/analyser_ALICE [data/data_ALICE/merged.root]
#include <cstdio>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TH1D.h"

static const int NCLS = 2, NSP = 5, NPT = 11, NIDPAIR = 2;
static const char* CLS_NAME[NCLS]  = {"HM", "LM"};
static const char* SP_NAME[NSP]    = {"pi", "K", "p", "K0S", "Lam"};
static const char* PAIR_NAME[NIDPAIR] = {"TPC_FMD12", "TPC_FMD3"};
static const double PT_EDGES[NPT + 1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5,
                                         1.75, 2.0, 3.0, 4.0, 5.0, 10.0};
static const int    V0M_HM_MIN = 157;     // must match gen_ALICE.cc
static const double LOW_STAT   = 1000.0;  // flag bins with fewer triggers than this

static TFile* gF = nullptr;
static TH1D* get(const std::string& name) {
    TH1D* h = (TH1D*)gF->Get(name.c_str());
    if (!h) fprintf(stderr, "  [missing histogram: %s]\n", name.c_str());
    return h;
}

int main(int argc, char* argv[]) {
    std::string inFile = (argc > 1) ? argv[1] : "data/data_ALICE/merged.root";
    gF = TFile::Open(inFile.c_str());
    if (!gF || gF->IsZombie()) { fprintf(stderr, "Cannot open %s\n", inFile.c_str()); return 1; }

    TH1D* hCounter = get("hCounter");
    TH1D* hNchCent = get("hNchCent");
    TH1D* hNchV0M  = get("hNchV0M");
    if (!hCounter) return 1;

    const double nGen = hCounter->GetBinContent(1);
    const double nReg = hCounter->GetBinContent(2);
    const double nSel = hCounter->GetBinContent(3);

    printf("File: %s\n\n", inFile.c_str());

    // ---- events ----
    printf("=== Events ===\n");
    printf("  generated             %12.0f\n", nGen);
    printf("  registered (V0A&&V0C) %12.0f   (%.4f%% of generated)\n", nReg, 100.0 * nReg / nGen);
    printf("  HM or LM              %12.0f   (%.4f%% of registered)\n", nSel, 100.0 * nSel / nReg);
    double nEv[NCLS] = {0, 0};
    for (int c = 0; c < NCLS; ++c) {
        TH1D* h = get(Form("hEvents_%s", CLS_NAME[c]));
        nEv[c] = h ? h->GetBinContent(1) : 0;
        printf("  %-3s                   %12.0f   (%.4f%% of registered)\n",
               CLS_NAME[c], nEv[c], 100.0 * nEv[c] / nReg);
    }
    printf("  neither HM nor LM     %12.0f\n\n", nReg - nEv[0] - nEv[1]);

    // ---- inclusive multiplicity of registered events ----
    printf("=== Inclusive (registered events) ===\n");
    if (hNchCent) printf("  <nchCent> = %.2f\n", hNchCent->GetMean());
    if (hNchV0M) {
        printf("  <V0M>     = %.2f\n", hNchV0M->GetMean());
        // value v sits in bin v+1 (bins centred on integers)
        const double above = hNchV0M->Integral(V0M_HM_MIN + 2, hNchV0M->GetNbinsX());
        printf("  fraction with V0M > %d: %.4f%%  (target 0.07%%)\n", V0M_HM_MIN, 100.0 * above / nReg);
    }
    printf("  (per-class <nchCent> is not stored in gen_ALICE: add hNchCent_HM / hNchCent_LM to get it)\n\n");

    // ---- reference arms ----
    printf("=== Reference particles per event (FMD, all charged, no PID) ===\n");
    for (int c = 0; c < NCLS; ++c) {
        TH1D* h12 = get(Form("hNassoc_%s_FMD12", CLS_NAME[c]));
        TH1D* h3  = get(Form("hNassoc_%s_FMD3",  CLS_NAME[c]));
        if (!h12 || !h3 || nEv[c] == 0) continue;
        printf("  %-3s  FMD1,2: %8.2f   FMD3: %8.2f\n", CLS_NAME[c],
               h12->GetBinContent(1) / nEv[c], h3->GetBinContent(1) / nEv[c]);
    }
    printf("\n");

    // ---- identified triggers ----
    printf("=== Identified triggers: totals and per event (all pT bins) ===\n");
    for (int c = 0; c < NCLS; ++c) {
        printf("  %s\n", CLS_NAME[c]);
        for (int s = 0; s < NSP; ++s) {
            TH1D* h = get(Form("hNtrig_%s_%s", CLS_NAME[c], SP_NAME[s]));
            if (!h) continue;
            const double tot = h->Integral();
            printf("     %-4s total = %12.0f   per event = %8.4f\n", SP_NAME[s], tot,
                   nEv[c] ? tot / nEv[c] : 0.0);
        }
    }
    printf("\n");

    // ---- N_trig per pT bin ----
    printf("=== N_trig per pT bin (* = fewer than %.0f) ===\n", LOW_STAT);
    for (int c = 0; c < NCLS; ++c) {
        printf("  %s\n", CLS_NAME[c]);
        printf("     %-12s", "pT bin");
        for (int s = 0; s < NSP; ++s) printf(" %12s", SP_NAME[s]);
        printf("\n");
        for (int b = 0; b < NPT; ++b) {
            printf("     %5.2f-%-5.2f ", PT_EDGES[b], PT_EDGES[b + 1]);
            for (int s = 0; s < NSP; ++s) {
                TH1D* h = get(Form("hNtrig_%s_%s", CLS_NAME[c], SP_NAME[s]));
                const double n = h ? h->GetBinContent(b + 1) : 0.0;
                printf(" %11.0f%c", n, n < LOW_STAT ? '*' : ' ');
            }
            printf("\n");
        }
    }
    printf("\n");

    // ---- pairs per pairing ----
    printf("=== Pairs filled (sum over pT bins) and pairs per trigger ===\n");
    for (int c = 0; c < NCLS; ++c) {
        printf("  %s\n", CLS_NAME[c]);
        for (int p = 0; p < NIDPAIR; ++p)
            for (int s = 0; s < NSP; ++s) {
                double pairs = 0, trig = 0;
                for (int b = 0; b < NPT; ++b) {
                    TH1D* h = get(Form("hCorr_%s_%s_%s_pt%d", CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], b));
                    if (h) pairs += h->Integral();
                }
                TH1D* ht = get(Form("hNtrig_%s_%s", CLS_NAME[c], SP_NAME[s]));
                if (ht) trig = ht->Integral();
                printf("     %-9s %-4s pairs = %14.0f   per trigger = %8.2f\n",
                       PAIR_NAME[p], SP_NAME[s], pairs, trig ? pairs / trig : 0.0);
            }
        TH1D* hff = get(Form("hCorr_%s_FMD12_FMD3", CLS_NAME[c]));
        if (hff) printf("     FMD12_FMD3 pairs = %14.0f   per event = %10.1f\n",
                        hff->Integral(), nEv[c] ? hff->Integral() / nEv[c] : 0.0);
    }

    gF->Close();
    return 0;
}
