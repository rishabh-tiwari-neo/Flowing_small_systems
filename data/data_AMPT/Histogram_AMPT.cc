// histogram_AMPT.cc
// Turns one AMPT tree file (tree "t", from AmptRootIO.cxx) into the same histograms that
// gen_ALICE.cc writes, with identical names, binning, acceptance windows, class logic,
// pairing and mixed-event pool, so post_ALICE and plot_ALICE read the output unchanged.
//
// Usage:
//   histogram_AMPT <input_tree.root> <output.root> [V0M_HM_MIN] [V0M_HM35_MIN] [NCH_HM_MIN] [NCH_LM_MAX] [input.ampt]
// Defaults (AMPT, from the 1.5M-event scan of the current input.ampt):
//   V0M_HM_MIN = 142, V0M_HM35_MIN = 72, NCH_HM_MIN = 25, NCH_LM_MAX = 20
// (gen_ALICE/PYTHIA uses 157, 45, 25, 20.)
// A class cut is "V0M > V0M_*_MIN" and "nchCent > NCH_HM_MIN", exactly as in gen_ALICE.
// Mixed events are formed within this one input file (one shard).
// Writes generator_info.txt next to the output file.
//
// Build: g++ -O2 -o histogram_AMPT histogram_AMPT.cc `root-config --cflags --libs`

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

// ---- acceptance windows (same as gen_ALICE.cc) -----
static const double FMD12_ETA_MIN = 1.7,  FMD12_ETA_MAX = 5.1;
static const double FMD3_ETA_MIN  = -3.1, FMD3_ETA_MAX  = -1.7;

static const double NCH_ETA_MAX = 0.8;
static const double NCH_PT_MIN  = 0.2, NCH_PT_MAX = 3.0;

static const double V0A_ETA_MIN = 2.8,  V0A_ETA_MAX = 5.1;
static const double V0C_ETA_MIN = -3.7, V0C_ETA_MAX = -1.7;

// ---- class cuts: set from the command line (defaults for AMPT) ----
static int V0M_HM_MIN   = 142;   // HM:   (nchV0A + nchV0C) > V0M_HM_MIN   (PYTHIA: 157)
static int V0M_HM35_MIN = 72;    // HM35: (nchV0A + nchV0C) > V0M_HM35_MIN (PYTHIA: 45)
static int NCH_HM_MIN   = 25;    // HM and HM35: nchCent > NCH_HM_MIN
static int NCH_LM_MAX   = 20;    // LM: nchCent < NCH_LM_MAX

// ---- identified-particle acceptance ----
static const double PID_ETA_MAX = 0.8;
static const double V0_Y_MAX    = 0.5;

// ---- common pT binning for all species ----
static const int    NpT_bins = 20;
static const double PT_EDGES[NpT_bins + 1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0,
                                         2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 10.0};

inline int ptBin(double pt) {
    if (pt < PT_EDGES[0] || pt >= PT_EDGES[NpT_bins]) return -1;
    return int(std::upper_bound(PT_EDGES, PT_EDGES + NpT_bins + 1, pt) - PT_EDGES) - 1;
}

// ---- event classes ----
enum EvClass { HM = 0, LM = 1, HM35 = 2, NCLS = 3 };
static const char* CLS_NAME[NCLS] = {"HM", "LM", "HM35"};

inline void classify(int nchV0A, int nchV0C, int nchCent, bool Ecl[NCLS]) {
    for (int c = 0; c < NCLS; ++c) Ecl[c] = false;
    if (!(nchV0A > 0 && nchV0C > 0)) return;                    // MB trigger
    const int v0m = nchV0A + nchV0C;
    Ecl[HM]   = (v0m > V0M_HM_MIN   && nchCent > NCH_HM_MIN);
    Ecl[HM35] = (v0m > V0M_HM35_MIN && nchCent > NCH_HM_MIN);
    Ecl[LM]   = (nchCent < NCH_LM_MAX);
}

// ---- identified species ----
enum Species { PI = 0, K, P, K0S, LAM, NSP };
static const char* SP_NAME[NSP] = {"pi", "K", "p", "K0S", "Lam"};

inline int speciesOf(int absId) {
    switch (absId) {
        case 211:  return PI;
        case 321:  return K;
        case 2212: return P;
        case 310:  return K0S;
        case 3122: return LAM;
        default:   return -1;
    }
}

// ---- pairings ----
enum Pairing { TPC_FMD12 = 0, TPC_FMD3 = 1, FMD12_FMD3 = 2, NPAIR = 3 };
static const int NIDPAIR = 2;
static const char* PAIR_NAME[NPAIR] = {"TPC_FMD12", "TPC_FMD3", "FMD12_FMD3"};

static const double DETA_MIN = 1.1;
static const double DETA_MAX = 6.4;
static const int NDETA = 53;

static const int POOL_DEPTH = 5;

static const int    NDPHI   = 20;
static const double DPHI_LO = -M_PI / 2.0;
static const double DPHI_HI =  3.0 * M_PI / 2.0;

inline double wrapDPhi(double d) {
    while (d <  DPHI_LO) d += 2.0 * M_PI;
    while (d >= DPHI_HI) d -= 2.0 * M_PI;
    return d;
}

struct Trig { double eta, phi; int sp, pb; };
struct Ref  { double eta, phi; };
struct PooledEvent { std::vector<Ref> f12, f3; };

// tree buffers (kMaxTracks of AmptRootIO.cxx)
static const int kMaxTracks = 70000;
static int   gNtrack;
static int   gPID[kMaxTracks], gCharge[kMaxTracks];
static float gPx[kMaxTracks], gPy[kMaxTracks], gPz[kMaxTracks], gE[kMaxTracks];

int main(int argc, char* argv[]) {

    if (argc < 3 || argc > 9) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_tree.root> <output.root> [V0M_HM_MIN] [V0M_HM35_MIN] [NCH_HM_MIN] [NCH_LM_MAX] [input.ampt] [chunk_events]"
                  << std::endl;
        return 1;
    }
    long chunk = 1000000;                       // events per block (0 = whole file in one block)
    if (argc > 8) chunk = std::atol(argv[8]);
    const std::string inFile  = argv[1];
    const std::string outFile = argv[2];
    if (argc > 3) V0M_HM_MIN   = std::atoi(argv[3]);
    if (argc > 4) V0M_HM35_MIN = std::atoi(argv[4]);
    if (argc > 5) NCH_HM_MIN   = std::atoi(argv[5]);
    if (argc > 6) NCH_LM_MAX   = std::atoi(argv[6]);
    const std::string inputAmpt = (argc > 7) ? argv[7] : "";

    TH1::AddDirectory(kFALSE);   // histograms are written explicitly, not owned by the input file

    // ---- generator_info.txt next to the output file (same content for every shard) ----
    {
        const size_t slash = outFile.find_last_of('/');
        const std::string dir = (slash == std::string::npos) ? "." : outFile.substr(0, slash);
        std::ofstream info(dir + "/generator_info.txt", std::ios::trunc);
        info << "Generator: AMPT (converted with histogram_AMPT.cc)\n"
             << "Class cuts: HM: V0M > " << V0M_HM_MIN << " and nchCent > " << NCH_HM_MIN
             << "; HM35: V0M > " << V0M_HM35_MIN << " and nchCent > " << NCH_HM_MIN
             << "; LM: nchCent < " << NCH_LM_MAX << "\n"
             << "(PYTHIA/gen_ALICE cuts: HM V0M > 157, HM35 V0M > 45, nchCent > 25, LM nchCent < 20)\n"
             << "No shoving strength g: AMPT has none, so no TParameter g is written by this program.\n";
        if (!inputAmpt.empty()) {
            std::ifstream in(inputAmpt);
            if (!in) {
                info << "input.ampt: could not open " << inputAmpt << "\n";
            } else {
                info << "input.ampt (" << inputAmpt << "):\n";
                std::string line;
                while (std::getline(in, line)) info << "  " << line << "\n";
            }
        } else {
            info << "input.ampt: not given (pass its path as the 7th argument to log the AMPT settings)\n";
        }
    }

    // ---- input tree ----
    TFile* fin = TFile::Open(inFile.c_str(), "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open " << inFile << std::endl; return 1; }
    TTree* tree = dynamic_cast<TTree*>(fin->Get("t"));
    if (!tree) { std::cerr << "No tree 't' in " << inFile << std::endl; return 1; }

    tree->SetBranchStatus("*", 0);
    for (const char* b : {"ntrack", "PID", "charge", "px", "py", "pz", "energy"}) tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("ntrack", &gNtrack);
    tree->SetBranchAddress("PID",    gPID);
    tree->SetBranchAddress("charge", gCharge);
    tree->SetBranchAddress("px",     gPx);
    tree->SetBranchAddress("py",     gPy);
    tree->SetBranchAddress("pz",     gPz);
    tree->SetBranchAddress("energy", gE);
    const long n_total = tree->GetEntries();
    if (chunk <= 0 || chunk > n_total) chunk = n_total;
    const long nChunks = (n_total + chunk - 1) / chunk;

    // events AMPT did not write (no final-state particles) still count as generated
    long long nNotWritten = 0;
    if (TNamed* te = dynamic_cast<TNamed*>(fin->Get("tEvents"))) {
        long long a = 0, b = 0;
        if (std::sscanf(te->GetTitle(), "written %lld, not written: no particles %lld", &a, &b) == 2) nNotWritten = b;
    }

    // ---- one block of 'chunk' events at a time: own histograms, own mixed-event pool, own output file ----
    for (long part = 0; part < nChunks; ++part) {
    const long first = part * chunk;
    const long last  = std::min(first + chunk, n_total);

    // ---- histograms and counters (same names and binning as gen_ALICE.cc) ----
    TH2D* hID[NCLS][NSP][NpT_bins][NIDPAIR];
    TH2D* hIDmix[NCLS][NSP][NpT_bins][NIDPAIR];
    TH2D* hFF[NCLS];
    TH2D* hFFmix[NCLS];
    TH1D* hNtrig[NCLS][NSP];
    TH1D* hNassoc[NCLS][2];
    TH1D* hEvents[NCLS];
    TH1D* hNchClass[NCLS];
    TH1D* hCounter = new TH1D("hCounter", "event counter;;events", 3, 0.5, 3.5);
    hCounter->GetXaxis()->SetBinLabel(1, "generated");
    hCounter->GetXaxis()->SetBinLabel(2, "registered (V0A&&V0C)");
    hCounter->GetXaxis()->SetBinLabel(3, "HM or LM or HM35");

    TH1D* hNchCent = new TH1D("hNchCent", "N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0);N_{ch};events",
                              301, -0.5, 300.5);
    TH1D* hNchV0M  = new TH1D("hNchV0M", "N_{ch} V0A+V0C;N_{ch};events", 401, -0.5, 400.5);

    for (int c = 0; c < NCLS; ++c) {
        for (int s = 0; s < NSP; ++s) {
            hNtrig[c][s] = new TH1D(Form("hNtrig_%s_%s", CLS_NAME[c], SP_NAME[s]),
                                    ";p_{T} (GeV/c);N_{trig}", NpT_bins, PT_EDGES);
            for (int b = 0; b < NpT_bins; ++b)
                for (int p = 0; p < NIDPAIR; ++p) {
                    hID[c][s][b][p] = new TH2D(
                        Form("hCorr_%s_%s_%s_pt%d", CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], b),
                        Form("%s %s %s, p_{T} %.2f-%.2f;#Delta#phi;#Delta#eta;pairs",
                             CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], PT_EDGES[b], PT_EDGES[b + 1]),
                        NDPHI, DPHI_LO, DPHI_HI, NDETA, DETA_MIN, DETA_MAX);
                    hIDmix[c][s][b][p] = new TH2D(
                        Form("hMix_%s_%s_%s_pt%d", CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], b),
                        Form("%s %s %s mixed, p_{T} %.2f-%.2f;#Delta#phi;#Delta#eta;pairs",
                             CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], PT_EDGES[b], PT_EDGES[b + 1]),
                        NDPHI, DPHI_LO, DPHI_HI, NDETA, DETA_MIN, DETA_MAX);
                }
        }

        hFF[c] = new TH2D(Form("hCorr_%s_FMD12_FMD3", CLS_NAME[c]),
                          Form("%s FMD12_FMD3;#Delta#phi;#Delta#eta;pairs", CLS_NAME[c]),
                          NDPHI, DPHI_LO, DPHI_HI, NDETA, DETA_MIN, DETA_MAX);
        hFFmix[c] = new TH2D(Form("hMix_%s_FMD12_FMD3", CLS_NAME[c]),
                             Form("%s FMD12_FMD3 mixed;#Delta#phi;#Delta#eta;pairs", CLS_NAME[c]),
                             NDPHI, DPHI_LO, DPHI_HI, NDETA, DETA_MIN, DETA_MAX);
        hNassoc[c][0] = new TH1D(Form("hNassoc_%s_FMD12", CLS_NAME[c]), "N_{assoc} FMD12", 1, 0, 1);
        hNassoc[c][1] = new TH1D(Form("hNassoc_%s_FMD3",  CLS_NAME[c]), "N_{assoc} FMD3",  1, 0, 1);
        hEvents[c]    = new TH1D(Form("hEvents_%s", CLS_NAME[c]), "events", 1, 0, 1);
        hNchClass[c]  = new TH1D(Form("hNchCent_%s", CLS_NAME[c]),
                                 Form("N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0), %s;N_{ch};events", CLS_NAME[c]),
                                 301, -0.5, 300.5);
    }

    // events AMPT did not write count as generated (never registered)
    if (part == 0) for (long long k = 0; k < nNotWritten; ++k) hCounter->Fill(1);

    // ---- event loop ----
    std::vector<Trig> trig;
    std::vector<Ref>  f12, f3;
    std::deque<PooledEvent> pool[NCLS];

    for (long i = first; i < last; ++i) {

        tree->GetEntry(i);
        hCounter->Fill(1);

        trig.clear(); f12.clear(); f3.clear();

        int nchCent = 0, nchV0A = 0, nchV0C = 0;
        const int entries = std::min(gNtrack, kMaxTracks);
        for (int j = 0; j < entries; ++j) {
            const double px = gPx[j], py = gPy[j], pz = gPz[j], E = gE[j];
            const double pt = std::sqrt(px * px + py * py);
            const double p  = std::sqrt(pt * pt + pz * pz);
            if (p == std::fabs(pz)) continue;                 // eta undefined (pT = 0)
            const double eta = 0.5 * std::log((p + pz) / (p - pz));
            const double phi = std::atan2(py, px);

            if (gCharge[j] != 0) {
                if (std::fabs(eta) < NCH_ETA_MAX && pt > NCH_PT_MIN && pt < NCH_PT_MAX) ++nchCent;
                if (eta > V0A_ETA_MIN && eta < V0A_ETA_MAX) ++nchV0A;
                if (eta > V0C_ETA_MIN && eta < V0C_ETA_MAX) ++nchV0C;
                if (eta > FMD12_ETA_MIN && eta < FMD12_ETA_MAX) f12.push_back({eta, phi});
                if (eta > FMD3_ETA_MIN  && eta < FMD3_ETA_MAX)  f3.push_back({eta, phi});
            }

            // identified triggers
            const int sp = speciesOf(std::abs(gPID[j]));
            if (sp < 0) continue;
            const bool isV0 = (sp == K0S || sp == LAM);
            bool inAcc;
            if (isV0) {
                if (E == std::fabs(pz)) continue;             // y undefined
                const double y = 0.5 * std::log((E + pz) / (E - pz));
                inAcc = (std::fabs(y) < V0_Y_MAX && std::fabs(eta) < PID_ETA_MAX);
            } else {
                inAcc = (std::fabs(eta) < PID_ETA_MAX);
            }
            const int pb = ptBin(pt);
            if (inAcc && pb >= 0) trig.push_back({eta, phi, sp, pb});
        }

        // MB trigger: V0A && V0C coincidence
        if (nchV0A > 0 && nchV0C > 0) {
            hCounter->Fill(2);
            hNchCent->Fill(nchCent);
            hNchV0M->Fill(nchV0A + nchV0C);
        }

        bool Ecl[NCLS];
        classify(nchV0A, nchV0C, nchCent, Ecl);
        if (Ecl[HM] || Ecl[LM] || Ecl[HM35]) hCounter->Fill(3);
        for (int c = 0; c < NCLS; ++c) {
            if (!Ecl[c]) continue;

            hEvents[c]->Fill(0.5);
            hNchClass[c]->Fill(nchCent);
            hNassoc[c][0]->Fill(0.5, (double)f12.size());
            hNassoc[c][1]->Fill(0.5, (double)f3.size());

            // same-event pairing
            for (const Trig& t : trig) {
                hNtrig[c][t.sp]->Fill(0.5 * (PT_EDGES[t.pb] + PT_EDGES[t.pb + 1]));

                for (const Ref& r : f12) {                       // TPC - FMD1,2
                    const double deta = r.eta - t.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hID[c][t.sp][t.pb][TPC_FMD12]->Fill(wrapDPhi(t.phi - r.phi), deta);
                }
                for (const Ref& r : f3) {                        // TPC - FMD3
                    const double deta = t.eta - r.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hID[c][t.sp][t.pb][TPC_FMD3]->Fill(wrapDPhi(t.phi - r.phi), deta);
                }
            }
            for (const Ref& a : f12)                             // FMD1,2 - FMD3
                for (const Ref& b : f3) {
                    const double deta = a.eta - b.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hFF[c]->Fill(wrapDPhi(a.phi - b.phi), deta);
                }

            // mixed-event pairing with the pooled events of this class
            for (const PooledEvent& m : pool[c]) {
                for (const Trig& t : trig) {
                    for (const Ref& r : m.f12) {
                        const double deta = r.eta - t.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hIDmix[c][t.sp][t.pb][TPC_FMD12]->Fill(wrapDPhi(t.phi - r.phi), deta);
                    }
                    for (const Ref& r : m.f3) {
                        const double deta = t.eta - r.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hIDmix[c][t.sp][t.pb][TPC_FMD3]->Fill(wrapDPhi(t.phi - r.phi), deta);
                    }
                }
                for (const Ref& a : f12)
                    for (const Ref& b : m.f3) {
                        const double deta = a.eta - b.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hFFmix[c]->Fill(wrapDPhi(a.phi - b.phi), deta);
                    }
            }

            // FIFO pool, pushed after pairing (an event is never mixed with itself)
            pool[c].push_back({f12, f3});
            if ((int)pool[c].size() > POOL_DEPTH) pool[c].pop_front();
        }
    }

    // ---- write (one file per block: <output>_part<k>.root; a single block keeps the plain name) ----
    std::string outName = outFile;
    if (nChunks > 1) {
        const size_t dot = outFile.rfind(".root");
        const std::string base = (dot == std::string::npos) ? outFile : outFile.substr(0, dot);
        outName = base + "_part" + std::to_string(part) + ".root";
    }
    TFile* fout = new TFile(outName.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create " << outName << std::endl; return 1; }
    for (int c = 0; c < NCLS; ++c) {
        for (int s = 0; s < NSP; ++s) {
            hNtrig[c][s]->Write();
            for (int b = 0; b < NpT_bins; ++b)
                for (int p = 0; p < NIDPAIR; ++p) { hID[c][s][b][p]->Write(); hIDmix[c][s][b][p]->Write(); }
        }
        hFF[c]->Write();
        hFFmix[c]->Write();
        hNassoc[c][0]->Write();
        hNassoc[c][1]->Write();
        hEvents[c]->Write();
        hNchClass[c]->Write();
    }
    hCounter->Write();
    hNchCent->Write();
    hNchV0M->Write();
    fout->Close();

    std::cout << "Done. Wrote " << outName
              << " (block " << (part + 1) << "/" << nChunks
              << ", events=" << ((last - first) + (part == 0 ? nNotWritten : 0))
              << ", registered=" << hCounter->GetBinContent(2)
              << ", HM=" << hEvents[HM]->GetEntries()
              << ", LM=" << hEvents[LM]->GetEntries()
              << ", HM35=" << hEvents[HM35]->GetEntries() << ")" << std::endl;

    // free this block's histograms before the next block
    for (int c = 0; c < NCLS; ++c) {
        for (int s = 0; s < NSP; ++s) {
            delete hNtrig[c][s];
            for (int b = 0; b < NpT_bins; ++b)
                for (int p = 0; p < NIDPAIR; ++p) { delete hID[c][s][b][p]; delete hIDmix[c][s][b][p]; }
        }
        delete hFF[c]; delete hFFmix[c];
        delete hNassoc[c][0]; delete hNassoc[c][1];
        delete hEvents[c]; delete hNchClass[c];
    }
    delete hCounter; delete hNchCent; delete hNchV0M;
    delete fout;
    }   // end of block loop

    fin->Close();
    return 0;
}