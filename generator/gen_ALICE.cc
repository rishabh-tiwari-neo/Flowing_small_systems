#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <string>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

// ---- acceptance windows -----
static const double FMD12_ETA_MIN = 1.7,  FMD12_ETA_MAX = 5.1;
static const double FMD3_ETA_MIN  = -3.1, FMD3_ETA_MAX  = -1.7;

static const double NCH_ETA_MAX = 0.8;
static const double NCH_PT_MIN  = 0.2, NCH_PT_MAX = 3.0;

// ---- V0 trigger / V0M estimator ------------------------------------
static const double V0A_ETA_MIN = 2.8,  V0A_ETA_MAX = 5.1;
static const double V0C_ETA_MIN = -3.7, V0C_ETA_MAX = -1.7;

// ---- class definitions ----
// HM   : Top ~0.07% of V0M (as in the paper's trigger), plus nchCent > 25
// LM   : V0 coincidence and nchCent < 20
// HM35 : looser V0M cut + nchCent > 25, chosen so <nchCent> ~ 35 as in the paper
//        (from a 1M-event run: V0M >= 46 with nchCent > 25 gives <nchCent> ~ 35, ~ top 13% of events)

static const int V0M_HM_MIN   = 157;   // HM:   (nchV0A + nchV0C) > 157
static const int V0M_HM35_MIN = 45;    // HM35: (nchV0A + nchV0C) > 45
static const int NCH_HM_MIN   = 25;    // HM and HM35: nchCent > 25
static const int NCH_LM_MAX   = 20;    // LM: nchCent < 20 


// ---- identified-particle acceptance ----
static const double PID_ETA_MAX = 0.8;   // pi, K, p: |eta| < 0.8
static const double V0_Y_MAX    = 0.5;   // K0S, Lambda: |y| < 0.5 and, simultaneously, |eta| < PID_ETA_MAX

// ---- common pT binning for all species ----
static const int    NpT_bins = 20;
static const double PT_EDGES[NpT_bins + 1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0,
                                         2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 10.0};

// Cataloging the pT bin of a particle
inline int ptBin(double pt) {
    if (pt < PT_EDGES[0] || pt >= PT_EDGES[NpT_bins]) return -1;
    return int(std::upper_bound(PT_EDGES, PT_EDGES + NpT_bins + 1, pt) - PT_EDGES) - 1;
}

// Event classifier block (EC)
// ---- event classes  ----
enum EvClass { HM = 0, LM = 1, HM35 = 2, NCLS = 3 };// NCLS = number of classes
static const char* CLS_NAME[NCLS] = {"HM", "LM", "HM35"};

// membership of one event in each class (an event can be in several)
inline void classify(int nchV0A, int nchV0C, int nchCent, bool Ecl[NCLS]) {
    for (int c = 0; c < NCLS; ++c) Ecl[c] = false;
    if (!(nchV0A > 0 && nchV0C > 0)) return;                    // MB trigger
    const int v0m = nchV0A + nchV0C;
    Ecl[HM]   = (v0m > V0M_HM_MIN   && nchCent > NCH_HM_MIN);
    Ecl[HM35] = (v0m > V0M_HM35_MIN && nchCent > NCH_HM_MIN);
    Ecl[LM]   = (nchCent < NCH_LM_MAX);
}
//end EC

// ---- identified species ----
enum Species { PI = 0, K, P, K0S, LAM, NSP };
static const char* SP_NAME[NSP] = {"pi", "K", "p", "K0S", "Lam"};

//function block to indentify particles
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
// TPC_FMD12  : identified trigger (|eta|<0.8 for p,K,pi, |y|<0.5 and |eta|<0.8 for k0s and Lambda) with a charged particle in FMD1,2
// TPC_FMD3   : identified trigger (|eta|<0.8 for p,K,pi, |y|<0.5 and |eta|<0.8 for k0s and Lambda) with a charged particle in FMD3
// FMD12_FMD3 : charged particle in FMD1,2 with one in FMD3 (no PID, no pT bins)
// only the first two are used for the identified-particle histograms (NIDPAIR)
enum Pairing { TPC_FMD12 = 0, TPC_FMD3 = 1, FMD12_FMD3 = 2, NPAIR = 3 };
static const int NIDPAIR = 2;
static const char* PAIR_NAME[NPAIR] = {"TPC_FMD12", "TPC_FMD3", "FMD12_FMD3"};

// |deta| window, as in the paper. Pairs outside it are discarded.
static const double DETA_MIN = 1.1;
static const double DETA_MAX = 6.4;

// delta-eta axis of the 2D histograms: 0.1-wide bins over [DETA_MIN, DETA_MAX)
static const int NDETA = 53;

// mixed-event pool: number of previous events (per class) that the triggers are paired with
static const int POOL_DEPTH = 5;

// delta-phi axis: 20 bins over [-pi/2, 3pi/2)
static const int    NDPHI   = 20;
static const double DPHI_LO = -M_PI / 2.0;
static const double DPHI_HI =  3.0 * M_PI / 2.0;

inline double wrapDPhi(double d) {
    while (d <  DPHI_LO) d += 2.0 * M_PI;
    while (d >= DPHI_HI) d -= 2.0 * M_PI;
    return d;
}

// Trigger and Reference particle properties
struct Trig { double eta, phi; int sp, pb; };
struct Ref  { double eta, phi; };
// reference particles of one past event, kept in the mixed-event pool
struct PooledEvent { std::vector<Ref> f12, f3; };

int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <n_events> <seed> <output_file> <g>" << std::endl;
        return 1;
    }

    long n_events        = std::atol(argv[1]);
    int  seed            = std::atoi(argv[2]);
    std::string outFile  = argv[3];
    double gShove        = std::atof(argv[4]);   // shoving strength, set in run_ALICE.sh
    std::string progFile = "logs/log_ALICE/progress_task" + std::to_string(seed) + ".txt"; // Progress monitoring and receiving event deatils

    Pythia8::Pythia pythia;
    pythia.readString("Print:quiet = on");

    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13000.0");

    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");

    pythia.readString("Ropewalk:RopeHadronization = on");
    pythia.readString("Ropewalk:doShoving = on");
    pythia.readString("Ropewalk:gAmplitude = " + std::to_string(gShove));   // from run_ALICE.sh
    pythia.readString("Ropewalk:r0 = 0.41");
    pythia.readString("Ropewalk:rCutOff = 10");
    pythia.readString("Ropewalk:deltay = 0.10");

    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("PartonVertex:ProtonRadius = 0.7");

    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    pythia.init();

    // generator_info.txt next to the output file: PYTHIA version and shoving strength g actually used

    {
        const size_t slash = outFile.find_last_of('/');
        const std::string dir = (slash == std::string::npos) ? "." : outFile.substr(0, slash);
        std::ofstream info(dir + "/generator_info.txt", std::ios::trunc);
        info << "PYTHIA version: " << pythia.settings.parm("Pythia:versionNumber") << "\n"
             << "Ropewalk:gAmplitude (g): " << pythia.settings.parm("Ropewalk:gAmplitude") << "\n";
    }

    // Histograms and counters
    // identified-trigger correlations, same event (S) and mixed event (B): [class][species][pT bin][TPC_FMD12 / TPC_FMD3]
    TH2D* hID[NCLS][NSP][NpT_bins][NIDPAIR];
    TH2D* hIDmix[NCLS][NSP][NpT_bins][NIDPAIR];
    // FMD12-FMD3 (no identified particle, no pT binning)
    TH2D* hFF[NCLS];
    TH2D* hFFmix[NCLS];
    // normalisation
    TH1D* hNtrig[NCLS][NSP];      // identified triggers vs pT (variable bins = PT_EDGES)
    TH1D* hNassoc[NCLS][2];       // associated particles in FMD12 / FMD3, 1 bin
    TH1D* hEvents[NCLS];          // events per class, 1 bin
    TH1D* hNchClass[NCLS];        // nchCent distribution per class (gives <nchCent> of HM, LM, HM35)
    TH1D* hCounter = new TH1D("hCounter", "event counter;;events", 3, 0.5, 3.5);
    hCounter->GetXaxis()->SetBinLabel(1, "generated");
    hCounter->GetXaxis()->SetBinLabel(2, "registered (V0A&&V0C)");
    hCounter->GetXaxis()->SetBinLabel(3, "HM or LM or HM35");

    // inclusive multiplicity distributions (all registered events), for monitoring the cuts
    TH1D* hNchCent = new TH1D("hNchCent", "N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0);N_{ch};events",
                              301, -0.5, 300.5);
    TH1D* hNchV0M  = new TH1D("hNchV0M", "N_{ch} V0A+V0C;N_{ch};events", 401, -0.5, 400.5);

    //creates histograms for each class, species, pT bin, and pairing
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

    // ---- event loop ----
    std::vector<Trig> trig;
    std::vector<Ref>  f12, f3;
    std::deque<PooledEvent> pool[NCLS];     // previous events of each class (mixed-event pool)

    for (long i = 0; i < n_events; ++i) {

        if (!pythia.next()) continue;
        hCounter->Fill(1);

        trig.clear(); f12.clear(); f3.clear();

        // multiplicities and candidate selection
        int nchCent = 0, nchV0A = 0, nchV0C = 0;
        const int entries = pythia.event.size();
        for (int j = 0; j < entries; ++j) {
            const Pythia8::Particle& p = pythia.event[j];
            if (!p.isFinal()) continue;

            const double eta = p.eta(), phi = p.phi(), pt = p.pT();

            if (p.isCharged()) {
                if (std::fabs(eta) < NCH_ETA_MAX && pt > NCH_PT_MIN && pt < NCH_PT_MAX) ++nchCent;
                if (eta > V0A_ETA_MIN && eta < V0A_ETA_MAX) ++nchV0A;
                if (eta > V0C_ETA_MIN && eta < V0C_ETA_MAX) ++nchV0C;
                // FMD reference particles: all charged, all pT, no PID
                if (eta > FMD12_ETA_MIN && eta < FMD12_ETA_MAX) f12.push_back({eta, phi});
                if (eta > FMD3_ETA_MIN  && eta < FMD3_ETA_MAX)  f3.push_back({eta, phi});
            }

            // identified triggers
            const int sp = speciesOf(std::abs(p.id()));
            if (sp < 0) continue;
            const bool isV0  = (sp == K0S || sp == LAM);
            const bool inAcc = isV0 ? (std::fabs(p.y()) < V0_Y_MAX && std::fabs(eta) < PID_ETA_MAX)
                                    : (std::fabs(eta)   < PID_ETA_MAX);
            const int pb = ptBin(pt);
            if (inAcc && pb >= 0) trig.push_back({eta, phi, sp, pb});
        }

        // MB trigger: V0A && V0C coincidence
        if (nchV0A > 0 && nchV0C > 0) {
            hCounter->Fill(2);
            hNchCent->Fill(nchCent);
            hNchV0M->Fill(nchV0A + nchV0C);
        }

        // Event classification
        bool Ecl[NCLS];
        classify(nchV0A, nchV0C, nchCent, Ecl);
        if (Ecl[HM] || Ecl[LM] || Ecl[HM35]) hCounter->Fill(3);
        for (int c = 0; c < NCLS; ++c) {
            if (!Ecl[c]) continue;

            hEvents[c]->Fill(0.5);
            hNchClass[c]->Fill(nchCent);
            hNassoc[c][0]->Fill(0.5, (double)f12.size());
            hNassoc[c][1]->Fill(0.5, (double)f3.size());


            // pairing loops
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

            // mixed events: current triggers / FMD1,2 particles with the references of the pooled events
            for (const PooledEvent& m : pool[c]) {
                for (const Trig& t : trig) {
                    for (const Ref& r : m.f12) {                 // TPC - FMD1,2
                        const double deta = r.eta - t.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hIDmix[c][t.sp][t.pb][TPC_FMD12]->Fill(wrapDPhi(t.phi - r.phi), deta);
                    }
                    for (const Ref& r : m.f3) {                  // TPC - FMD3
                        const double deta = t.eta - r.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hIDmix[c][t.sp][t.pb][TPC_FMD3]->Fill(wrapDPhi(t.phi - r.phi), deta);
                    }
                }
                for (const Ref& a : f12)                         // FMD1,2 - FMD3
                    for (const Ref& b : m.f3) {
                        const double deta = a.eta - b.eta;
                        if (deta < DETA_MIN || deta > DETA_MAX) continue;
                        hFFmix[c]->Fill(wrapDPhi(a.phi - b.phi), deta);
                    }
            }

            // FIFO pool 
            // pushed after pairing, so mixing starts with the 2nd event of the class and saturates at POOL_DEPTH = n
            pool[c].push_back({f12, f3});
            if ((int)pool[c].size() > POOL_DEPTH) pool[c].pop_front();
        }

        // progress monitoring every 100th event
        if (i % 100 == 0 || i == n_events - 1) {
            std::ofstream pf(progFile);
            pf << (i + 1);
            pf.close();
        }
    }

    // ---- write ----
    TFile* fout = new TFile(outFile.c_str(), "RECREATE");
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

    std::cout << "[seed " << seed << "] Done. Wrote " << outFile
              << " (events=" << n_events
              << ", HM=" << hEvents[HM]->GetEntries()
              << ", LM=" << hEvents[LM]->GetEntries()
              << ", HM35=" << hEvents[HM35]->GetEntries() << ")" << std::endl;

    return 0;
}