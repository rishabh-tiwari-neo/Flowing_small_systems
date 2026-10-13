#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

// ---- acceptance windows -------------------------------------------
// symmetric windows are used as fabs(eta) < X, so one constant covers both edges
static const double FMD12_ETA_MIN = 1.7,  FMD12_ETA_MAX = 5.1;
static const double FMD3_ETA_MIN  = -3.1, FMD3_ETA_MAX  = -1.7;

static const double NCH_ETA_MAX = 0.8;
static const double NCH_PT_MIN  = 0.2, NCH_PT_MAX = 3.0;

// ---- V0 trigger / V0M estimator ------------------------------------
static const double V0A_ETA_MIN = 2.8,  V0A_ETA_MAX = 5.1;
static const double V0C_ETA_MIN = -3.7, V0C_ETA_MAX = -1.7;
static const int    V0M_HM_MIN  = 157;   // HM: (nchV0A + nchV0C) > 157 (~top 0.07% in the proxy)

// ---- class definitions (as in the paper) ---------------------------
static const int NCH_HM_MIN = 25;        // HM: V0M cut && nchCent > 25
static const int NCH_LM_MAX = 20;        // LM: V0 coincidence && nchCent < 20

// ---- identified-particle acceptance --------------------------------
static const double PID_ETA_MAX = 0.8;   // pi, K, p: |eta| < 0.8
static const double V0_Y_MAX    = 0.5;   // K0S, Lambda: |y| < 0.5

// ---- common pT binning for all species -----------------------------
static const int    NPT = 11;
static const double PT_EDGES[NPT + 1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5,
                                         1.75, 2.0, 3.0, 4.0, 5.0, 10.0};

inline int ptBin(double pt) {
    if (pt < PT_EDGES[0] || pt >= PT_EDGES[NPT]) return -1;
    return int(std::upper_bound(PT_EDGES, PT_EDGES + NPT + 1, pt) - PT_EDGES) - 1;
}

// ---- event classes (integer tags, used as array indices) -------------
// NCLS = number of classes; NONE = neither HM nor LM (event is skipped)
enum EvClass { HM = 0, LM = 1, NCLS = 2, NONE = -1 };
static const char* CLS_NAME[NCLS] = {"HM", "LM"};

inline int classify(int nchV0A, int nchV0C, int nchCent) {
    if (!(nchV0A > 0 && nchV0C > 0)) return NONE;               // MB trigger
    if ((nchV0A + nchV0C) > V0M_HM_MIN && nchCent > NCH_HM_MIN) return HM;
    if (nchCent < NCH_LM_MAX) return LM;
    return NONE;
}

// ---- species (identified triggers) ---------------------------------
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

// ---- pairings (integer tags, used as array indices) --------------------
// TPC_FMD12  : identified trigger (|eta|<0.8) with a charged particle in FMD1,2
// TPC_FMD3   : identified trigger (|eta|<0.8) with a charged particle in FMD3
// FMD12_FMD3 : charged particle in FMD1,2 with one in FMD3 (no PID, no pT bins)
// only the first two are used for the identified-particle histograms (NIDPAIR)
enum Pairing { TPC_FMD12 = 0, TPC_FMD3 = 1, FMD12_FMD3 = 2, NPAIR = 3 };
static const int NIDPAIR = 2;
static const char* PAIR_NAME[NPAIR] = {"TPC_FMD12", "TPC_FMD3", "FMD12_FMD3"};

// |deta| window, as in the paper (pp). Pairs outside it are discarded.
static const double DETA_MIN = 1.1;
static const double DETA_MAX = 6.4;

// delta-phi axis: 20 bins over [-pi/2, 3pi/2)
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

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <n_events> <seed> <output_file>" << std::endl;
        return 1;
    }

    long n_events        = std::atol(argv[1]);
    int  seed            = std::atoi(argv[2]);
    std::string outFile  = argv[3];
    std::string progFile = "logs/log_ALICE/progress_task" + std::to_string(seed) + ".txt";

    Pythia8::Pythia pythia;
    pythia.readString("Print:quiet = on");

    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13000.0");

    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");

    pythia.readString("Ropewalk:RopeHadronization = on");
    pythia.readString("Ropewalk:doShoving = on");
    pythia.readString("Ropewalk:gAmplitude = 10");   // change per run: 3 / 10 / 40
    pythia.readString("Ropewalk:r0 = 0.41");
    pythia.readString("Ropewalk:rCutOff = 10");
    pythia.readString("Ropewalk:deltay = 0.10");

    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("PartonVertex:ProtonRadius = 0.7");

    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    pythia.init();

    // ---- book histograms (all 1D in delta-phi) ----------------------
    // identified-trigger correlations: [class][species][pT bin][TPC_FMD12 / TPC_FMD3]
    TH1D* hID[NCLS][NSP][NPT][NIDPAIR];
    // FMD12-FMD3 (no identified particle, no pT binning)
    TH1D* hFF[NCLS];
    // normalisation
    TH1D* hNtrig[NCLS][NSP];      // identified triggers vs pT (variable bins = PT_EDGES)
    TH1D* hNtrigFF[NCLS];         // FMD12 particles (trigger arm of FMD12-FMD3), 1 bin
    TH1D* hNassoc[NCLS][2];       // associated particles in FMD12 / FMD3, 1 bin
    TH1D* hEvents[NCLS];          // events per class, 1 bin
    TH1D* hCounter = new TH1D("hCounter", "event counter;;events", 3, 0.5, 3.5);
    hCounter->GetXaxis()->SetBinLabel(1, "generated");
    hCounter->GetXaxis()->SetBinLabel(2, "registered (V0A&&V0C)");
    hCounter->GetXaxis()->SetBinLabel(3, "HM or LM");

    // inclusive multiplicity distributions (all registered events), for monitoring the HM/LM cuts
    TH1D* hNchCent = new TH1D("hNchCent", "N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0);N_{ch};events",
                              301, -0.5, 300.5);
    TH1D* hNchV0M  = new TH1D("hNchV0M", "N_{ch} V0A+V0C;N_{ch};events", 401, -0.5, 400.5);
    hNchCent->Sumw2();
    hNchV0M->Sumw2();

    for (int c = 0; c < NCLS; ++c) {
        for (int s = 0; s < NSP; ++s) {
            hNtrig[c][s] = new TH1D(Form("hNtrig_%s_%s", CLS_NAME[c], SP_NAME[s]),
                                    ";p_{T} (GeV/c);N_{trig}", NPT, PT_EDGES);
            hNtrig[c][s]->Sumw2();
            for (int b = 0; b < NPT; ++b)
                for (int p = 0; p < NIDPAIR; ++p) {
                    hID[c][s][b][p] = new TH1D(
                        Form("hCorr_%s_%s_%s_pt%d", CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], b),
                        Form("%s %s %s, p_{T} %.2f-%.2f;#Delta#phi;pairs",
                             CLS_NAME[c], PAIR_NAME[p], SP_NAME[s], PT_EDGES[b], PT_EDGES[b + 1]),
                        NDPHI, DPHI_LO, DPHI_HI);
                    hID[c][s][b][p]->Sumw2();
                }
        }
        hFF[c] = new TH1D(Form("hCorr_%s_FMD12_FMD3", CLS_NAME[c]),
                          Form("%s FMD12_FMD3;#Delta#phi;pairs", CLS_NAME[c]),
                          NDPHI, DPHI_LO, DPHI_HI);
        hFF[c]->Sumw2();
        hNtrigFF[c]   = new TH1D(Form("hNtrig_%s_FMD12",  CLS_NAME[c]), "N_{trig} FMD12",  1, 0, 1);
        hNassoc[c][0] = new TH1D(Form("hNassoc_%s_FMD12", CLS_NAME[c]), "N_{assoc} FMD12", 1, 0, 1);
        hNassoc[c][1] = new TH1D(Form("hNassoc_%s_FMD3",  CLS_NAME[c]), "N_{assoc} FMD3",  1, 0, 1);
        hEvents[c]    = new TH1D(Form("hEvents_%s", CLS_NAME[c]), "events", 1, 0, 1);
        hNtrigFF[c]->Sumw2();
        hNassoc[c][0]->Sumw2();
        hNassoc[c][1]->Sumw2();
        hEvents[c]->Sumw2();
    }

    // ---- event loop ---------------------------------------------------
    std::vector<Trig> trig;
    std::vector<Ref>  f12, f3;

    for (long li = 0; li < n_events; ++li) {

        if (!pythia.next()) continue;
        hCounter->Fill(1);

        trig.clear(); f12.clear(); f3.clear();

        // pass 1: multiplicities + candidate collection
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
            const bool inAcc = isV0 ? (std::fabs(p.y()) < V0_Y_MAX)
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

        const int c = classify(nchV0A, nchV0C, nchCent);
        if (c != NONE) {
            hCounter->Fill(3);
            hEvents[c]->Fill(0.5);
            hNtrigFF[c]->Fill(0.5, (double)f12.size());
            hNassoc[c][0]->Fill(0.5, (double)f12.size());
            hNassoc[c][1]->Fill(0.5, (double)f3.size());

            // pass 2: pairs
            for (const Trig& t : trig) {
                hNtrig[c][t.sp]->Fill(0.5 * (PT_EDGES[t.pb] + PT_EDGES[t.pb + 1]));

                for (const Ref& r : f12) {                       // TPC - FMD1,2
                    const double deta = r.eta - t.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hID[c][t.sp][t.pb][TPC_FMD12]->Fill(wrapDPhi(t.phi - r.phi));
                }
                for (const Ref& r : f3) {                        // TPC - FMD3
                    const double deta = t.eta - r.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hID[c][t.sp][t.pb][TPC_FMD3]->Fill(wrapDPhi(t.phi - r.phi));
                }
            }
            for (const Ref& a : f12)                             // FMD1,2 - FMD3
                for (const Ref& b : f3) {
                    const double deta = a.eta - b.eta;
                    if (deta < DETA_MIN || deta > DETA_MAX) continue;
                    hFF[c]->Fill(wrapDPhi(a.phi - b.phi));
                }
        }

        if (li % 100 == 0 || li == n_events - 1) {
            std::ofstream pf(progFile);
            pf << (li + 1);
            pf.close();
        }
    }

    // ---- write --------------------------------------------------------
    TFile* fout = new TFile(outFile.c_str(), "RECREATE");
    for (int c = 0; c < NCLS; ++c) {
        for (int s = 0; s < NSP; ++s) {
            hNtrig[c][s]->Write();
            for (int b = 0; b < NPT; ++b)
                for (int p = 0; p < NIDPAIR; ++p) hID[c][s][b][p]->Write();
        }
        hFF[c]->Write();
        hNtrigFF[c]->Write();
        hNassoc[c][0]->Write();
        hNassoc[c][1]->Write();
        hEvents[c]->Write();
    }
    hCounter->Write();
    hNchCent->Write();
    hNchV0M->Write();
    fout->Close();

    std::cout << "[seed " << seed << "] Done. Wrote " << outFile
              << " (events=" << n_events
              << ", HM=" << hEvents[HM]->GetEntries()
              << ", LM=" << hEvents[LM]->GetEntries() << ")" << std::endl;

    return 0;
}