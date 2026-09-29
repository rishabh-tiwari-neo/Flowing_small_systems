#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

// =====================================================================
// GENERATION + SAME-EVENT PAIRING, three-subevent method, two
// independent event-activity classifications (Nch, Nmpi), all in
// memory, single process per task.
//
// Eta windows (trigger/associated arms) follow the ALICE identified-
// hadron paper's three-subevent method:
//   TPC   : |eta| < 0.8            (trigger arm)
//   FMD12 : 1.7 < eta < 5.1        (associated arm, forward)
//   FMD3  : -3.1 < eta < -1.7      (associated arm, backward)
//
// Three pairings are built per event:
//   pairing 0: TPC   - FMD12   (deta = eta_FMD12 - eta_TPC,   ~[0.9,5.9])
//   pairing 1: TPC   - FMD3    (deta = eta_TPC   - eta_FMD3,  ~[0.9,3.9])
//   pairing 2: FMD12 - FMD3    (deta = eta_FMD12 - eta_FMD3,  ~[3.4,8.2])
// deta is defined so it is always positive for the natural pairing
// direction; axis ranges below are booked with a small margin.
//
// Event classification follows Ortiz, Sahu & Bencedi (arXiv:2512.09195):
// each event is classified TWICE, independently, by
//   - Nch  : charged particles, |eta| < 0.8, 0.2 < pT < 3.0 GeV/c
//   - Nmpi : MC-truth number of multiparton interactions
// into 6 classes each (I-VI, low -> high activity) -> 12 class-slots
// total. The same event fills BOTH its Nch-class histograms and its
// Nmpi-class histograms (these are two separate, parallel binnings of
// the same event sample, not a 6x6 cross product).
//
// *** CLASS EDGES BELOW ARE PLACEHOLDERS ***
// Ortiz et al. define classes so each one holds the same FRACTION of
// events (percentile binning), which requires knowing the inclusive
// Nch/Nmpi distributions beforehand. This generator fills hNchIncl and
// hNmpiIncl for exactly that purpose: run a short calibration pass,
// derive equal-fraction edges from those two histograms, then update
// kNchEdges/kNmpiEdges below before a full production run.
//
// No mixed-event background is built here: at generator (truth) level
// there is no detector acceptance to correct for. Non-flow subtraction
// (template fit against the lowest-activity class) is done offline in
// post-processing.
//
// No particle identification yet -- all final charged particles are
// treated identically (extension point for later: tag pid per particle
// and split histograms by species, as in the ALICE paper's pi/K/p/K0S/
// Lambda breakdown).
//
// argv[1] = n_events
// argv[2] = seed
// argv[3] = output_file
// =====================================================================

// ---- acceptance windows -------------------------------------------
static const double TPC_ETA_MAX    = 0.8;
static const double FMD12_ETA_MIN  = 1.7,  FMD12_ETA_MAX = 5.1;
static const double FMD3_ETA_MIN   = -3.1, FMD3_ETA_MAX  = -1.7;

static const double NCH_ETA_MAX = 0.8;
static const double NCH_PT_MIN  = 0.2, NCH_PT_MAX = 3.0;

// ---- event classification (PLACEHOLDER edges, see note above) -----
static const int N_CLASSES = 6;
// upper-inclusive edges; last edge is a large sentinel (no upper bound)
static const double kNchEdges[N_CLASSES]  = {2, 5, 9, 14, 22, 1.0e9};
static const double kNmpiEdges[N_CLASSES] = {1, 2, 3, 5, 8, 1.0e9};

enum Estimator { kEstNch = 0, kEstNmpi = 1, kNEstimators = 2 };
static const char* kEstName[kNEstimators] = {"Nch", "Nmpi"};

enum Pairing { kPairTpcFmd12 = 0, kPairTpcFmd3 = 1, kPairFmd12Fmd3 = 2, kNPairings = 3 };
static const char* kPairName[kNPairings] = {"TPC_FMD12", "TPC_FMD3", "FMD12_FMD3"};

// deta axis ranges per pairing (lo, hi), with margin around the
// geometric min/max quoted in the header comment
static const double kDEtaLo[kNPairings] = {0.5, 0.5, 3.0};
static const double kDEtaHi[kNPairings] = {6.5, 4.5, 9.0};
static const int    kDEtaBins[kNPairings] = {24, 16, 24};

static const int    kNDPhiBins = 20;
static const double kDPhiLo = -M_PI / 2.0;
static const double kDPhiHi =  3.0 * M_PI / 2.0;

inline int ClassIndex(double val, const double* edges) {
    for (int c = 0; c < N_CLASSES; ++c) if (val <= edges[c]) return c;
    return N_CLASSES - 1;
}

inline int GlobalClass(int estimator, int cls) { return estimator * N_CLASSES + cls; }
static const int N_GLOBAL_CLASSES = kNEstimators * N_CLASSES; // 12

inline double WrapDPhi(double dphi) {
    while (dphi < kDPhiLo)  dphi += 2.0 * M_PI;
    while (dphi >= kDPhiHi) dphi -= 2.0 * M_PI;
    return dphi;
}

// --- Manually growable flat array: (eta, phi), doubling capacity on
//     overflow -- same style/logic as the earlier two-arm generator. ---
struct DynArr {
    double *eta = nullptr;
    double *phi = nullptr;
    long capacity = 0;
    long size     = 0;

    void reserve(long newCap) {
        double *neweta = new double[newCap];
        double *newphi = new double[newCap];
        if (eta) {
            std::copy(eta, eta + size, neweta);
            std::copy(phi, phi + size, newphi);
            delete[] eta; delete[] phi;
        }
        eta = neweta; phi = newphi;
        capacity = newCap;
    }

    inline void push(double e, double p) {
        if (size >= capacity) reserve(capacity == 0 ? 64 : capacity * 2);
        eta[size] = e; phi[size] = p;
        size++;
    }

    ~DynArr() { delete[] eta; delete[] phi; }
};

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <n_events> <seed> <output_file>" << std::endl;
        return 1;
    }

    long n_events        = std::atol(argv[1]);
    int  seed            = std::atoi(argv[2]);
    std::string outFile  = argv[3];
    std::string progFile = "logs/progress_task" + std::to_string(seed) + ".txt";

    Pythia8::Pythia pythia;
    pythia.readString("Print:quiet = on");

    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13000.0");

    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");

    // rope hadronization + shoving (Table I, Ortiz/Sahu/Bencedi; g=10
    // primary/conservative choice)
    pythia.readString("Ropewalk:RopeHadronization = on");
    pythia.readString("Ropewalk:doShoving = on");
    pythia.readString("Ropewalk:gAmplitude = 10");
    pythia.readString("Ropewalk:r0 = 0.41");
    pythia.readString("Ropewalk:rCutOff = 10");
    pythia.readString("Ropewalk:deltay = 0.10");

    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("PartonVertex:ProtonRadius = 0.7");

    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    pythia.init();

    // ---- book histograms -------------------------------------------
    TH2D* hCorr[kNPairings][N_GLOBAL_CLASSES];
    for (int p = 0; p < kNPairings; ++p) {
        for (int g = 0; g < N_GLOBAL_CLASSES; ++g) {
            int est = g / N_CLASSES;
            int cls = g % N_CLASSES;
            TString name  = TString::Format("hCorr_%s_%s_c%d",
                                             kPairName[p], kEstName[est], cls);
            TString title = TString::Format("%s, %s class %d;#Delta#eta;#Delta#phi",
                                             kPairName[p], kEstName[est], cls);
            hCorr[p][g] = new TH2D(name, title,
                                    kDEtaBins[p], kDEtaLo[p], kDEtaHi[p],
                                    kNDPhiBins, kDPhiLo, kDPhiHi);
            hCorr[p][g]->Sumw2();
        }
    }

    // per-pairing trigger-arm counters, indexed by global class, used
    // for normalization at the post-processing stage
    TH1D* hNtrig[kNPairings];
    for (int p = 0; p < kNPairings; ++p) {
        TString name = TString::Format("hNtrig_%s", kPairName[p]);
        hNtrig[p] = new TH1D(name, "trigger-arm particles per class;global class;N_{trig}",
                              N_GLOBAL_CLASSES, -0.5, N_GLOBAL_CLASSES - 0.5);
        hNtrig[p]->Sumw2();
    }

    // event counters per estimator (6 bins each)
    TH1D* hEvents[kNEstimators];
    for (int e = 0; e < kNEstimators; ++e) {
        TString name = TString::Format("hEvents_%s", kEstName[e]);
        hEvents[e] = new TH1D(name, "events per class;class;N_{ev}",
                               N_CLASSES, -0.5, N_CLASSES - 0.5);
        hEvents[e]->Sumw2();
    }

    // inclusive distributions, for deriving equal-fraction class edges
    TH1D* hNchIncl  = new TH1D("hNchIncl",  "N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0);N_{ch};events",
                                301, -0.5, 300.5);
    TH1D* hNmpiIncl = new TH1D("hNmpiIncl", "N_{mpi};N_{mpi};events", 61, -0.5, 60.5);
    hNchIncl->Sumw2();
    hNmpiIncl->Sumw2();

    // ---- event loop ---------------------------------------------------
    DynArr tpcArr, fmd12Arr, fmd3Arr;
    tpcArr.reserve(64);
    fmd12Arr.reserve(64);
    fmd3Arr.reserve(64);

    long long total_ntrig[kNPairings] = {0, 0, 0};

    for (long li = 0; li < n_events; ++li) {
        tpcArr.size   = 0;
        fmd12Arr.size = 0;
        fmd3Arr.size  = 0;

        if (!pythia.next()) continue;

        int nch = 0;
        int entries = pythia.event.size();
        for (int j = 0; j < entries; j++) {
            if (!pythia.event[j].isFinal()) continue;
            if (pythia.event[j].charge() == 0) continue;

            double eta = pythia.event[j].eta();
            double phi = pythia.event[j].phi();
            double pt  = pythia.event[j].pT();

            if (std::fabs(eta) < NCH_ETA_MAX && pt > NCH_PT_MIN && pt < NCH_PT_MAX)
                ++nch;

            if (std::fabs(eta) < TPC_ETA_MAX)
                tpcArr.push(eta, phi);
            if (eta > FMD12_ETA_MIN && eta < FMD12_ETA_MAX)
                fmd12Arr.push(eta, phi);
            if (eta > FMD3_ETA_MIN && eta < FMD3_ETA_MAX)
                fmd3Arr.push(eta, phi);
        }

        int nmpi = pythia.info.nMPI();
        hNchIncl->Fill(nch);
        hNmpiIncl->Fill(nmpi);

        int clsNch  = ClassIndex(nch,  kNchEdges);
        int clsNmpi = ClassIndex(nmpi, kNmpiEdges);
        int gNch    = GlobalClass(kEstNch,  clsNch);
        int gNmpi   = GlobalClass(kEstNmpi, clsNmpi);

        hEvents[kEstNch]->Fill(clsNch);
        hEvents[kEstNmpi]->Fill(clsNmpi);

        // trigger-arm sizes for this event, per pairing
        long nTrigThisEvent[kNPairings] = {tpcArr.size, tpcArr.size, fmd12Arr.size};
        for (int p = 0; p < kNPairings; ++p) {
            hNtrig[p]->Fill(gNch,  (double)nTrigThisEvent[p]);
            hNtrig[p]->Fill(gNmpi, (double)nTrigThisEvent[p]);
            total_ntrig[p] += nTrigThisEvent[p];
        }

        // pairing 0: TPC - FMD12
        for (long a = 0; a < tpcArr.size; a++) {
            for (long b = 0; b < fmd12Arr.size; b++) {
                double deta = fmd12Arr.eta[b] - tpcArr.eta[a];
                if (deta < kDEtaLo[kPairTpcFmd12] || deta > kDEtaHi[kPairTpcFmd12]) continue;
                double dphi = WrapDPhi(tpcArr.phi[a] - fmd12Arr.phi[b]);
                hCorr[kPairTpcFmd12][gNch]->Fill(deta, dphi);
                hCorr[kPairTpcFmd12][gNmpi]->Fill(deta, dphi);
            }
        }

        // pairing 1: TPC - FMD3
        for (long a = 0; a < tpcArr.size; a++) {
            for (long b = 0; b < fmd3Arr.size; b++) {
                double deta = tpcArr.eta[a] - fmd3Arr.eta[b];
                if (deta < kDEtaLo[kPairTpcFmd3] || deta > kDEtaHi[kPairTpcFmd3]) continue;
                double dphi = WrapDPhi(tpcArr.phi[a] - fmd3Arr.phi[b]);
                hCorr[kPairTpcFmd3][gNch]->Fill(deta, dphi);
                hCorr[kPairTpcFmd3][gNmpi]->Fill(deta, dphi);
            }
        }

        // pairing 2: FMD12 - FMD3
        for (long a = 0; a < fmd12Arr.size; a++) {
            for (long b = 0; b < fmd3Arr.size; b++) {
                double deta = fmd12Arr.eta[a] - fmd3Arr.eta[b];
                if (deta < kDEtaLo[kPairFmd12Fmd3] || deta > kDEtaHi[kPairFmd12Fmd3]) continue;
                double dphi = WrapDPhi(fmd12Arr.phi[a] - fmd3Arr.phi[b]);
                hCorr[kPairFmd12Fmd3][gNch]->Fill(deta, dphi);
                hCorr[kPairFmd12Fmd3][gNmpi]->Fill(deta, dphi);
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
    for (int p = 0; p < kNPairings; ++p) {
        for (int g = 0; g < N_GLOBAL_CLASSES; ++g) hCorr[p][g]->Write();
        hNtrig[p]->Write();
    }
    for (int e = 0; e < kNEstimators; ++e) hEvents[e]->Write();
    hNchIncl->Write();
    hNmpiIncl->Write();
    fout->Close();

    std::cout << "[seed " << seed << "] Done. Wrote " << outFile
              << " (events=" << n_events
              << ", ntrig[TPC-FMD12]=" << total_ntrig[kPairTpcFmd12]
              << ", ntrig[TPC-FMD3]=" << total_ntrig[kPairTpcFmd3]
              << ", ntrig[FMD12-FMD3]=" << total_ntrig[kPairFmd12Fmd3]
              << ")" << std::endl;

    return 0;
}