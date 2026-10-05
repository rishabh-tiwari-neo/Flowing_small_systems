#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// ---- acceptance windows -------------------------------------------
static const double TPC_ETA_MAX    = 0.8;
static const double FMD12_ETA_MIN  = 1.7,  FMD12_ETA_MAX = 5.1;
static const double FMD3_ETA_MIN   = -3.1, FMD3_ETA_MAX  = -1.7;

static const double V0A_ETA_MIN    = 2.8,  V0A_ETA_MAX   = 5.1;
static const double V0C_ETA_MIN    = -3.7, V0C_ETA_MAX   = -1.7;

static const double V0_Y_MAX = 0.5;   // K0S, Lambda/Lambdabar accepted by rapidity |y| < 0.5

static const double NCH_ETA_MAX = 0.8;
static const double NCH_PT_MIN  = 0.2, NCH_PT_MAX = 3.0;

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

    // ---- output tree ---------------------------------------------------
    TFile* fout = new TFile(outFile.c_str(), "RECREATE");
    TTree* tree = new TTree("events", "pp 13 TeV, unclassified");

    // event-level
    long  evtNo    = 0;
    int   nMPI     = 0;
    int   nchCent  = 0;  // |eta|<0.8, 0.2<pT<3.0, charged final state
    int   nchV0A   = 0;  // charged final state in V0A range
    int   nchV0C   = 0;  // charged final state in V0C range
    int   nchFMD12 = 0;  // charged final state in FMD1,2 range
    int   nchFMD3  = 0;  // charged final state in FMD3 range

    // particle-level (final-state charged + K0S + Lambda/Lambdabar
    // inside a detector window; K0S/Lambda: |y| < 0.5)
    std::vector<float> pt, eta, phi, rap;
    std::vector<int>   pid;
    std::vector<short> charge;

    tree->Branch("evtNo",    &evtNo);
    tree->Branch("nMPI",     &nMPI);
    tree->Branch("nchCent",  &nchCent);
    tree->Branch("nchV0A",   &nchV0A);
    tree->Branch("nchV0C",   &nchV0C);
    tree->Branch("nchFMD12", &nchFMD12);
    tree->Branch("nchFMD3",  &nchFMD3);

    tree->Branch("pt",     &pt);
    tree->Branch("eta",    &eta);
    tree->Branch("phi",    &phi);
    tree->Branch("y",      &rap);
    tree->Branch("pid",    &pid);
    tree->Branch("charge", &charge);

    // ---- histograms ----------------------------------------------------
    // nchCent in bins of 10: [0,10), [10,20), ... [90,100)
    TH1D* hNch10 = new TH1D("hNch10",
        "N_{ch} (|#eta|<0.8, 0.2<p_{T}<3.0);N_{ch};events", 10, 0, 100);
    hNch10->Sumw2();

    // forward-detector charged multiplicities, same bins of 10: [0,10) ... [90,100)
    TH1D* hNchV0A = new TH1D("hNchV0A",
        "N_{ch} in V0A (2.8<#eta<5.1);N_{ch};events", 10, 0, 100);
    TH1D* hNchV0C = new TH1D("hNchV0C",
        "N_{ch} in V0C (-3.7<#eta<-1.7);N_{ch};events", 10, 0, 100);
    TH1D* hNchFMD12 = new TH1D("hNchFMD12",
        "N_{ch} in FMD1,2 (1.7<#eta<5.1);N_{ch};events", 10, 0, 100);
    TH1D* hNchFMD3 = new TH1D("hNchFMD3",
        "N_{ch} in FMD3 (-3.1<#eta<-1.7);N_{ch};events", 10, 0, 100);
    hNchV0A->Sumw2(); hNchV0C->Sumw2();
    hNchFMD12->Sumw2(); hNchFMD3->Sumw2();

    // nMPI variable bins: [0,3), [3,7), [7,11), [11,15), [15,19), [19,23), [23,34)
    const double mpiEdges[] = {0, 3, 7, 11, 15, 19, 23, 34};
    TH1D* hNmpi = new TH1D("hNmpi",
        "N_{MPI};N_{MPI};events", 7, mpiEdges);
    hNmpi->Sumw2();

    // ---- event loop ---------------------------------------------------
    for (long li = 0; li < n_events; ++li) {

        if (!pythia.next()) continue;

        pt.clear(); eta.clear(); phi.clear(); rap.clear();
        pid.clear(); charge.clear();
        nchCent = nchV0A = nchV0C = nchFMD12 = nchFMD3 = 0;

        const int entries = pythia.event.size();
        for (int j = 0; j < entries; j++) {
            const Pythia8::Particle& p = pythia.event[j];
            if (!p.isFinal()) continue;

            const int  id  = p.id();
            const int  q   = p.charge();   // in units of e
            const bool isCharged = (p.isCharged());
            const bool isK0S     = (id == 310);
            const bool isLambda  = (std::abs(id) == 3122);

            if (!isCharged && !isK0S && !isLambda) continue;

            const double e  = p.eta();
            const double pp = p.pT();

            if (isCharged) {
                if (std::fabs(e) < NCH_ETA_MAX && pp > NCH_PT_MIN && pp < NCH_PT_MAX)
                    ++nchCent;
                if (e > V0A_ETA_MIN   && e < V0A_ETA_MAX)   ++nchV0A;
                if (e > V0C_ETA_MIN   && e < V0C_ETA_MAX)   ++nchV0C;
                if (e > FMD12_ETA_MIN && e < FMD12_ETA_MAX) ++nchFMD12;
                if (e > FMD3_ETA_MIN  && e < FMD3_ETA_MAX)  ++nchFMD3;
            }

            // K0S and Lambda/Lambdabar: stored by rapidity, |y| < V0_Y_MAX
            // charged particles: stored if inside at least one detector window (eta)
            bool store;
            if (isK0S || isLambda) {
                store = std::fabs(p.y()) < V0_Y_MAX;
            } else {
                store = (e > V0C_ETA_MIN   && e < V0C_ETA_MAX)   ||
                        (std::fabs(e) < TPC_ETA_MAX)             ||
                        (e > FMD12_ETA_MIN && e < FMD12_ETA_MAX) ||
                        (e > FMD3_ETA_MIN  && e < FMD3_ETA_MAX)  ||
                        (e > V0A_ETA_MIN   && e < V0A_ETA_MAX);
            }
            if (!store) continue;

            pt.push_back(pp);
            eta.push_back(e);
            phi.push_back(p.phi());
            rap.push_back(p.y());
            pid.push_back(id);
            charge.push_back((short)q);
        }

        nMPI  = pythia.info.nMPI();
        evtNo = li;

        hNch10->Fill(nchCent);
        hNchV0A->Fill(nchV0A);
        hNchV0C->Fill(nchV0C);
        hNchFMD12->Fill(nchFMD12);
        hNchFMD3->Fill(nchFMD3);
        hNmpi->Fill(nMPI);

        tree->Fill();

        if (li % 100 == 0 || li == n_events - 1) {
            std::ofstream pf(progFile);
            pf << (li + 1);
            pf.close();
        }
    }

    // ---- write --------------------------------------------------------
    fout->cd();
    tree->Write();
    hNch10->Write();
    hNchV0A->Write();
    hNchV0C->Write();
    hNchFMD12->Write();
    hNchFMD3->Write();
    hNmpi->Write();
    fout->Close();

    std::cout << "[seed " << seed << "] Done. Wrote " << outFile
              << " (events=" << n_events << ")" << std::endl;

    return 0;
}