// AmptRootIO.cxx
// ROOT output of AMPT v1.26t9b-v2.26t9b. The Fortran code (main.f, linana.f)
// calls the extern "C" functions below, which fill the tree "t" while AMPT runs
// (ZSTD level 5 = compression 505).
//
// Tree "t", one entry per event (events without final-state particles are not
// written):
//   gImp                        impact parameter b (fm)
//   Npart                       participant nucleons (AMPT npart1 + npart2)
//   Ncoll                       binary nucleon-nucleon collisions (HIJING NCOLT)
//   psiRP                       reaction-plane angle (0 unless iphirp = 1)
//   Nab                         number of nucleons, A_projectile + A_target
//   Stat[Nab]                   HIJING status of each nucleon, 0 = no interaction
//   Nx, Ny, Nz [Nab]            nucleon positions (fm), projectile first, shifted
//                               by +-b/2 along psiRP
//   ntrack                      number of stored particles
//   PID, charge [ntrack]        PDG code, charge (e)
//   px, py, pz, energy [ntrack] momentum (GeV/c), energy (GeV)
// Stored particles: every charged particle, and K0S, phi, Lambda, Xi, Omega,
// at any eta and pT.
//
// Environment: AMPT_ROOT_OUT = output file (default ana/ampt.root);
// AMPT_ENABLE_TXT=1 / AMPT_ENABLE_ZPC=1 also write ana/ampt.dat / ana/zpc.dat.

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParticlePDG.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unordered_map>

namespace {

constexpr int kMaxNucl   = 600;   // nucleons per event (HIJING: at most 300 per nucleus)
constexpr int kMaxTracks = 70000; // stored particles per event (Pb-Pb: ~60k)

const char* kBranches =
    "gImp = b (fm); Npart = AMPT npart1 + npart2; Ncoll = HIJING NCOLT (binary NN "
    "collisions); psiRP = HIJING reaction-plane angle; Nab = A_proj + A_targ nucleons, "
    "projectile first; Stat = HIJING NFP(i,5) / NFT(i,5), 0: no interaction (as "
    "ana/npart-xy.dat); Nx, Ny, Nz in fm, shifted by +-b/2 along psiRP; tracks PID, charge, "
    "px, py, pz (GeV/c), energy (GeV) of every charged particle (charge != 0, e and mu "
    "included) and K0S, phi, Lambda, Xi, Omega, at any eta and pT";

TFile* outFile = nullptr;
TTree* tree    = nullptr;

// Branch buffers
float gImp  = 0;
int   Npart = 0;
int   Ncoll = 0;
float psiRP = 0;
int   Nab   = 0;
int   Stat[kMaxNucl];
float Nx[kMaxNucl], Ny[kMaxNucl], Nz[kMaxNucl];
int   ntrack = 0;
int   PID[kMaxTracks], charge[kMaxTracks];
float px[kMaxTracks], py[kMaxTracks], pz[kMaxTracks], energy[kMaxTracks];

int       eventId    = 0;
int       nParticles = 0; // final-state particles of the event, stored or not
long long nWritten   = 0;
long long nEmpty     = 0;
int       seedHijing = 0;
int       seedZpc    = 0;

bool EnvFlag(const char* name) {
    const char* v = std::getenv(name);
    return v && (v[0] == '1' || v[0] == 'y' || v[0] == 'Y' || v[0] == 't' || v[0] == 'T');
}

// Charge in units of e: from TDatabasePDG, plus the AMPT deuteron (code 42) and
// nuclei (PDG code 10LZZZAAAI); unknown codes are taken as neutral
int Charge(int pid) {
    static std::unordered_map<int, int> cache;
    if (auto it = cache.find(pid); it != cache.end())
        return it->second;
    const int apid = std::abs(pid);
    double q = 0;
    if (apid == 42)
        q = pid > 0 ? 1 : -1;
    else if (apid > 1000000000)
        q = (pid > 0 ? 1 : -1) * ((apid / 10000) % 1000);
    else if (TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pid))
        q = p->Charge() / 3.0;
    else
        std::cerr << "[AmptRootIO] PDG code " << pid << " unknown, taken as neutral" << std::endl;
    return cache[pid] = std::lround(q);
}

// Stored: charged particles and K0S, phi, Lambda, Xi, Omega (and antiparticles)
bool IsStored(int pid, int q) {
    const int apid = std::abs(pid);
    return q != 0 || pid == 310 || pid == 333 || apid == 3122 || apid == 3312 ||
           apid == 3322 || apid == 3334;
}

} // namespace

extern "C" {

void root_close_();

int iroot_is_txt_enabled_() { return EnvFlag("AMPT_ENABLE_TXT") ? 1 : 0; }

int iroot_is_zpc_enabled_() { return EnvFlag("AMPT_ENABLE_ZPC") ? 1 : 0; }

void root_set_seeds_(int* nseed, int* iseedp) {
    seedHijing = *nseed;
    seedZpc    = *iseedp;
}

void root_init_() {
    const char* outName = std::getenv("AMPT_ROOT_OUT");
    if (!outName || !*outName)
        outName = "ana/ampt.root";

    outFile = new TFile(outName, "RECREATE", "", 505);
    if (outFile->IsZombie()) {
        std::cerr << "[AmptRootIO] ERROR: cannot open the output file " << outName << std::endl;
        std::exit(1);
    }

    tree = new TTree("t", "AMPT events");
    tree->Branch("gImp", &gImp, "gImp/F");
    tree->Branch("Npart", &Npart, "Npart/I");
    tree->Branch("Ncoll", &Ncoll, "Ncoll/I");
    tree->Branch("psiRP", &psiRP, "psiRP/F");
    tree->Branch("Nab", &Nab, "Nab/I");
    tree->Branch("Stat", Stat, "Stat[Nab]/I");
    tree->Branch("Nx", Nx, "Nx[Nab]/F");
    tree->Branch("Ny", Ny, "Ny[Nab]/F");
    tree->Branch("Nz", Nz, "Nz[Nab]/F");
    tree->Branch("ntrack", &ntrack, "ntrack/I");
    tree->Branch("PID", PID, "PID[ntrack]/I");
    tree->Branch("charge", charge, "charge[ntrack]/I");
    tree->Branch("px", px, "px[ntrack]/F");
    tree->Branch("py", py, "py[ntrack]/F");
    tree->Branch("pz", pz, "pz[ntrack]/F");
    tree->Branch("energy", energy, "energy[ntrack]/F");

    // A STOP inside AMPT still closes the file with the events written so far
    std::atexit(root_close_);
}

// Called after HIJING (main.f pass_nucleons_to_root): projectile nucleons first
void root_set_nucleons_(int* na, int* nb, int* stat, float* x, float* y, float* z) {
    Nab = *na + *nb;
    if (Nab > kMaxNucl) {
        std::cerr << "[AmptRootIO] ERROR: " << Nab << " nucleons, more than " << kMaxNucl << std::endl;
        std::exit(1);
    }
    for (int i = 0; i < Nab; i++) {
        Stat[i] = stat[i];
        Nx[i]   = x[i];
        Ny[i]   = y[i];
        Nz[i]   = z[i];
    }
}

// Spectator nucleon counts are passed by main.f but not stored
void root_set_header_(int* ncoll, int* /*nsp1n*/, int* /*nsp1p*/, int* /*nsp2n*/, int* /*nsp2p*/) {
    Ncoll = *ncoll;
}

void root_begin_event_(int* iaevt, int* /*iarun*/, int* /*natt*/, float* bimp, int* npart1,
                       int* npart2, int* /*nelp*/, int* /*ninp*/, int* /*nelt*/, int* /*nint*/,
                       float* phirp) {
    eventId    = *iaevt;
    gImp       = *bimp;
    Npart      = *npart1 + *npart2;
    psiRP      = *phirp;
    ntrack     = 0;
    nParticles = 0;
}

// One final-state particle; its freeze-out position x, y, z, t is not stored
void root_add_particle_(int* ityp, float* ppx, float* ppy, float* ppz, float* m, float* /*x*/,
                        float* /*y*/, float* /*z*/, float* /*t*/) {
    nParticles++;
    const int q = Charge(*ityp);
    if (!IsStored(*ityp, q))
        return;
    if (ntrack >= kMaxTracks) {
        std::cerr << "[AmptRootIO] ERROR: event " << eventId << " has more than " << kMaxTracks
                  << " particles to store" << std::endl;
        std::exit(1);
    }
    const double dpx = *ppx, dpy = *ppy, dpz = *ppz, dm = *m;
    PID[ntrack]    = *ityp;
    charge[ntrack] = q;
    px[ntrack]     = *ppx;
    py[ntrack]     = *ppy;
    pz[ntrack]     = *ppz;
    energy[ntrack] = std::sqrt(dpx * dpx + dpy * dpy + dpz * dpz + dm * dm);
    ntrack++;
}

void root_end_event_() {
    if (nParticles == 0) {
        nEmpty++;
        return;
    }
    tree->Fill();
    nWritten++;
}

void root_close_() {
    if (!outFile)
        return;
    outFile->cd();
    TNamed("source", Form("AMPT v1.26t9b-v2.26t9b (AmptRootIO.cxx), HIJING seed %d, ZPC seed %d",
                          seedHijing, seedZpc))
        .Write();
    TNamed("branches", kBranches).Write();
    TNamed("tEvents", Form("written %lld, not written: no particles %lld", nWritten, nEmpty)).Write();
    outFile->Write();
    outFile->Close();
    delete outFile;
    outFile = nullptr;
    // summary line read by the production scripts
    std::cout << "[AmptRootIO] tree t events: " << nWritten << " (no particles: " << nEmpty << ")"
              << std::endl;
}

} // extern "C"
