#include <map>
#include <cmath>
void see_AMPT(const char* fn = "mb100k.root") {
  TFile f(fn); TTree* t = (TTree*)f.Get("t");
  const int N = 70000;
  static int PID[N], charge[N]; static float px[N], py[N], pz[N], en[N];
  int ntrack;
  t->SetBranchStatus("*", 0);
  for (auto b : {"ntrack","PID","charge","px","py","pz","energy"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("ntrack", &ntrack); t->SetBranchAddress("PID", PID);
  t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("px", px); t->SetBranchAddress("py", py);
  t->SetBranchAddress("pz", pz); t->SetBranchAddress("energy", en);

  std::map<int,long> cnt;                 // all stored PIDs
  long nEv = t->GetEntries(), nReg = 0;
  double sCent=0, sA=0, sC=0, sF12=0, sF3=0;           // sums over registered events
  double sp[5] = {0,0,0,0,0};                          // pi K p K0S Lam in acceptance, pT 0.5-10
  const char* spn[5] = {"pi","K","p","K0S","Lam"};
  long nHM=0, nHM35=0, nLM=0;

  for (long i = 0; i < nEv; ++i) {
    t->GetEntry(i);
    int cent=0, a=0, c=0, f12=0, f3=0; int s[5] = {0,0,0,0,0};
    for (int j = 0; j < ntrack && j < N; ++j) {
      cnt[PID[j]]++;
      double pt = std::hypot(px[j], py[j]);
      double p  = std::sqrt(pt*pt + pz[j]*pz[j]);
      if (p == std::fabs(pz[j])) continue;
      double eta = 0.5*std::log((p + pz[j])/(p - pz[j]));
      if (charge[j] != 0) {
        if (std::fabs(eta) < 0.8 && pt > 0.2 && pt < 3.0) ++cent;
        if (eta > 2.8  && eta < 5.1)  ++a;
        if (eta > -3.7 && eta < -1.7) ++c;
        if (eta > 1.7  && eta < 5.1)  ++f12;
        if (eta > -3.1 && eta < -1.7) ++f3;
      }
      int ap = std::abs(PID[j]), k = -1;
      if (ap==211) k=0; else if (ap==321) k=1; else if (ap==2212) k=2;
      else if (ap==310) k=3; else if (ap==3122) k=4;
      if (k < 0 || pt < 0.5 || pt >= 10.0) continue;
      bool acc = std::fabs(eta) < 0.8;
      if (k >= 3) {
        double y = 0.5*std::log((en[j] + pz[j])/(en[j] - pz[j]));
        acc = acc && std::fabs(y) < 0.5;
      }
      if (acc) ++s[k];
    }
    if (!(a > 0 && c > 0)) continue;       // MB trigger
    ++nReg;
    sCent += cent; sA += a; sC += c; sF12 += f12; sF3 += f3;
    for (int k = 0; k < 5; ++k) sp[k] += s[k];
    int v0m = a + c;
    if (v0m > 157 && cent > 25) ++nHM;
    if (v0m > 45  && cent > 25) ++nHM35;
    if (cent < 20) ++nLM;
  }

  printf("\n== stored PIDs (all events, all eta, pT) ==\n%8s %12s %10s\n", "PID", "count", "per event");
  for (auto& x : cnt) printf("%8d %12ld %10.4f\n", x.first, x.second, (double)x.second/nEv);

  printf("\n== events: total %ld, registered (V0A&&V0C) %ld (%.1f%%) ==\n", nEv, nReg, 100.*nReg/nEv);
  printf("== mean per registered event ==\n");
  printf("nchCent %.2f  V0A %.2f  V0C %.2f  V0M %.2f  FMD12 %.2f  FMD3 %.2f\n",
         sCent/nReg, sA/nReg, sC/nReg, (sA+sC)/nReg, sF12/nReg, sF3/nReg);
  printf("== triggers per registered event (acceptance, 0.5<pT<10) ==\n");
  for (int k = 0; k < 5; ++k) printf("%-4s %.4f\n", spn[k], sp[k]/nReg);
  printf("== class sizes with the PYTHIA cuts (V0M>157, V0M>45, nchCent>25 / <20) ==\n");
  printf("HM %ld  HM35 %ld (%.1f%% of registered)  LM %ld\n", nHM, nHM35, 100.*nHM35/nReg, nLM);
}