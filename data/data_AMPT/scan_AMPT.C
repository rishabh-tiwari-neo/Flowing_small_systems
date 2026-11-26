#include <vector>
#include <algorithm>
#include <cmath>
void scan_AMPT(const char* fn = "mb100k.root") {
  TFile f(fn); TTree* t = (TTree*)f.Get("t");
  const int N = 70000;
  static int charge[N]; static float px[N], py[N], pz[N];
  int ntrack;
  t->SetBranchStatus("*", 0);
  for (auto b : {"ntrack","charge","px","py","pz"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("ntrack", &ntrack); t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("px", px); t->SetBranchAddress("py", py); t->SetBranchAddress("pz", pz);

  std::vector<int> cent, v0m;                       // registered events only
  long nEv = t->GetEntries();
  for (long i = 0; i < nEv; ++i) {
    t->GetEntry(i);
    int c = 0, a = 0, cc = 0;
    for (int j = 0; j < ntrack && j < N; ++j) {
      if (charge[j] == 0) continue;
      double pt = std::hypot(px[j], py[j]);
      double p  = std::sqrt(pt*pt + pz[j]*pz[j]);
      if (p == std::fabs(pz[j])) continue;
      double eta = 0.5*std::log((p + pz[j])/(p - pz[j]));
      if (std::fabs(eta) < 0.8 && pt > 0.2 && pt < 3.0) ++c;
      if (eta > 2.8  && eta < 5.1)  ++a;
      if (eta > -3.7 && eta < -1.7) ++cc;
    }
    if (!(a > 0 && cc > 0)) continue;               // MB: V0A && V0C
    cent.push_back(c); v0m.push_back(a + cc);
  }
  const long nReg = cent.size();
  printf("events %ld, registered %ld\n", nEv, nReg);

  auto line = [&](const char* lab, int vcut, int ncut) {
    long n = 0; double s = 0;
    for (long i = 0; i < nReg; ++i)
      if (v0m[i] > vcut && cent[i] > ncut) { ++n; s += cent[i]; }
    if (n) printf("%-18s N=%8ld  %7.3f%% of reg  <nchCent>=%6.2f\n", lab, n, 100.*n/nReg, s/n);
    else   printf("%-18s N=0\n", lab);
  };

  printf("\n-- nchCent > 25 and V0M > cut --\n");
  for (int v = 20; v <= 120; v += 5) line(Form("V0M>%d", v), v, 25);

  printf("\n-- V0M > cut only (no nchCent cut) --\n");
  for (int v = 20; v <= 120; v += 10) line(Form("V0M>%d", v), v, -1);

  printf("\n-- top x%% of V0M (registered events), with nchCent > 25 --\n");
  std::vector<int> sv(v0m); std::sort(sv.begin(), sv.end());
  for (double x : {0.07, 0.2, 0.5, 1., 2., 5., 10., 13., 20.}) {
    int thr = sv[(size_t)((1. - x/100.) * nReg)];  // top x% starts at this V0M
    line(Form("top %.2f%% (V0M>%d)", x, thr - 1), thr - 1, 25);
  }

  printf("\n-- nchCent only (no V0M cut) --\n");
  for (int n : {20, 25, 30, 35, 40}) line(Form("nchCent>%d", n), -1, n);
}
