// Run from data/test_run/:  root -l -b -q proton_stats.C
// Proton statistics in HM events: counts per pT bin, equal-statistics edges,
// and a rough estimate of the pair statistics available for V2Delta.
void proton_stats() {
    TFile *f = TFile::Open("merged.root");
    TTree *t = (TTree*)f->Get("events");

    // V0M threshold, top 0.07% of registered events (same as the analyser)
    t->Draw("nchV0A+nchV0C>>hV(400,0,400)", "nchV0A>0&&nchV0C>0", "goff");
    TH1 *hV = (TH1*)gDirectory->Get("hV");
    double q = 0.9993, thr = 0;
    hV->GetQuantiles(1, &thr, &q);
    thr = std::round(thr * 10.0) / 10.0;

    TString hm   = Form("(nchV0A+nchV0C)>%.1f && nchCent>25", thr);
    TString prot = "abs(pid)==2212 && abs(eta)<0.8 && pt>0.5 && pt<10";

    long nEvHM = t->GetEntries(hm);

    // mean number of forward reference particles per HM event
    t->Draw("nchFMD12>>h12(2000,0,2000)", hm, "goff");
    double m12 = ((TH1*)gDirectory->Get("h12"))->GetMean();
    t->Draw("nchFMD3>>h3(2000,0,2000)", hm, "goff");
    double m3 = ((TH1*)gDirectory->Get("h3"))->GetMean();
    double refs = m12 + m3;

    // protons in HM, 0.1 GeV/c bins from 0.5 to 10
    t->Draw("pt>>hp(95,0.5,10)", hm + "&&" + prot, "goff");
    TH1 *hp = (TH1*)gDirectory->Get("hp");
    double nProt = hp->Integral();

    printf("V0M threshold = %.1f\n", thr);
    printf("HM events     = %ld\n", nEvHM);
    printf("protons (p+pbar, |eta|<0.8, 0.5<pT<10) in HM = %.0f\n", nProt);
    printf("mean forward refs per HM event: FMD1,2 = %.1f  FMD3 = %.1f  sum = %.1f\n\n",
           m12, m3, refs);

    // ---- fixed pT edges ----
    const int NE = 11;
    double e[NE] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10};
    printf("=== fixed pT bins ===\n");
    printf("%-10s %10s %14s %16s\n", "pT bin", "protons", "pairs(approx)", "sigma(V2D)~");
    for (int i = 0; i < NE - 1; ++i) {
        int b1 = hp->FindBin(e[i] + 0.05), b2 = hp->FindBin(e[i+1] - 0.05);
        double n = hp->Integral(b1, b2);
        double pairs = n * refs;
        double sig = pairs > 0 ? 1.0 / std::sqrt(2.0 * pairs) : 0;
        printf("%4.1f-%-4.1f %10.0f %14.0f %16.4f\n", e[i], e[i+1], n, pairs, sig);
    }

    // ---- equal-statistics edges ----
    const int NQ = 6;
    double probs[NQ+1], qs[NQ+1];
    for (int i = 0; i <= NQ; ++i) probs[i] = (double)i / NQ;
    hp->GetQuantiles(NQ + 1, qs, probs);
    printf("\n=== %d equal-statistics bins (%.0f protons each) ===\n", NQ, nProt / NQ);
    for (int i = 0; i < NQ; ++i)
        printf("  %5.2f - %5.2f GeV/c\n", qs[i], qs[i+1]);

    printf("\nsigma(V2D)~ = 1/sqrt(2*pairs): the error of <cos 2dphi> if pairs were independent.\n");
    printf("Pairs share reference particles and the dEta gap removes some, so the real error is larger.\n");
}
