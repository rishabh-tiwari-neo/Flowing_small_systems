// compare_hists.C
// Compares every histogram of file A with the same-named histogram of file B, bin by bin
// (all bins including under/overflow). Mixed-event histograms (hMix_*) are reported separately,
// because they depend on where the pool resets (block boundaries); everything else must be
// exactly equal for ANY number of cores.
//
// Usage:  root -l -b -q 'compare_hists.C("A.root","B.root")'
// Read-only: it does not write or change any file.

#include <cstdio>
#include <cstring>
#include <cmath>
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TCollection.h"
#include <algorithm>

void compare_hists(const char* fileA, const char* fileB)
{
    TFile* a = TFile::Open(fileA, "READ");
    TFile* b = TFile::Open(fileB, "READ");
    if (!a || a->IsZombie() || !b || b->IsZombie()) { printf("Cannot open the files\n"); return; }

    int nOther = 0, nOtherBad = 0, nMix = 0, nMixBad = 0, nMissing = 0;
    double maxOther = 0, mixSumA = 0, mixSumB = 0, maxMixRel = 0;

    TIter next(a->GetListOfKeys());
    TKey* k;
    while ((k = (TKey*)next())) {
        TObject* oa = k->ReadObj();
        TH1* ha = dynamic_cast<TH1*>(oa);
        if (!ha) continue;
        TH1* hb = dynamic_cast<TH1*>(b->Get(ha->GetName()));
        if (!hb || hb->GetNcells() != ha->GetNcells()) {
            printf("MISSING or different binning in B: %s\n", ha->GetName());
            ++nMissing;
            continue;
        }
        const bool isMix = (std::strncmp(ha->GetName(), "hMix_", 5) == 0);
        double maxd = 0, sa = 0, sb = 0;
        for (int i = 0; i < ha->GetNcells(); ++i) {
            const double va = ha->GetBinContent(i), vb = hb->GetBinContent(i);
            maxd = std::max(maxd, std::fabs(va - vb));
            sa += va; sb += vb;
        }
        if (isMix) {
            ++nMix; mixSumA += sa; mixSumB += sb;
            if (maxd > 0) ++nMixBad;
            if (sa > 0) maxMixRel = std::max(maxMixRel, std::fabs(sa - sb) / sa);
        } else {
            ++nOther;
            if (maxd > 0) { ++nOtherBad; maxOther = std::max(maxOther, maxd); printf("DIFFERENT: %s (max |diff| = %g)\n", ha->GetName(), maxd); }
        }
    }

    // histograms present in B but not in A
    TIter nextB(b->GetListOfKeys());
    while ((k = (TKey*)nextB())) {
        if (!a->Get(k->GetName())) { printf("Only in B: %s\n", k->GetName()); ++nMissing; }
    }

    printf("\n--- non-mixed histograms: %d compared, %d differ (largest |diff| = %g)\n", nOther, nOtherBad, maxOther);
    printf("--- mixed histograms (hMix_*): %d compared, %d differ; total content A = %.0f, B = %.0f; largest per-histogram relative difference of the total = %g\n",
           nMix, nMixBad, mixSumA, mixSumB, maxMixRel);
    printf("--- missing / mismatched: %d\n", nMissing);
    printf("%s\n", (nOtherBad == 0 && nMissing == 0) ? "RESULT: non-mixed histograms are identical." : "RESULT: NOT identical, see lines above.");
    a->Close(); b->Close();
}
