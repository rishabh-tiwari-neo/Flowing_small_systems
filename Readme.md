# Long-range correlations and flow signatures in small systems

Welcome! This study is on long-range correlations and flow signatures in small systems (pp collisions at 13 TeV) using PYTHIA8 with string shoving and rope hadronisation. It follows the analysis strategy of the ALICE papers on ultra-long-range correlations and partonic flow.

## Layout

```
main/
├── run_ALICE.sh            # Runs generation, merges, then fits and plots automatically
├── run_test.sh             # test run script (just a test script)
├── Makefile                # Auto-detects ROOT via root-config and PYTHIA (uses $PYTHIA8 if set, otherwise searches standard folders)
├── generator/
│   ├── gen_ALICE.cc        # PYTHIA8 generator for ALICE-like generation (the keyword is "Like")
│   └── stat_check.cc       # PYTHIA8 generator for checking statistics in the above given configuration (was needed pre-run, for stats)
├── post/
│   ├── analyser_ALICE.cc   # Prints event, multiplicity and trigger statistics of a merged ALICE file (used pre-run again)
│   ├── analyser.cc         # Same kind of pre-run stats, for the separate nchnmpi study (see data/data_nchnmpi below)
│   ├── post_ALICE.cc       # Template fit, V2Delta extraction, 3x2PC combination -> v2(pT) and plain v2
│   └── plot_ALICE.cc       # Makes all plots, including v2 and the NCQ-scaling plots, from the post-processing output
├── data/
│   ├── data_ALICE/         # ROOT files and v2_results_<class>.{root,csv} from gen_ALICE + post_ALICE
│   ├── data_nchnmpi/       # Output of the separate nchnmpi study
│   └── test_run/           # Small test-run output, incl. proton_stats.C
├── logs/
│   ├── log_ALICE/          # Seeds, per-task logs, progress files, summaries, fit tables (post_all.txt)
│   ├── log_nchnmpi/        # Logs for the nchnmpi study
│   └── test_run/           # Logs for a small test run
├── plots/
│   ├── plot_ALICE_HM/      # Finished HM-class plots (PDF and PNG)
│   ├── plot_ALICE_HM35/    # Finished HM35-class plots (PDF and PNG)
│   └── plot_nchnmpi/       # Plots for the nchnmpi study
└── Readme.md               # You are here
```

Note: `diag_ALICE.cc`, `dphi_ALICE.cc` and `slices_ALICE.cc` (the diagnostic tools described in "3b. Diagnosing a bad or noisy fit" below) are not currently part of this layout — add them to `post/` yourself if you need them; see that section for what they do and how to build them.

## What the generator does (`gen_ALICE.cc`)

- PYTHIA8: Monash tune (`Tune:pp = 14`), non-diffractive events, string shoving and rope hadronisation (`Ropewalk:gAmplitude = 10` by default, other Ropewalk settings follow Table I of Ortiz, Sahu and Bencedi).
- Per event, one pass over the final-state particles counts:
  - central charged multiplicity `nchCent` (|eta| < 0.8, 0.2 < pT < 3 GeV/c),
  - charged particles in the V0A and V0C windows (MB trigger and V0M estimator),
  - charged particles in FMD1,2 (1.7 < eta < 5.1) and FMD3 (-3.1 < eta < -1.7), used as reference particles (all pT, no PID).
- Identified triggers: pi, K, p (|eta| < 0.8) and K0S, Lambda (|y| < 0.5 and, at the same time, |eta| < 0.8), in 20 pT bins with edges
  0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 10.0 GeV/c.
- Event classes (an event can be in more than one):
  - HM: V0M > 157 and nchCent > 25 (about the top 0.07% of events)
  - HM35: V0M > 45 and nchCent > 25 (tuned so that mean nchCent is about 35)
  - LM: nchCent < 20 (the non-flow template)
- Same-event pairing (S), filled into 2D histograms of Delta-phi x Delta-eta (0.1-wide Delta-eta bins from 1.1 to 6.4):
  - TPC-FMD1,2 and TPC-FMD3 (identified trigger with a charged reference particle),
  - FMD1,2-FMD3 (two unidentified reference sets, all pT, no pT bins, no species dependence),
  - only pairs with 1.1 < |Delta-eta| < 6.4 are kept.
- Mixed-event pairing (B), same pairings and same 2D binning, as in the ALICE papers:
  - for each class the generator keeps the reference particles (eta, phi in FMD1,2 and FMD3) of the previous `POOL_DEPTH` events of that class (5, set at the top of `gen_ALICE.cc`),
  - the triggers (and FMD1,2 particles) of the current event are paired with the references of every pooled event, then the current event is added to the pool, so an event is never mixed with itself,
  - events are only mixed within the same class, which is the multiplicity-class matching used by ALICE. There is a single primary vertex in PYTHIA, so no vertex binning is needed.
- Normalisation counts are stored too: identified triggers per pT bin (`hNtrig`) and the number of FMD1,2 and FMD3 particles (`hNassoc`).
- Generator info: every run writes `generator_info.txt` into the data folder (next to the output files) with the PYTHIA version and the shoving strength g (`Ropewalk:gAmplitude`), read back from the PYTHIA settings that were used.

## What the post-processing does

- `analyser_ALICE`: prints event counts, class fractions (HM, LM, HM35), mean multiplicities per class, trigger yields, and same-event and mixed-event pair counts of a merged `gen_ALICE` file. The pT bins are read from the file. Again for testing and stats. 
- `post_ALICE`: for each class (HM, HM35), pairing, species and pT bin, builds the per-trigger yield Y(Delta-phi) = S/B (see Methodology), fits it with the LM template (template fit: weighted linear least squares, errors from the fit covariance matrix), and extracts V2Delta. The three pairings are then combined (3x2PC):
  - v2^2 = V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3), keeping the sign. No square root is taken, since v2^2 can come out negative when the combination fails.
  - plain v2 = sqrt(v2^2) is also written, defined only where v2^2 > 0. v2 is always positive, so it loses the sign of the underlying V2Delta — see Methodology, point 7.
  Any species can be analysed in coarser pT bins via the `REBIN` table in `post_ALICE.cc` (one entry per species, each an inclusive range of generator pT-bin indices to merge; an empty entry keeps the 20 generator bins). By default K0S and Lambda are merged into 0.5-0.8, 0.8-1.25, 1.25-1.75, 1.75-2.25, 2.25-2.75, 2.75-3.25, 3.25-4.0, 4.0-5.0, 5.0-6.0, 6.0-10 GeV/c, since they are sparse; pi, K and p keep the 20 generator bins by default. The merging sums the 2D S and B histograms and the trigger counts over the grouped generator bins before the S/B division, for both the class and the LM template.
  Each of the three pairings also sums Delta-eta only over its own kinematically-correct window (see Methodology, point 6), not the full 1.1-6.4 generator axis.
  Output per class: a ROOT file (V2Delta, v2^2 and v2 as TGraphErrors, one per species and pairing) and a CSV file (same numbers, plus the Delta-eta windows used, recorded as `#` header lines).
- `plot_ALICE`: for each class makes seven plots (PDF and PNG), all from `v2_results_<class>.root`, in one call — no separate tool needed:
  - `v2sq_vs_pt`, `v2_vs_pt`: v2^2 and plain v2 vs pT, all species.
  - `V2D_TPC_FMD12`, `V2D_TPC_FMD3`: the two pT-dependent pairings, all species.
  - `V2D_FMD12_FMD3`: the single reference point.
  - `NCQ_pT`, `NCQ_KET`: number-of-constituent-quark scaling, v2/nq vs pT/nq and v2/nq vs KET/nq (KET = sqrt(pT^2+m^2) - m), with nq = 2 for mesons (pi, K, K0S) and 3 for baryons (p, Lambda). Built from the `v2_<sp>` graphs `post_ALICE` writes. Points are kept only if v2/err(v2) exceeds `MIN_SIG_V2` in `plot_ALICE.cc` (default 4). Since v2 is always positive (see above), these two plots lose the sign of V2Delta the same way `v2_vs_pt` does.
- `diag_ALICE`, `dphi_ALICE`, `slices_ALICE`: diagnostic tools, not part of the normal pipeline. See "Diagnosing a bad or noisy fit" below.

## Methodology

This follows the ALICE papers (arXiv:2504.02359 and the Nature Commun. partonic-flow paper) and the template-fit setup of Ortiz, Sahu and Bencedi.

1. **Signal and background.** S is the same-event pair distribution in (Delta-phi, Delta-eta). B is the mixed-event distribution: the same pairing, but the two particles come from different events of the same class, so it carries the single-particle shapes and pair counting but no real correlation. As in the Nature paper, B is normalised with one constant: its Delta-phi average at the Delta-eta bin where B is largest. The corrected correlation is S/B, taken bin by bin in (Delta-phi, Delta-eta). In this perfect-acceptance generator B is expected to be nearly flat in Delta-phi, so S/B mainly reweights the Delta-eta slices.
2. **Per-trigger yield.** S is divided by the number of triggers: identified particles in the pT bin for TPC pairings, FMD1,2 particles for FMD1,2-FMD3. S/B is summed over Delta-eta (1.1 to 6.4) to give Y(Delta-phi), with B's statistical error included. No unit-area normalisation is applied.
3. **Template fit.** Y_HM(Delta-phi) = F x Y_LM(Delta-phi) + G x [1 + 2 V2D cos(2 Delta-phi) + 2 V3D cos(3 Delta-phi)], with F and G free. The LM class (nchCent < 20) is processed with the same S/B steps and acts as the non-flow template. V2Delta comes from a2/(2G); the error from the covariance of the linear least-squares fit.
4. **Three-subevent combination (3x2PC).** v2^2 = V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3).
5. **pT merging.** S and B are summed over the merged generator pT bins (with the trigger counts) before the division, for the HM class and for the LM template alike. By default this applies to K0S and Lambda only; any species can be merged via the `REBIN` table in `post_ALICE.cc`.
6. **Selection differences from ALICE.** HM is the top 0.07% of V0M with nchCent > 25 (mean nchCent about 54 here, against about 35 in ALICE). HM35 relaxes V0M so the mean is about 35; it is a different selection from ALICE's 0.07%. `gen_ALICE.cc` fills all three pairings into the same 1.1-6.4 Delta-eta axis, wider than the paper's 5.0-6.0 for the FMD pair, but the three pairings do not all cover 1.1-6.4 kinematically: TPC-FMD1,2 reaches 0.9-5.9, TPC-FMD3 reaches 0.9-3.9, and FMD1,2-FMD3 reaches 3.4-8.2 (from the trigger |eta|<0.8, FMD1,2 1.7-5.1 and FMD3 -3.1 to -1.7 acceptances). Outside its own range a pairing's Delta-eta slice is exactly empty, but just inside the edge the number of (eta_trig, eta_ref) combinations that can reach a Delta-eta that close to the kinematic boundary tapers toward zero, so a few slices there can have only a handful of raw pairs even with very high total statistics. Because the per-slice variance scales like 1/B^4, those thin edge slices can dominate a point's error or make F (the template scale) and a2 (the cos(2 Delta-phi) amplitude) indistinguishable in the fit (seen as a large |correlation coefficient rho| in `diag_ALICE`), which inflates the error and destabilises the central value — this hit K0S and Lambda hardest, at high pT, in the TPC-FMD3 pairing specifically. `post_ALICE.cc` now sums each pairing only over its own kinematic Delta-eta window, trimmed `DETA_EDGE_MARGIN` (default 0.4) off each edge to stay clear of the taper. Confirmed by running the fit before and after: the affected K0S/Lambda high-pT bins went from large, sign-flipping, or order-of-magnitude-outlier V2Delta values (and |rho| up to about 0.96) to smooth, small-error values consistent with their neighbours.
7. **Sign of V2Delta.** With shoving on and this g, V2Delta(TPC-FMD1,2) and V2Delta(TPC-FMD3) come out negative for every species once the deta-window fix above is applied, growing smoothly in magnitude with pT. This is consistent with string shoving's expected negative elliptic flow (repulsion pushes along the major axis of the initial geometry, the opposite of hydrodynamic flow) noted in the ALICE partonic-flow paper and in Ortiz, Sahu and Bencedi. v2^2 (point 4) is positive because two negative numbers divide by a positive V2Delta(FMD1,2-FMD3) in this class — so v2 = sqrt(v2^2), as plotted in `v2_vs_pt`, `NCQ_pT` and `NCQ_KET`, loses this sign and can look like ordinary positive flow. Check the signed V2Delta directly (`V2D_TPC_FMD12`/`V2D_TPC_FMD3`, or the CSV columns) before drawing conclusions from v2^2 or v2 alone. A Monash run with shoving off, through the same pipeline, is the direct test of whether this negative signal is shoving-driven: V2Delta should collapse to about zero there if so.

## User manual

### 1. Requirements

- ROOT (with `root-config` available in your shell).
- PYTHIA8. The Makefile uses `$PYTHIA8` if it is set (the folder that contains `include/` and `lib/`). Otherwise it tries `pythia8-config --prefix`, then `~/Pythia/pythia8*`, `~/pythia8*`, `/opt/pythia8*` (newest version first) and `/usr/local`. Each build prints the folder and version it used. If nothing is found, set `export PYTHIA8=/path/to/pythia8`.
- A C++17 compiler (g++).

### 2. Generate events, fit and plot (one command)

From the project root (`main/`):

```
chmod +x run_ALICE.sh     # first time only
./run_ALICE.sh
```

The script:

1. compiles `generator/gen_ALICE.cc`,
2. launches the parallel tasks with distinct seeds (saved in `logs/log_ALICE/seeds.log`) and shows a progress bar per task,
3. merges the task files with `hadd` into `data/data_ALICE/merged.root` (stops if the merge fails),
4. recompiles `post_ALICE` and `plot_ALICE`,
5. runs `post_ALICE` (HM and HM35) and saves the printed fit tables to `logs/log_ALICE/post_all.txt`,
6. runs `plot_ALICE` (HM and HM35).

Settings at the top of `run_ALICE.sh`:

- `TOTAL`: total number of events.
- `NCORES`: number of parallel tasks.

Note: a new run writes into the output folder and removes the old task output there. To keep earlier results, copy or rename the folder first.

### 3. Run the post-processing by hand

Compile (first time, or after editing the programs):

In the main folder run
```
make 
```

Run from the project root, no arguments needed:

```
./post/post_ALICE         # fits HM and HM35 from data/data_ALICE/merged.root
./post/plot_ALICE         # makes all seven plots (incl. v2 and NCQ scaling) for HM and HM35
./post/plot_ALICE HM35    # plots one class only
```

Files and folders used (set as constants at the top of the programs):

- `post_ALICE` reads `data/data_ALICE/merged.root` and writes `data/data_ALICE/v2_results_<class>.root` and `.csv`.
- `plot_ALICE` reads `data/data_ALICE/v2_results_<class>.root` and writes into `plots/plot_ALICE_<class>/` (the folder is created if missing).
- The class list (HM, HM35) is defined in both `post_ALICE.cc` and `plot_ALICE.cc`; add a class in both.

### 3b. Diagnosing a bad or noisy fit

If a bin's V2Delta looks like an outlier (large error, sign flip, or a value far from its neighbours), these three tools look inside the fit that produced it. None of them are part of the normal pipeline; build and run them individually, same pattern as `post_ALICE`:

```
./post/diag_ALICE                                    # per class/species/bin: nTrig, nTrigLM, npts, chi2/ndf, rho (F-a2 correlation), both V2Delta
./post/dphi_ALICE   <class> <species> <b0> <b1> [pairing]   # Y(dphi) vs Y_LM(dphi), 20 points, for one chosen bin
./post/slices_ALICE <class> <species> <b0> <b1> [pairing]   # inside one dphi bin, the individual deta slices that get summed
```

`<b0> <b1>` are inclusive generator pT-bin indices (0-19); `[pairing]` is `TPC_FMD12` (default) or `TPC_FMD3`. A large `|rho|` (close to +-1) means the fit cannot separate the template scale from the flow amplitude in that bin; a high chi2/ndf with small `|rho|` means the template simply does not describe the shape there. `slices_ALICE`'s `worst%` column shows whether a single thin deta slice is driving a point's error — this is what the Delta-eta-window fix (Methodology, point 6) addresses.

Both `diag_ALICE.cc` and `post_ALICE.cc` must be kept in sync by hand for the `REBIN` table and the per-pairing Delta-eta windows, since `diag_ALICE` is a standalone copy, not a shared include.

### 4. Monitoring a run

- Progress files: `logs/log_ALICE/progress_task<seed>.txt` (one number per task: events done).
- Per-task logs and the merged summary: `logs/log_ALICE/`.

### 5. How to read the results

- The main result is V2Delta per species and pT bin, for TPC-FMD1,2 and TPC-FMD3, plus the single V2Delta(FMD1,2-FMD3) per class.

### 6. Things to keep in mind

- Record the PYTHIA version used for each run (the build prints it). Results from different versions (for example 8.312, 8.316, 8.317) can differ slightly even with the same seeds.
- Output files from before the mixed-event step (1D Delta-phi histograms) are not compatible with the current `post_ALICE`, and old and new files cannot be merged.
- The mixed-event step multiplies the pairing work by roughly the pool depth. Check the run time on a small run first.
- Errors are statistical only (from the fit covariance). They treat the LM template as exact and the pairs as independent, so they are somewhat too small (chi2/ndf of the fits is typically about 1.2 to 2.2). There are no systematic uncertainties yet.
- The cosines in the fit are evaluated at the Delta-phi bin centres, which biases V2 low by about 1.6%.
- The shoving strength g is set inside `gen_ALICE.cc` (`Ropewalk:gAmplitude`) and is written, with the PYTHIA version, to `generator_info.txt` in the data folder. Keep different g values in different folders anyway.
- Histogram names must stay identical in `gen_ALICE.cc` and `post_ALICE.cc`; species names must also match in `plot_ALICE.cc`.
- Fits fail (printed as "fit failed") when a bin has too few pairs. This is expected in small test runs, especially for HM at high pT.
- The HM class's V2Delta(FMD1,2-FMD3) reference point can come out consistent with zero (sign-ambiguous) even in a large run, since HM is a narrow top-0.07% selection. When that happens v2^2 for every species in that class is unreliable (dividing by a near-zero, sign-ambiguous denominator), even though the per-pairing V2Delta(TPC-FMD1,2)/V2Delta(TPC-FMD3) values for that class are still fine on their own. Check V2Delta(FMD1,2-FMD3)'s own significance before trusting a class's v2^2 column.
- `NCQ_pT` and `NCQ_KET` previously looked identical when made by an earlier, standalone NCQ script, even though its KET = sqrt(pT^2+m^2) - m formula differs from plain pT. That script is now folded into `plot_ALICE.cc` (`drawNCQ`), with the KET formula unchanged — but `drawNCQ`'s call sites had an extra `outdir` argument that doesn't match its signature, so the file failed to compile at all until that was fixed (fixed; `drawNCQ(f, kind, cls, g)`, no `outdir`). Since `plot_ALICE` with the NCQ plots never successfully built before that fix, the two issues are unrelated: the identical-plots report is from the old standalone script, not this one. Still worth a pixel-by-pixel diff of the two new PNGs on the next successful build, since the root cause of the original report was never confirmed.