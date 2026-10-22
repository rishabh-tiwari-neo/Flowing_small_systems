# Long-range correlations and flow signatures in small systems

Welcome! This study is on long-range correlations and flow signatures in small systems (pp collisions at 13 TeV) using PYTHIA8 with string shoving and rope hadronisation. It follows the analysis strategy of the ALICE papers on ultra-long-range correlations and partonic flow.

## Layout

```
main/
├── run_ALICE.sh            # Runs generation, merges, then fits and plots automatically
├── run_test.sh             # test run script
├── Makefile                # Auto-detects ROOT via root-config and PYTHIA (uses $PYTHIA8 if set, otherwise searches standard folders)
├── generator/
│   └── gen_ALICE.cc        # PYTHIA8 generator for ALICE-like generation
│   └── stat_check.cc       # PYTHIA8 generator for checking statistics in the above given configuration (need for pre-run)
├── post/
│   ├── analyser_ALICE.cc   # Prints event, multiplicity and trigger statistics of a merged file
│   ├── post_ALICE.cc       # Template fit, V2Delta extraction, 3x2PC combination -> v2(pT)
│   └── plot_ALICE.cc       # Makes the plots from the post-processing output
├── data/                   # ROOT files (one sub-folder per generation run)
├── logs/                   # Seeds, per-task logs, progress files, summaries, fit tables
├── plots/                  # Finished plots 
└── Readme.md               # You are here
```

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

- `analyser_ALICE`: prints event counts, class fractions (HM, LM, HM35), mean multiplicities per class, trigger yields, and same-event and mixed-event pair counts of a merged `gen_ALICE` file. The pT bins are read from the file.
- `post_ALICE`: for each class (HM, HM35), pairing, species and pT bin, builds the per-trigger yield Y(Delta-phi) = S/B (see Methodology), fits it with the LM template (template fit: weighted linear least squares, errors from the fit covariance matrix), and extracts V2Delta. The three pairings are then combined (3x2PC):
  - v2^2 = V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3), keeping the sign. No square root is taken, since v2^2 can come out negative when the combination fails.
  K0S and Lambda are sparse, so for them the pT bins are merged before the S/B division: the 2D S and B histograms and the trigger counts of the generator bins are summed, giving the bins 0.5-0.8, 0.8-1.25, 1.25-1.75, 1.75-2.25, 2.25-2.75, 2.75-3.25, 3.25-4.0, 4.0-5.0, 5.0-6.0, 6.0-10 GeV/c (table `ptGroups` in `post_ALICE.cc`). pi, K and p keep the 20 generator bins.
  Output per class: a ROOT file and a CSV file.
- `plot_ALICE`: for each class makes four plots (PDF and PNG): v2^2 vs pT, V2Delta(TPC-FMD1,2) vs pT, V2Delta(TPC-FMD3) vs pT, and the single V2Delta(FMD1,2-FMD3) point.

## Methodology

This follows the ALICE papers (arXiv:2504.02359 and the Nature Commun. partonic-flow paper) and the template-fit setup of Ortiz, Sahu and Bencedi.

1. **Signal and background.** S is the same-event pair distribution in (Delta-phi, Delta-eta). B is the mixed-event distribution: the same pairing, but the two particles come from different events of the same class, so it carries the single-particle shapes and pair counting but no real correlation. As in the Nature paper, B is normalised with one constant: its Delta-phi average at the Delta-eta bin where B is largest. The corrected correlation is S/B, taken bin by bin in (Delta-phi, Delta-eta). In this perfect-acceptance generator B is expected to be nearly flat in Delta-phi, so S/B mainly reweights the Delta-eta slices.
2. **Per-trigger yield.** S is divided by the number of triggers: identified particles in the pT bin for TPC pairings, FMD1,2 particles for FMD1,2-FMD3. S/B is summed over Delta-eta (1.1 to 6.4) to give Y(Delta-phi), with B's statistical error included. No unit-area normalisation is applied.
3. **Template fit.** Y_HM(Delta-phi) = F x Y_LM(Delta-phi) + G x [1 + 2 V2D cos(2 Delta-phi) + 2 V3D cos(3 Delta-phi)], with F and G free. The LM class (nchCent < 20) is processed with the same S/B steps and acts as the non-flow template. V2Delta comes from a2/(2G); the error from the covariance of the linear least-squares fit.
4. **Three-subevent combination (3x2PC).** v2^2 = V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3).
5. **pT merging for K0S and Lambda.** S and B are summed over the merged generator pT bins (with the trigger counts) before the division, for the HM class and for the LM template alike.
6. **Selection differences from ALICE.** HM is the top 0.07% of V0M with nchCent > 25 (mean nchCent about 54 here, against about 35 in ALICE). HM35 relaxes V0M so the mean is about 35; it is a different selection from ALICE's 0.07%. The Delta-eta window is 1.1 to 6.4, wider than the paper's 5.0 to 6.0 for the FMD pair.

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
./post/plot_ALICE         # plots HM and HM35
./post/plot_ALICE HM35    # plots one class only
```

Files and folders used (set as constants at the top of the programs):

- `post_ALICE` reads `data/data_ALICE/merged.root` and writes `data/data_ALICE/v2_results_<class>.root` and `.csv`.
- `plot_ALICE` reads `data/data_ALICE/v2_results_<class>.root` and writes into `plots/plot_ALICE_<class>/` (the folder is created if missing).
- The class list (HM, HM35) is defined in both `post_ALICE.cc` and `plot_ALICE.cc`; add a class in both.

### 4. Monitoring a run

- Progress files: `logs/log_ALICE/progress_task<seed>.txt` (one number per task: events done).
- Per-task logs and the merged summary: `logs/log_ALICE/`.

### 5. How to read the results

- The main result is V2Delta per species and pT bin, for TPC-FMD1,2 and TPC-FMD3, plus the single V2Delta(FMD1,2-FMD3) per class.
- v2^2 is only meaningful where V2D(TPC-FMD1,2), V2D(TPC-FMD3) and V2D(FMD1,2-FMD3) have the same sign. A negative v2^2 means the combination failed (noise or factorisation breaking), not negative flow.

### 6. Things to keep in mind

- Record the PYTHIA version used for each run (the build prints it). Results from different versions (for example 8.312, 8.316, 8.317) can differ slightly even with the same seeds.
- Output files from before the mixed-event step (1D Delta-phi histograms) are not compatible with the current `post_ALICE`, and old and new files cannot be merged.
- The mixed-event step multiplies the pairing work by roughly the pool depth. Check the run time on a small run first.
- Errors are statistical only (from the fit covariance). They treat the LM template as exact and the pairs as independent, so they are somewhat too small (chi2/ndf of the fits is typically about 1.2 to 2.2). There are no systematic uncertainties yet.
- The cosines in the fit are evaluated at the Delta-phi bin centres, which biases V2 low by about 1.6%.
- The shoving strength g is set inside `gen_ALICE.cc` (`Ropewalk:gAmplitude`) and is written, with the PYTHIA version, to `generator_info.txt` in the data folder. Keep different g values in different folders anyway.
- Files generated with and without the K0S/Lambda |eta| < 0.8 cut (it changes the Lambda 0.5-0.6 GeV/c bin) must not be merged.
- Histogram names must stay identical in `gen_ALICE.cc` and `post_ALICE.cc`; species names must also match in `plot_ALICE.cc`.
- Fits fail (printed as "fit failed") when a bin has too few pairs. This is expected in small test runs, especially for HM at high pT.

## Known behaviour of the results (250M events, S/B)

- Reference: V2D(FMD1,2-FMD3) = 2.52e-4 +- 8.6e-6 for HM35 (chi2/ndf 1.20) and 1.5e-4 +- 5.0e-5 for HM (1.46).
- Pions in HM35 have a negative V2Delta above about 0.8 GeV/c. It comes almost entirely from the Delta-eta 1.1-2.0 slice (short-range non-flow next to the Delta-eta cut) where the template fit also has a poor chi2/ndf; the other Delta-eta slices are close to zero. The resulting v2^2 is positive only because both TPC V2Delta values are negative, so it is not a flow signal.
- K0S and Lambda have few triggers at high pT. In S/B the B cells with few counts get a large weight (c/B) and give large and irregular errors. Merging the pT bins fixed this for Lambda; for K0S the values at 3.25-5 GeV/c are still dominated by noise and should not be read.
- v2^2 plots: the large error bars of the noisy K0S and Lambda points stretch the y-axis.
