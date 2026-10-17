# Long-range correlations and flow signatures in small systems

Welcome! This study is on long-range correlations and flow signatures in small systems (pp collisions at 13 TeV) using PYTHIA8 with string shoving and rope hadronisation. It follows the analysis strategy of the ALICE papers on ultra-long-range correlations and partonic flow.

## Layout

```
main/
├── run_ALICE.sh            # Runs generation, merges, then fits and plots automatically
├── run_test.sh             # test run script
├── Makefile                # Auto-detects ROOT via root-config; requires the PYTHIA8 environment variable
├── generator/
│   └── gen_ALICE.cc        # PYTHIA8 generator for ALICE-like generation
│   └── stat_.cc            # PYTHIA8 generator for checking statistics in the above given configuration (need for pre-run)
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
- Identified triggers: pi, K, p (|eta| < 0.8) and K0S, Lambda (|y| < 0.5), in 20 pT bins with edges
  0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 10.0 GeV/c.
- Event classes (an event can be in more than one):
  - HM: V0M > 157 and nchCent > 25 (about the top 0.07% of events)
  - HM35: V0M > 45 and nchCent > 25 (tuned so that mean nchCent is about 35)
  - LM: nchCent < 20 (the non-flow template)
- Same-event pairing, filling Delta-phi histograms directly (no mixed-event step; the generator has perfect acceptance):
  - TPC-FMD1,2 and TPC-FMD3 (identified trigger with a charged reference particle),
  - FMD1,2-FMD3 (two unidentified reference sets, all pT, no pT bins: one number per class),
  - only pairs with 1.1 < |Delta-eta| < 6.4 are kept.

## What the post-processing does

- `analyser_ALICE`: prints event counts, class fractions, mean multiplicities, trigger yields and pair counts.
- `post_ALICE`: for each class (HM, HM35), pairing, species and pT bin, fits the Delta-phi projection with the LM template (template fit: weighted linear least squares, errors from the fit covariance matrix), and extracts V2Delta. The three pairings are then combined (3x2PC):
  - v2^2 = V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3), keeping the sign,
  - v2 = sign(v2^2) x sqrt(|v2^2|).
  Output per class: a ROOT file and a CSV file.
- `plot_ALICE`: for each class makes five plots (PDF and PNG): v2 vs pT, v2^2 vs pT, V2Delta(TPC-FMD1,2) vs pT, V2Delta(TPC-FMD3) vs pT, and the single V2Delta(FMD1,2-FMD3) point.

## User manual

### 1. Requirements

- ROOT (with `root-config` available in your shell).
- PYTHIA8 (the run script uses the include and lib paths of your PYTHIA8 installation; the Makefile needs the `$PYTHIA8` environment variable).
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

In the main file run
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

The analyser is still run on its own:

```
./post/analyser_ALICE data/data_ALICE/merged.root
```

### 4. Monitoring a run

- Progress files: `logs/log_ALICE/progress_task<seed>.txt` (one number per task: events done).
- Per-task logs and the merged summary: `logs/log_ALICE/`.

### 5. How to read the results

- The main result is V2Delta per species and pT bin, for TPC-FMD1,2 and TPC-FMD3, plus the single V2Delta(FMD1,2-FMD3) per class.
- v2 from 3x2PC is only meaningful where both TPC V2Delta values have the same sign as the reference. A negative v2^2 means the combination failed (statistical noise or factorisation breaking), not negative flow. A positive v2 from two negative V2Delta is not flow either.
- When the reference is consistent with zero, or the TPC V2Delta have the opposite sign, report V2Delta per species with errors instead of v2.

### 6. Things to keep in mind

- Errors are statistical only (from the fit covariance). They treat the LM template as exact and the pairs as independent, so they are somewhat too small (typical chi2/ndf is about 1.3 to 1.6 for the species fits). There are no systematic uncertainties yet.
- The cosines in the fit are evaluated at the Delta-phi bin centres, which biases V2 low by about 1.6%.
- The shoving strength g is set inside `gen_ALICE.cc` (`Ropewalk:gAmplitude`) and is not yet stored in the output file. Keep different g values in different folders and note the value in your log.
- Histogram names must stay identical in `gen_ALICE.cc` and `post_ALICE.cc`; species names must also match in `plot_ALICE.cc`.
- Fits fail (printed as "fit failed") when a bin has too few pairs. This is expected in small test runs, especially for HM at high pT.
- `analyser_ALICE` still needs to be updated for the 20 pT bins and the current histogram names.
