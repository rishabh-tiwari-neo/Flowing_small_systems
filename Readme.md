# V2delta: long-range correlations and flow signatures in small systems

Welcome! This repository studies long-range correlations and flow signatures in small systems (pp collisions at 13 TeV) using PYTHIA8 with string shoving and rope hadronisation. It follows the analysis strategy of the ALICE papers on ultra-long-range correlations and partonic flow.

## Layout

```
main/
├── run_ALICE.sh            # Runs generation (parallel tasks) and merges the output
├── run.sh, run_test.sh     # Older / test run scripts
├── Makefile                # Auto-detects ROOT via root-config; requires the PYTHIA8 environment variable
├── generator/
│   └── gen_ALICE.cc        # PYTHIA8 generator + event loop (see below)
├── post/
│   ├── analyser_ALICE.cc   # Prints event, multiplicity and trigger statistics of a merged file
│   ├── post_ALICE.cc       # Template fit, V2Delta extraction, 3x2PC combination -> v2(pT)
│   └── plot_ALICE.cc       # Makes the plots from the post-processing output
├── data/                   # ROOT files (one sub-folder per generation run)
├── logs/                   # Seeds, per-task logs, progress files, summaries
├── plots/                  # Finished plots (one sub-folder per post-processing run)
└── Readme.md               # You are here
```

## What the generator does (`gen_ALICE.cc`)

- PYTHIA8: Monash tune (`Tune:pp = 14`), non-diffractive events, string shoving and rope hadronisation (`Ropewalk:gAmplitude = 10` by default, other Ropewalk settings follow Table I of Ortiz, Sahu and Bencedi).
- Per event, one pass over the final-state particles counts:
  - central charged multiplicity `nchCent` (|eta| < 0.8, 0.2 < pT < 3 GeV/c),
  - charged particles in the V0A and V0C windows (MB trigger and V0M estimator),
  - charged particles in FMD1,2 (1.7 < eta < 5.1) and FMD3 (-3.1 < eta < -1.7), used as reference particles.
- Identified triggers: pi, K, p (|eta| < 0.8) and K0S, Lambda (|y| < 0.5), in 16 pT bins from 0.5 to 10 GeV/c.
- Event classes (an event can be in more than one):
  - HM: V0M > 157 and nchCent > 25
  - HM35: V0M > 45 and nchCent > 25 (tuned so that mean nchCent is about 35)
  - LM: nchCent < 20 (the non-flow template)
- Same-event pairing, filling Delta-phi histograms directly (no mixed-event step for now):
  - TPC-FMD1,2 and TPC-FMD3 (identified trigger with a charged reference particle),
  - FMD1,2-FMD3 (two unidentified reference sets),
  - only pairs with 1.1 < |Delta-eta| < 6.4 are kept.

## What the post-processing does

- `analyser_ALICE`: prints event counts, class fractions, mean multiplicities, trigger yields and pair counts.
- `post_ALICE`: for each pairing, species and pT bin, fits the HM (or HM35) Delta-phi projection with the LM template (template fit), extracts V2Delta, and combines the three pairings into v2(pT) = sqrt(V2D(TPC-FMD1,2) x V2D(TPC-FMD3) / V2D(FMD1,2-FMD3)). Output: a ROOT file and a CSV file.
- `plot_ALICE`: v2 vs pT, v2^2 vs pT, and V2Delta of the pairings.

## User manual

### 1. Requirements

- ROOT (with `root-config` available in your shell).
- PYTHIA8 (the run script uses the include and lib paths of your PYTHIA8 installation; the Makefile needs the `$PYTHIA8` environment variable).
- A C++17 compiler (g++).

### 2. Generate events

From the project root (`main/`):

```
chmod +x run_ALICE.sh     # first time only
./run_ALICE.sh
```

The script compiles `generator/gen_ALICE.cc`, launches the parallel tasks with distinct seeds (saved in `logs/log_ALICE/seeds.log`), shows a progress bar per task, and merges the task files with `hadd` into `data/data_ALICE/merged.root`.

Settings at the top of `run_ALICE.sh`:

- `TOTAL`: total number of events.
- `NCORES`: number of parallel tasks.

Note: a new run writes into the output folder and removes the old output there. To keep earlier results, copy or rename the folder first, and never merge files made with different pT binnings or different settings.

### 3. Compile the post-processing programs (first time, or after editing them)

```
g++ -O2 -o post/analyser_ALICE post/analyser_ALICE.cc `root-config --cflags --libs`
g++ -O2 -o post/post_ALICE     post/post_ALICE.cc     `root-config --cflags --libs`
g++ -O2 -o post/plot_ALICE     post/plot_ALICE.cc     `root-config --cflags --libs`
```

### 4. Run the post-processing

```
./post/analyser_ALICE data/data_ALICE/merged.root
./post/post_ALICE data/data_ALICE/merged.root data/data_ALICE/v2_results_HM35 HM35
./post/plot_ALICE data/data_ALICE/v2_results_HM35.root plots/plot_ALICE_HM35
```

- The third argument of `post_ALICE` is the event class (`HM` or `HM35`). The LM class is always used as the template.
- `post_ALICE` writes `<out_base>.root` and `<out_base>.csv`.
- `plot_ALICE` writes PDF and PNG files into the given plots folder.
- Use a different output name for each class or run so earlier results are not overwritten.

### 5. Monitoring a run

- Progress files: `logs/log_ALICE/progress_task<seed>.txt` (one number per task: events done).
- Per-task logs and the merged summary: `logs/log_ALICE/`.

### 6. Things to keep in mind

- The shoving strength g is set inside `gen_ALICE.cc` (`Ropewalk:gAmplitude`) and is not yet stored in the output file. Keep different g values in different folders and note the value in your log.
- Histogram names must stay identical in `gen_ALICE.cc` and `post_ALICE.cc`.
- `analyser_ALICE` still needs to be updated for the 16 pT bins and the current histogram names.
