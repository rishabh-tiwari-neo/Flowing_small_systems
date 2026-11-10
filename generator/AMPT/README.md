# AMPT v1.26t9b-v2.26t9b for pp collisions — fast ROOT edition

The official AMPT model `v1.26t9b-v2.26t9b` (Z.-W. Lin et al.,
<https://myweb.ecu.edu/linz/ampt/>), set up for **proton-proton collisions**
in the string-melting mode, with two additions:

* it writes a **ROOT tree `t`** directly while it runs (impact parameter,
  participants, nucleon positions and the final-state particles, see
  [Output tree `t`](#output-tree-t)); the text files `ana/ampt.dat` and
  `ana/zpc.dat` are optional;
* **exact speed-ups** of the transport code: the same seeds give the same
  events as the official code, only faster (see [Changes](#changes-to-the-official-code)).

The changes in the Fortran sources are marked with `cfast` comments; the ROOT
output is in `AmptRootIO.cxx`. The meaning of every line of `input.ampt` is
explained in the `readme` of the official release (link above).

---

## Files

| File | Content |
| :--- | :--- |
| `main.f`, `amptsub.f`, `art1f.f`, `hijing1.383_ampt.f`, `hipyset1.35.f`, `linana.f`, `zpc.f` | AMPT Fortran sources |
| `AmptRootIO.cxx` | ROOT output (tree `t`) |
| `Makefile` | build (`make`, `make pgo`, `make clean`) |
| `input.ampt` | input of the run: pp at √s = 5.36 TeV |
| `examples/input.pp_*GeV.ampt` | pp inputs at √s = 5.36, 13.0 and 13.6 TeV |
| `exec` | one run with random seeds |
| `run_fast_ampt.sh` | parallel run on the local machine, merged into one ROOT file |
| `production/runAMPT_fast.sh` | large productions: many files of fixed size, resumable |
| `ana/` | output directory (created by `exec` and the run scripts) |

---

## Build

Needs `gfortran`, `g++` and ROOT (`root-config` in the `PATH`).

```bash
make -j$(nproc)      # or: make pgo   (~10% faster, trains on input.ampt)
```

* `-march=native`: the binary is made for the CPU type it was built on; build it
  again on another kind of machine.
* If the compiler has no LTO linker plugin, build with `make LTO=`.
* The build is IEEE-exact: no `-ffast-math`, and `-ffp-contract=off` so that
  `-march=native` does not fuse `a*b+c` into FMA instructions (they round
  differently and would make every event different from the official code).

---

## Physics settings (`input.ampt`)

* **System**: p+p (`PROJ = P`, `TARG = P`, `IAP = IZP = IAT = IZT = 1`), √s = 5.36 TeV (`EFRM`, `FRAME = CMS`)
* **Impact parameter**: minimum bias, `BMIN = 0`, `BMAX = 3 fm`
* **AMPT mode**: string melting (`ISOFT = 4`)
* **Lund string fragmentation**: a = 0.30 (`PARJ(41)`), b = 0.15 GeV⁻² (`PARJ(42)`)
* **Parton cascade**: screening mass μ = 2.265 fm⁻¹, αs = 0.33 (σ ≈ 3 mb)
* **Hadron cascade**: `NTMAX = 150` timesteps of `DT = 0.2 fm/c`
* **Popcorn mechanism**: on (`PARJ(5) = 1.0`)
* **Decays at the end of the hadron cascade**: K0S on, φ on, π0 off. So K0S and φ
  do not appear in the output, only their decay products; switch these flags
  to 0 to keep them.
* **Reaction plane**: random orientation (last line, so `psiRP` is random)
* **Seeds**: read at run time (`ihjsed = 11`), see [Seeds](#seeds)

Other energies: copy an input of `examples/` to `input.ampt`, e.g.
`cp examples/input.pp_13600GeV.ampt input.ampt`.

---

## Running

### One run

```bash
./exec               # NEVNT events of input.ampt, random seeds -> ana/ampt.root, log in nohup.out
```

or directly, with your own seeds:

```bash
mkdir -p ana
printf '12345\n67890\n' | ./ampt > nohup.out     # -> ana/ampt.root
```

### Parallel on the local machine

```bash
# ./run_fast_ampt.sh [N_WORKERS] [TOTAL_EVENTS] [OUTPUT] [EVENTS_PER_SHARD]
./run_fast_ampt.sh 8 100000 ana/pp_5360GeV_100k.root
```

The events are split into shards (default `TOTAL/(4·N_WORKERS)` events each)
that the workers take one after another. Every shard gets its own random seeds;
the shards are merged with `hadd`, and the merged file lists all seeds in
`source` and the summed event counts in `tEvents`. `MERGE=0` keeps one file per
shard instead, `KEEP_SHARDS=1` keeps the run directories `work_shards/`.

### Large productions

```bash
cd production
./runAMPT_fast.sh -l                                         # systems: pp_5360GeV, pp_13600GeV, pp_13000GeV
./runAMPT_fast.sh -s pp_5360GeV -n 8 -m 10 -e 100000         # 8 workers x 10 files x 100k events
./runAMPT_fast.sh -h                                         # all options
```

Writes `rootfiles/pp/<energy>/AMPT_<system>_<events>_MB<i>.root` (next to
`production/`) and logs in `Simulation/<tag>/`. Files already there are skipped,
so rerunning the same command fills in missing files after a failure or an
interruption. For a long production started in the background:

```bash
setsid nohup ./runAMPT_fast.sh -s pp_5360GeV -n 8 -m 10 > /dev/null 2>&1 < /dev/null &
tail -f ../Simulation/AMPT_pp_5360GeV_100k/production.log
```

### Seeds

With `ihjsed = 11` in `input.ampt`, AMPT reads the seeds from standard input:
line 1 = HIJING seed (0..1073741823; AMPT uses 2·seed+1), optional line 2 = ZPC
(parton cascade) seed (> 0; without it the value in `input.ampt` is used). The
scripts draw random seeds; the seeds are stored in the output file (`source`).

### Environment variables

| Variable | Effect |
| :--- | :--- |
| `AMPT_ROOT_OUT` | output file (default `ana/ampt.root`) |
| `AMPT_ENABLE_TXT=1` | also write `ana/ampt.dat` (the official text output) |
| `AMPT_ENABLE_ZPC=1` | also write `ana/zpc.dat` (final partons) |

---

## Output tree `t`

One entry per event; events without final-state particles are not written.
All branches are `Float_t` (F) or `Int_t` (I).

| Branch | Type | Content |
| :--- | :--- | :--- |
| `gImp` | F | impact parameter b (fm) |
| `Npart` | I | participant nucleons (AMPT npart1 + npart2) |
| `Ncoll` | I | binary nucleon-nucleon collisions (HIJING NCOLT) |
| `psiRP` | F | reaction-plane angle |
| `Nab` | I | number of nucleons, A_proj + A_targ (2 for pp) |
| `Stat[Nab]` | I | HIJING status of each nucleon, 0 = no interaction |
| `Nx`, `Ny`, `Nz` `[Nab]` | F | nucleon positions (fm), projectile first, shifted by ±b/2 along psiRP |
| `ntrack` | I | number of stored particles |
| `PID[ntrack]` | I | PDG code |
| `charge[ntrack]` | I | charge (e) |
| `px`, `py`, `pz` `[ntrack]` | F | momentum (GeV/c) |
| `energy[ntrack]` | F | energy (GeV) |

Stored particles: every charged particle (charge != 0, e and mu included) and
K0S, φ, Λ, Ξ, Ω (and antiparticles), at any η and pT.

Metadata (TNamed): `branches` (the table above), `source` (AMPT version and
the seeds), `tEvents` (events written, and not written because they had no
particles). If AMPT stops on an internal error, the file is still closed with
the events written so far.

Reading the tree, e.g. in a ROOT macro:

```cpp
TFile f("ana/ampt.root");
TTree *t = (TTree *)f.Get("t");
t->Draw("ntrack");                                  // multiplicity distribution
t->Draw("sqrt(px*px+py*py)", "abs(PID)==211");      // charged-pion pT
```

---

## Changes to the official code

No physics algorithm, parameter or random-number sequence is changed.

Check against the unmodified official code (built with `gfortran -O`), 1000 pp
events at 5.36 TeV with `input.ampt` and the same seeds, single core:

| | official code | this code (`make`) |
| :--- | ---: | ---: |
| time | 93.9 s | 19.1 s (4.9× faster) |
| `ana/ampt.dat` (121050 lines) | | identical except 3 particles, each off by 1 in the last printed digit of pz (10⁻⁴ GeV/c): floating-point rounding, the same events |

Small systems such as pp gain most, because the per-timestep work that was
removed dominates there. For more statistics, run events in parallel
(`run_fast_ampt.sh`, `production/runAMPT_fast.sh`): the events are independent,
so the throughput grows with the number of cores.

What was changed (each change keeps the result of the official code):

| Where | Official code | Now |
| :--- | :--- | :--- |
| `art1f.f` RELCOL (hadron cascade) | every hadron pair tested at every timestep | `nxtpar` skips the pairs rejected by the first tests (E=0, LB, distance > 5/10 fm) in a tight loop, same arithmetic |
| `art1f.f` RELCOL | 3 arrays of 150001 entries reset every timestep | reset up to the highest index used |
| `art1f.f` DENS | 11 grids of 41×41×49 cells cleared and scanned every timestep | only the box of cells holding hadrons |
| `linana.f` hbtout | O(N²) search of each hadron in the freeze-out record | hash on (px, ART code), same first match |
| `linana.f` locldr (coalescence) | 4 Lorentz boosts with 7 divisions each per quark pair | boost coefficients once per pair, unused boost skipped |
| `zpc.f` getict (parton cascade) | O(N) scan for the next collision at every collision | tournament tree, same parton (incl. ties) |
| `zpc.f` collision check chain | static (`SAVE`) locals | automatic locals (all set before use) |
| `zpc.f` isco1 (pair collision test) | 8 separate parton arrays (8 cache lines per parton) | packed copy, 1 cache line per parton, synchronized at every change |
| `main.f`, `linana.f` output | formatted `ana/ampt.dat`, `ana/zpc.dat` | ROOT tree `t`; text files optional |
| `main.f` seeds | HIJING seed at run time only | also the ZPC seed (line 2 of standard input) |
