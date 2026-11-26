#!/bin/bash
# run_AMPT_hist_maker.sh
# Turns AMPT tree files into the histogram file that post/plot read, using NCORES processes.
#
# Put this script and Histogram_AMPT_par.cc in any folder (e.g. next to the trees), edit the
# settings below, and run:   ./run_AMPT_hist_maker.sh
#
# What it does:
#   1. builds Histogram_AMPT_par with the plain g++ line (no Makefile and no $PYTHIA8 needed)
#   2. lists the tree files (INPUT_DIR/INPUT_PATTERN, sorted) = one long chain of events
#   3. splits the events into NCORES equal, non-overlapping ranges, one process per range
#   4. hadds the block files of each core, then the cores -> $OUT_DIR/$MERGED_NAME
# It never deletes the input trees and never touches an existing merged.root.

# ================= settings =================
NCORES=1

INPUT_DIR="."                # folder holding the AMPT tree files
INPUT_PATTERN="mb1p5M.root" #Write the file name here
OUT_DIR="hist_out_1core"           # results go here (made if missing)
MERGED_NAME="merged.root"

V0M_HM_MIN=142
V0M_HM35_MIN=72
NCH_HM_MIN=25
NCH_LM_MAX=20
INPUT_AMPT=""                # path to the input.ampt used for the run (logged in generator_info.txt); may be empty
CHUNK=1000000                # events per block (= per mixed-event pool); 0 = one block per core

KEEP_PARTS=0                 # 1 = keep the per-block and per-core files after the final merge

SRC="Histogram_AMPT.cc"
BIN="./Histogram_AMPT"
# ============================================

cd "$(dirname "$0")" || exit 1

for tool in g++ root-config hadd; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Missing tool: $tool (is ROOT set up?)"; exit 1; }
done
[ -f "$SRC" ] || { echo "Cannot find $SRC next to this script."; exit 1; }

echo "Building $BIN ..."
g++ -O2 -o "$BIN" "$SRC" $(root-config --cflags --libs) || { echo "Build failed. Aborting."; exit 1; }

mkdir -p "$OUT_DIR"
if ls "$OUT_DIR"/h_core*_part*.root >/dev/null 2>&1; then
    echo "$OUT_DIR already holds h_core*_part*.root from an earlier run."
    echo "Move or remove them first (this script does not delete old results). Aborting."
    exit 1
fi

LIST="$OUT_DIR/filelist.txt"
ls -1 "$INPUT_DIR"/$INPUT_PATTERN 2>/dev/null | sort -V > "$LIST"
NFILES=$(wc -l < "$LIST")
if [ "$NFILES" -eq 0 ]; then
    echo "No files match $INPUT_DIR/$INPUT_PATTERN. Edit INPUT_DIR / INPUT_PATTERN. Aborting."
    exit 1
fi
echo "$NFILES tree file(s) listed in $LIST (first: $(head -n 1 "$LIST"))"

echo "Launching $NCORES process(es) ..."
START=$(date +%s)
PIDS=()
for ((i=0; i<NCORES; i++)); do
    $BIN $i $NCORES "$LIST" "$OUT_DIR/h" $V0M_HM_MIN $V0M_HM35_MIN $NCH_HM_MIN $NCH_LM_MAX "$INPUT_AMPT" $CHUNK \
        > "$OUT_DIR/core$i.log" 2>&1 &
    PIDS+=($!)
done

# simple progress: how many block files exist
while true; do
    ALIVE=0
    for p in "${PIDS[@]}"; do kill -0 "$p" 2>/dev/null && ALIVE=1; done
    NB=$(ls "$OUT_DIR"/h_core*_part*.root 2>/dev/null | wc -l)
    printf "\r[%4ds] %d block file(s) written" $(( $(date +%s) - START )) "$NB"
    [ $ALIVE -eq 0 ] && break
    sleep 10
done
echo ""

FAILED=0
for ((i=0; i<NCORES; i++)); do
    wait "${PIDS[$i]}" || { echo "core $i FAILED, see $OUT_DIR/core$i.log"; FAILED=1; }
done
[ $FAILED -eq 0 ] || { echo "Not merging because a core failed. Aborting."; exit 1; }
echo "All cores finished in $(( $(date +%s) - START ))s. Merging ..."

for ((i=0; i<NCORES; i++)); do
    hadd -f "$OUT_DIR/core_$i.root" "$OUT_DIR"/h_core${i}_part*.root > "$OUT_DIR/hadd_core$i.log" 2>&1 \
        || { echo "hadd of core $i failed, see $OUT_DIR/hadd_core$i.log. Aborting."; exit 1; }
done
hadd -f "$OUT_DIR/$MERGED_NAME" "$OUT_DIR"/core_*.root > "$OUT_DIR/hadd_final.log" 2>&1 \
    || { echo "Final hadd failed, see $OUT_DIR/hadd_final.log. Aborting."; exit 1; }

if [ "$KEEP_PARTS" -eq 0 ]; then
    rm -f "$OUT_DIR"/h_core*_part*.root "$OUT_DIR"/core_*.root
fi

echo "Done: $OUT_DIR/$MERGED_NAME  and  $OUT_DIR/generator_info.txt"
echo "To use with post/plot, copy them by hand into data/data_AMPT/ (this script does not overwrite merged.root)."
