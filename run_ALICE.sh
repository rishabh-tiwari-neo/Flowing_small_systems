# Hardcoded run configuration
TOTAL=10000

NCORES=4

g=10   #gAmplitude

GEN_BIN="generator/gen_ALICE"
DATA_DIR="data/data_ALICE"
LOG_DIR="logs/log_ALICE"  

echo "Building..."
rm -f "$GEN_BIN"
make "$GEN_BIN"
if [ $? -ne 0 ]; then
    echo "Build failed. Aborting."
    exit 1
fi

if [ ! -x "$GEN_BIN" ]; then
    echo "Expected binary $GEN_BIN not found after build. Aborting."
    exit 1
fi

mkdir -p "$DATA_DIR" "$LOG_DIR"
rm -f "$LOG_DIR"/progress_task*.txt "$LOG_DIR"/task*.log "$LOG_DIR"/summary.log
rm -f "$DATA_DIR"/output_task*.root

BASE=$((TOTAL / NCORES))
REM=$((TOTAL % NCORES))

SEEDS=()
SEED_LOG="$LOG_DIR/seeds.log"
echo "Seeds for this run ($(date)):" > "$SEED_LOG"
for ((i=0; i<NCORES; i++)); do
    S=$(( (RANDOM * RANDOM + RANDOM + i) % 900000000 + 1 ))
    SEEDS+=($S)
    echo "task $i -> seed $S" >> "$SEED_LOG"
done
echo "Seeds saved to $SEED_LOG"

echo "Launching $NCORES tasks for $TOTAL total events (~$BASE events/task)..."
echo "(each task: generate -> same-event pairing, per event)"

START_TIME=$(date +%s)

PIDS=()
TOTALS=()
for ((i=0; i<NCORES; i++)); do
    N=$BASE
    if [ $i -lt $REM ]; then
        N=$((N+1))
    fi
    TOTALS+=($N)
    SEED=${SEEDS[$i]}
    OUT="$DATA_DIR/output_task${i}.root"

    (
        TSTART=$(date +%s)
        echo "[task $i] started $(date), seed=$SEED, events=$N"
        ./$GEN_BIN $N $SEED $OUT $g
        TEND=$(date +%s)
        echo "[task $i] finished $(date), runtime=$((TEND - TSTART))s"
    ) > "$LOG_DIR/task${i}.log" 2>&1 &
    PIDS+=($!)
done

# --- multi-line live progress monitor ---
BAR_WIDTH=30

for ((i=0; i<NCORES; i++)); do echo ""; done

while true; do
    printf "\033[%dA" "$NCORES"

    ALIVE=0
    for ((i=0; i<NCORES; i++)); do
        SEED=${SEEDS[$i]}
        F="$LOG_DIR/progress_task${SEED}.txt"
        N=${TOTALS[$i]}

        V=0
        if [ -f "$F" ]; then
            V=$(cat "$F" 2>/dev/null)
            V=${V:-0}
        fi

        FRAC=$(awk -v v=$V -v n=$N 'BEGIN { if (n>0) printf "%.4f", v/n; else print 0 }')
        FILLED=$(awk -v f=$FRAC -v w=$BAR_WIDTH 'BEGIN { printf "%d", f*w }')
        PCT=$(awk -v f=$FRAC 'BEGIN { printf "%d", f*100 }')

        BAR=$(printf '%*s' "$FILLED" '' | tr ' ' '=')
        SPACES=$(printf '%*s' "$((BAR_WIDTH - FILLED))" '')

        printf "\rtask%-2d [gen ]: [%s%s] %3d%% (%d/%d)\033[K\n" "$i" "$BAR" "$SPACES" "$PCT" "$V" "$N"

        if kill -0 "${PIDS[$i]}" 2>/dev/null; then
            ALIVE=1
        fi
    done

    if [ $ALIVE -eq 0 ]; then
        break
    fi

    sleep 1
done

wait
END_TIME=$(date +%s)
TOTAL_RUNTIME=$((END_TIME - START_TIME))

echo "All tasks finished in ${TOTAL_RUNTIME}s. Merging..."

MERGED="$DATA_DIR/merged.root"
hadd -f "$MERGED" "$DATA_DIR"/output_task*.root
if [ $? -ne 0 ]; then
    echo "hadd failed. Aborting."
    exit 1
fi
echo "Merged file: $MERGED"

# store g in the merged file (read by plot_ALICE)
cat > "$LOG_DIR/add_g.C" <<MACRO
#include "TFile.h"
#include "TParameter.h"
void add_g(const char* fname, double g) {
    TFile f(fname, "UPDATE");
    TParameter<double> p("g", g);
    p.Write();
    f.Close();
}
MACRO
root -l -b -q "$LOG_DIR/add_g.C(\"$MERGED\", $g)"
if [ $? -ne 0 ]; then
    echo "Could not store g in $MERGED. Aborting."
    exit 1
fi

# --- write summary log ---
{
    echo "=== Run Summary ==="
    echo "Date: $(date)"
    echo "Total events requested: $TOTAL"
    echo "Cores used: $NCORES"
    echo "Total wall-clock runtime: ${TOTAL_RUNTIME}s"
    echo ""
    echo "--- Seeds used ---"
    cat "$SEED_LOG"
    echo ""
    echo "--- Per-task details ---"
    for ((i=0; i<NCORES; i++)); do
        echo "[task $i]"
        grep -E "started|finished|Done" "$LOG_DIR/task${i}.log"
        echo ""
    done
} > "$LOG_DIR/summary.log"

echo "Summary written to $LOG_DIR/summary.log"

# clean up progress files
rm -f "$LOG_DIR"/progress_task*.txt

# --- post-processing: template fit and plots ---
if [ ! -f "$MERGED" ]; then
    echo "Merged file $MERGED not found. Aborting."
    exit 1
fi

echo "Compiling post_ALICE and plot_ALICE..."
g++ -O2 -o post/post_ALICE post/post_ALICE.cc `root-config --cflags --libs`
if [ $? -ne 0 ]; then
    echo "post_ALICE build failed. Aborting."
    exit 1
fi
g++ -O2 -o post/plot_ALICE post/plot_ALICE.cc `root-config --cflags --libs`
if [ $? -ne 0 ]; then
    echo "plot_ALICE build failed. Aborting."
    exit 1
fi

echo ""
echo "=== Fitting ==="
./post/post_ALICE 2>&1 | tee "$LOG_DIR/post_all.txt"
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "post_ALICE failed. Aborting."
    exit 1
fi

echo "=== Plotting ==="
./post/plot_ALICE
if [ $? -ne 0 ]; then
    echo "plot_ALICE failed. Aborting."
    exit 1
fi

echo ""
echo "Done. Results in $DATA_DIR, plots in plots/plot_ALICE_HM and plots/plot_ALICE_HM35."