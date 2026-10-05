#!/bin/bash
# Hardcoded run configuration
TOTAL=10000
NCORES=5

GEN_BIN="generator/stat_check"
SRC="generator/stat_check.cc"
OUT_DIR="data/test_run"
LOG_DIR="logs/test_run"

# --- clear old binary and recompile ---
echo "Removing old binary..."
rm -f "$GEN_BIN"

echo "Compiling $GEN_BIN..."
make "$GEN_BIN"
if [ $? -ne 0 ] || [ ! -x "$GEN_BIN" ]; then
    echo "Build failed. Aborting."
    exit 1
fi
echo "Build OK."

mkdir -p "$OUT_DIR" "$LOG_DIR" logs
rm -f "$OUT_DIR"/output_task*.root "$OUT_DIR"/merged.root
rm -f "$LOG_DIR"/task*.log "$LOG_DIR"/summary.log

BASE=$((TOTAL / NCORES))
REM=$((TOTAL % NCORES))

SEEDS=()
SEED_LOG="$LOG_DIR/seeds.log"
echo "Seeds for this run ($(date)):" > $SEED_LOG
for ((i=0; i<NCORES; i++)); do
    S=$(( (RANDOM * RANDOM + RANDOM + i) % 900000000 + 1 ))
    SEEDS+=($S)
    echo "task $i -> seed $S" >> $SEED_LOG
done
echo "Seeds saved to $SEED_LOG"

echo "Launching $NCORES tasks for $TOTAL total events (~$BASE events/task)..."

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
    OUT="$OUT_DIR/output_task${i}.root"

    (
        TSTART=$(date +%s)
        echo "[task $i] started $(date), seed=$SEED, events=$N"
        ./$GEN_BIN $N $SEED $OUT
        TEND=$(date +%s)
        echo "[task $i] finished $(date), runtime=$((TEND - TSTART))s"
    ) > $LOG_DIR/task${i}.log 2>&1 &
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
        F="logs/progress_task${SEED}.txt"
        N=${TOTALS[$i]}

        V=0
        if [ -f "$F" ]; then
            V=$(cat "$F" 2>/dev/null)
            V=${V:-0}
        fi

        FRAC=$(awk -v v=$V -v n=$N 'BEGIN { if (n>0) printf "%.4f", v/n; else print 0 }')
        FILLED=$(awk -v f=$FRAC -v w=$BAR_WIDTH 'BEGIN { printf "%d", f*w }')
        PCT=$(awk -v f=$FRAC 'BEGIN { printf "%d", f*100 }')

        BAR=$(printf "%0.s=" $(seq 1 $FILLED 2>/dev/null))
        SPACES=$(printf "%0.s " $(seq 1 $((BAR_WIDTH - FILLED)) 2>/dev/null))

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

hadd -f $OUT_DIR/merged.root $OUT_DIR/output_task*.root
echo "Merged file: $OUT_DIR/merged.root"

# --- write summary log ---
{
    echo "=== Run Summary ==="
    echo "Date: $(date)"
    echo "Binary: $GEN_BIN"
    echo "Total events requested: $TOTAL"
    echo "Cores used: $NCORES"
    echo "Total wall-clock runtime: ${TOTAL_RUNTIME}s"
    echo ""
    echo "--- Seeds used ---"
    cat $SEED_LOG
    echo ""
    echo "--- Per-task details ---"
    for ((i=0; i<NCORES; i++)); do
        echo "[task $i]"
        grep -E "started|finished|Done" $LOG_DIR/task${i}.log
        echo ""
    done
} > $LOG_DIR/summary.log

echo "Summary written to $LOG_DIR/summary.log"

# clean up this run's progress files
for S in "${SEEDS[@]}"; do
    rm -f "logs/progress_task${S}.txt"
done

echo "Done."