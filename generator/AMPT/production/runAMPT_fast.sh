#!/bin/bash
# Mass production of pp events with the fast-root AMPT: NJOBS workers, one core each, produce
# NJOBS*NCYCLES files of NEVENTS events:
#   1. run AMPT with its own random HIJING and ZPC seeds; AMPT writes the ROOT tree t itself
#      (no ampt.dat / zpc.dat / conversion step),
#   2. check the event count and move the ROOT file to OUTDIR,
#   3. delete the run directory, keeping a short log.
# Output files are <TAG>_MB<i>.root with i = 1 .. NJOBS*NCYCLES; files already in OUTDIR are
# skipped, so rerunning the same command fills in only the missing ones, and rerunning with a
# larger NCYCLES extends the production. The missing indices form a queue that the workers take
# from in increasing order. A worker stops at its first failure and keeps that run directory;
# the others go on. A second production of the same OUTDIR and TAG refuses to start while the
# first runs (lock file OUTDIR/.TAG.lock).
#
# Usage:
#   runAMPT_fast.sh -s SYSTEM [-n NJOBS] [-m NCYCLES] [-e NEVENTS] [-t TAG] [-b BUILD]
#                   [-o OUTDIR] [-w WORKDIR] [-r RUNBASE] [-i INPUT_AMPT]
#   runAMPT_fast.sh -l      list the systems with their input and OUTDIR
#   SYSTEM          pp_5360GeV, pp_13600GeV or pp_13000GeV: sets INPUT_AMPT (from examples/)
#                   and OUTDIR (rootfiles/pp/<energy> next to production/); both can be
#                   overridden by their option
#   NJOBS           parallel workers (default: number of cores, nproc)
#   NCYCLES         files per worker (default 1)
#   NEVENTS         events per file (default 100000)
#   TAG             default AMPT_<SYSTEM>_<NEVENTS>
#   BUILD           directory with the ampt executable (default: the directory above
#                   production/, after make)
#   RUNBASE         local disk for the run directories (default /var/tmp/$USER/ampt_fast)
#   WORKDIR         logs (default: <base>/Simulation/<TAG> for OUTDIR = <base>/rootfiles/...,
#                   else OUTDIR/Simulation/<TAG>)
# Detached long production:
#   setsid nohup ./runAMPT_fast.sh ... > /dev/null 2>&1 < /dev/null &
# Progress:  tail -f <WORKDIR>/production.log
# Stop all:  kill -- -$(cat <WORKDIR>/production_fast.pgid)

TOP=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PPW=${AMPT_PP_OUTDIR:-$TOP/rootfiles/pp}

# system: input.ampt in examples/, output directory
SYSTEMS="
pp_5360GeV   input.pp_5360GeV.ampt    $PPW/5360GeV
pp_13600GeV  input.pp_13600GeV.ampt   $PPW/13600GeV
pp_13000GeV  input.pp_13000GeV.ampt   $PPW/13000GeV
"

NJOBS=$(nproc)
NCYCLES=1
NEVENTS=100000
SYSTEM=
TAG=
OUTDIR=
BUILD=
INPUT=
RUNBASE=/var/tmp/$USER/ampt_fast
WORK=

usage() { sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts "s:n:m:e:t:o:b:i:r:w:lh" opt; do
    case $opt in
        s) SYSTEM=$OPTARG ;;
        n) NJOBS=$OPTARG ;;
        m) NCYCLES=$OPTARG ;;
        e) NEVENTS=$OPTARG ;;
        t) TAG=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        b) BUILD=$(realpath "$OPTARG") ;;
        i) INPUT=$(realpath "$OPTARG") ;;
        r) RUNBASE=$OPTARG ;;
        w) WORK=$OPTARG ;;
        l) echo "$SYSTEMS" | sed '/^$/d'; exit 0 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
[[ -n $SYSTEM ]] || { echo "need -s SYSTEM"; usage 1; }
read -r _ sinput soutdir < <(echo "$SYSTEMS" | awk -v s="$SYSTEM" '$1 == s')
[[ -n $sinput ]] || { echo "unknown system $SYSTEM (runAMPT_fast.sh -l lists them)"; exit 1; }
INPUT=${INPUT:-$TOP/examples/$sinput}
OUTDIR=${OUTDIR:-$soutdir}
BUILD=${BUILD:-$TOP}
if [[ -z $TAG ]]; then
    if ((NEVENTS % 1000000 == 0)); then nev=$((NEVENTS / 1000000))M
    elif ((NEVENTS % 1000 == 0)); then nev=$((NEVENTS / 1000))k
    else nev=$NEVENTS; fi
    TAG=AMPT_${SYSTEM}_$nev
fi
# logs next to the rootfiles directory: <base>/rootfiles/... -> <base>/Simulation/<TAG>
WORK=${WORK:-${OUTDIR%%/rootfiles*}/Simulation/$TAG}
LOGDIR=$WORK/logs
mkdir -p "$LOGDIR" || exit 1

log() { echo "$(date '+%F %T') $*" | tee -a "$WORK/production.log"; }
die() { log "ERROR: $*"; exit 1; }

set -uo pipefail

[[ $NJOBS =~ ^[1-9][0-9]*$ && $NCYCLES =~ ^[1-9][0-9]*$ && $NEVENTS =~ ^[1-9][0-9]*$ ]] \
    || die "NJOBS, NCYCLES and NEVENTS must be positive integers"
[[ -x $BUILD/ampt ]] || die "no ampt executable in $BUILD (run make first)"
[[ -s $INPUT ]] || die "input.ampt template not found: $INPUT"
mkdir -p "$RUNBASE/$TAG" || die "cannot create $RUNBASE/$TAG"
mkdir -p "$OUTDIR" && [[ -w $OUTDIR ]] || die "output directory not writable: $OUTDIR"
# One production per OUTDIR and TAG: a second one would delete the run
# directories of the first. The workers and AMPT inherit fd 8, so the lock is held until they end.
LOCK=$OUTDIR/.$TAG.lock
exec 8>> "$LOCK" || die "cannot open the lock file $LOCK"
flock -n 8 || die "another production of $TAG is running ($(cat "$LOCK")), nothing started"
echo "$(hostname), pid $$, since $(date '+%F %T')" > "$LOCK"
ps -o pgid= -p $$ | tr -d ' ' > "$WORK/production_fast.pgid"

log "=== fast-root start on $(hostname): $SYSTEM, $NJOBS jobs x $NCYCLES cycles x $NEVENTS events"
log "ampt $BUILD/ampt ($(md5sum < "$BUILD/ampt" | cut -c1-32)), run dirs $RUNBASE/$TAG, output $OUTDIR/${TAG}_MB<i>.root"
log "input $INPUT ($(md5sum < "$INPUT" | cut -c1-32))"

# input.ampt: the template with NEVNT set and the seeds read at run time
sed -E -e "s/^[0-9]+([[:space:]]+! NEVNT)/${NEVENTS}\1/" -e "s/^[0-9]+([[:space:]]+! ihjsed)/11\1/" \
    "$INPUT" > "$WORK/input_fast.ampt"
grep -qE "^${NEVENTS}[[:space:]]+! NEVNT" "$WORK/input_fast.ampt" || die "could not set NEVNT in input.ampt"
grep -qE "^11[[:space:]]+! ihjsed" "$WORK/input_fast.ampt" || die "could not set ihjsed in input.ampt"

# Distinct HIJING seeds for all runs (AMPT uses 2*seed+1 as a 32-bit integer: below 1e9) and
# positive ZPC seeds
mapfile -t SEEDS < <(shuf -i 1-999999999 -n $((NJOBS * NCYCLES)) --random-source=/dev/urandom)
mapfile -t ZSEEDS < <(shuf -i 1-999999999 -n $((NJOBS * NCYCLES)) --random-source=/dev/urandom)

# Checks the summary line of tree t in the AMPT output $1: events written + empty must be
# NEVENTS, at least one written. Prints the numbers.
check_tree() {
    local out=$1 summary nstored nempty
    summary=$(grep -a '^\[AmptRootIO\] tree t events:' "$out" | tail -1)
    nstored=$(sed -n 's/.*events: \([0-9]*\).*/\1/p' <<< "$summary")
    nempty=$(sed -n 's/.*no particles: \([0-9]*\).*/\1/p' <<< "$summary")
    if [[ -z $nstored || -z $nempty ]] || ((nstored == 0 || nstored + nempty != NEVENTS)); then
        echo "'$summary'"
        return 1
    fi
    echo "$nstored events in tree t, $nempty empty not written"
}

run_one() {
    local idx=$1 seed=$2 zseed=$3
    local name=${TAG}_MB${idx}
    local out=$OUTDIR/$name.root
    local rundir=$RUNBASE/$TAG/run_MB$idx

    if [[ -e $out ]]; then
        log "MB$idx: $out already exists, skipping"
        return 0
    fi
    rm -rf "$rundir"
    mkdir -p "$rundir/ana"
    cp "$BUILD/ampt" "$rundir/"
    cp "$WORK/input_fast.ampt" "$rundir/input.ampt"
    echo "MB$idx $seed $zseed" >> "$WORK/seeds_fast.txt"

    log "MB$idx: AMPT started (seeds $seed $zseed)"
    (cd "$rundir" && AMPT_ROOT_OUT="$rundir/$name.root" ./ampt > nohup.out 2>&1 <<< "$seed
$zseed")
    local rc=$?
    tail -n 30 "$rundir/nohup.out" > "$LOGDIR/${name}_ampt.log"
    local result
    if ! result=$(check_tree "$rundir/nohup.out") || ((rc != 0)) || [[ ! -s $rundir/$name.root ]]; then
        log "MB$idx: AMPT failed (exit code $rc, summary $result), run directory kept: $rundir"
        return 1
    fi
    mv -n "$rundir/$name.root" "$out"
    if [[ -e $rundir/$name.root ]]; then
        log "MB$idx: could not move the output to $out, run directory kept: $rundir"
        return 1
    fi
    rm -rf "$rundir"
    log "MB$idx: done, $out ($(du -h "$out" | cut -f1)): $result"
}

# Work queue on the local disk: the indices without an output file and a read position, shared
# by the workers under a lock
QUEUE=$RUNBASE/$TAG/queue_fast.txt
for ((i = 1; i <= NJOBS * NCYCLES; i++)); do
    [[ -e $OUTDIR/${TAG}_MB$i.root ]] || echo "$i"
done > "$QUEUE"
echo 0 > "$QUEUE.pos"
log "$(wc -l < "$QUEUE") of $((NJOBS * NCYCLES)) files to produce"

next_index() {
    (
        flock 9 || exit 1
        local n
        n=$(< "$QUEUE.pos")
        echo $((n + 1)) > "$QUEUE.pos"
        sed -n "$((n + 1))p" "$QUEUE"
    ) 9> "$QUEUE.lock"
}

worker() {
    local j=$1 idx
    while idx=$(next_index) && [[ -n $idx ]]; do
        if ! run_one "$idx" "${SEEDS[idx - 1]}" "${ZSEEDS[idx - 1]}"; then
            log "worker $j: stopping after failure in MB$idx"
            return 1
        fi
    done
}

pids=()
for ((j = 1; j <= NJOBS; j++)); do
    worker "$j" &
    pids+=($!)
done

nfail=0
for pid in "${pids[@]}"; do
    wait "$pid" || ((nfail++))
done

nout=$(ls "$OUTDIR/${TAG}"_MB*.root 2> /dev/null | wc -l)
if ((nfail > 0)); then
    log "=== finished with $nfail failed worker(s); $nout ${TAG} files in $OUTDIR. Rerun the same command to fill the gaps."
    exit 1
fi
log "=== all done: $nout ${TAG} files in $OUTDIR"
