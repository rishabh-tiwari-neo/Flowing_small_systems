#!/usr/bin/env bash
set -euo pipefail

# Multi-core AMPT runner: the events are split into shards of EVENTS_PER_SHARD
# events, which N_WORKERS AMPT processes work through (a new shard starts as
# soon as one finishes, so slow central events do not leave cores idle at the
# end); the shard ROOT files are merged with hadd.
#
# Usage: ./run_fast_ampt.sh [N_WORKERS] [TOTAL_EVENTS] [OUTPUT_ROOT_FILE] [EVENTS_PER_SHARD]
#   N_WORKERS         default: nproc
#   TOTAL_EVENTS      default: NEVNT of input.ampt
#   OUTPUT_ROOT_FILE  default: ana/ampt_merged.root
#   EVENTS_PER_SHARD  default: TOTAL_EVENTS / (4 * N_WORKERS), at least 1
# Environment: AMPT_ENABLE_TXT=1 / AMPT_ENABLE_ZPC=1 also keep ana/ampt.dat /
# ana/zpc.dat of each shard (in work_shards/, kept with KEEP_SHARDS=1);
# MERGE=0 keeps one file per shard, OUTPUT_ROOT_FILE without .root + _MB<i>.root.

NJOBS=${1:-$(nproc)}
CONFIG_EVENTS=$(awk '/! *NEVNT/ {print $1; exit}' input.ampt 2>/dev/null || echo "10")
TOTAL_EVENTS=${2:-${CONFIG_EVENTS}}
FINAL_OUT=${3:-"ana/ampt_merged.root"}
PER_SHARD=${4:-$(( TOTAL_EVENTS / (4 * NJOBS) ))}
[[ "${PER_SHARD}" -lt 1 ]] && PER_SHARD=1
NSHARDS=$(( (TOTAL_EVENTS + PER_SHARD - 1) / PER_SHARD ))
[[ "${NJOBS}" -gt "${NSHARDS}" ]] && NJOBS="${NSHARDS}"

if [[ ! -x "./ampt" ]]; then
    echo "Executable './ampt' not found. Running make ..."
    make -j"$(nproc)"
fi

echo "=========================================================="
echo " Starting Parallel AMPT Generation"
echo " Workers         : ${NJOBS}"
echo " Total Events    : ${TOTAL_EVENTS} (${NSHARDS} shards of <= ${PER_SHARD})"
echo " Output Target   : ${FINAL_OUT}"
echo "=========================================================="

mkdir -p ana work_shards
rm -rf work_shards/shard_*

# Random seeds: HIJING seed in 1..10^9 (AMPT uses 2*NSEED+1, which must stay
# below 2^31) and a positive ZPC seed (zpc.f ran1 needs iseedp > 0); the
# HIJING seeds of the shards are all different
declare -A USED_SEEDS=()
random_seed() {
    local r
    if [[ -r /dev/urandom ]]; then
        r=$(od -vAn -N4 -tu4 < /dev/urandom | tr -d ' ')
    else
        r=$(( RANDOM * 32768 + RANDOM ))
    fi
    echo $(( r % 1000000000 + 1 ))
}

PIDS=()
cleanup() {
    if [[ ${#PIDS[@]} -gt 0 ]]; then
        echo ""
        echo "Terminating running AMPT workers..."
        for pid in "${PIDS[@]}"; do
            kill -TERM "${pid}" 2>/dev/null || true
        done
        wait 2>/dev/null || true
    fi
    exit 1
}
trap cleanup SIGINT SIGTERM

run_shard() {   # shard_id n_events hijing_seed zpc_seed
    local dir="work_shards/shard_$1"
    mkdir -p "${dir}/ana"
    sed -e "s/^[0-9]*\(\s*! NEVNT\)/$2\1/" \
        -e 's/^[0-9]*\(\s*! ihjsed\)/11\1/' \
        input.ampt > "${dir}/input.ampt"
    cd "${dir}"
    AMPT_ROOT_OUT="ana/shard_$1.root" ../../ampt > ampt.log 2>&1 <<EOF
$3
$4
EOF
}

START_TIME=$(date +%s)
FAILURES=0
REMAINING=${TOTAL_EVENTS}
for id in $(seq 1 "${NSHARDS}"); do
    NEV=$(( REMAINING < PER_SHARD ? REMAINING : PER_SHARD ))
    REMAINING=$(( REMAINING - NEV ))
    SEED_HIJING=$(random_seed)
    while [[ -n "${USED_SEEDS[${SEED_HIJING}]:-}" ]]; do SEED_HIJING=$(random_seed); done
    USED_SEEDS[${SEED_HIJING}]=1
    SEED_ZPC=$(random_seed)

    # at most NJOBS shards at a time
    while [[ ${#PIDS[@]} -ge ${NJOBS} ]]; do
        if ! wait -n; then FAILURES=$(( FAILURES + 1 )); fi
        ALIVE=()
        for pid in "${PIDS[@]}"; do
            kill -0 "${pid}" 2>/dev/null && ALIVE+=("${pid}")
        done
        PIDS=("${ALIVE[@]+"${ALIVE[@]}"}")
    done
    echo "  [shard ${id}/${NSHARDS}] ${NEV} events (HIJING seed ${SEED_HIJING}, ZPC seed ${SEED_ZPC})"
    ( run_shard "${id}" "${NEV}" "${SEED_HIJING}" "${SEED_ZPC}" ) &
    PIDS+=("$!")
done

for pid in "${PIDS[@]+"${PIDS[@]}"}"; do
    if ! wait "${pid}"; then FAILURES=$(( FAILURES + 1 )); fi
done
PIDS=()

if [[ "${FAILURES}" -gt 0 ]]; then
    echo "ERROR: ${FAILURES} AMPT shards exited with an error. Check work_shards/shard_*/ampt.log" >&2
    exit 1
fi

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
[[ "${ELAPSED}" -eq 0 ]] && ELAPSED=1

echo ""
mkdir -p "$(dirname "${FINAL_OUT}")"
if [[ "${MERGE:-1}" == "0" ]]; then
    # one file per shard
    for f in work_shards/shard_*/ana/shard_*.root; do
        id=${f##*/shard_}
        mv -n "${f}" "${FINAL_OUT%.root}_MB${id%.root}.root"
    done
    echo "Shard files: ${FINAL_OUT%.root}_MB<1..${NSHARDS}>.root"
else
echo "Merging ${NSHARDS} ROOT shards into ${FINAL_OUT} ..."
if ! hadd -f505 -j "${NJOBS}" "${FINAL_OUT}" work_shards/shard_*/ana/shard_*.root > work_shards/hadd.log 2>&1; then
    cat work_shards/hadd.log >&2
    exit 1
fi

# hadd keeps the metadata (TNamed) of every shard as extra cycles: keep one
# of each, one "source" listing the seeds of all shards, and the event counts
# "tEvents" summed over the shards
cat > work_shards/merge_meta.C <<'EOF'
void merge_meta(const char *fname) {
  TFile f(fname, "UPDATE");
  std::map<TString, TString> title;
  TString sources;
  long long nWritten = 0, nEmpty = 0;
  std::vector<TString> cycles;
  for (TObject *o : *f.GetListOfKeys()) {
    TKey *k = (TKey *)o;
    if (TString(k->GetClassName()) != "TNamed") continue;
    TNamed *n = (TNamed *)k->ReadObj();
    if (TString(k->GetName()) == "source")
      sources += (sources.IsNull() ? "" : "; ") + TString(n->GetTitle());
    else if (TString(k->GetName()) == "tEvents") {
      long long a = 0, b = 0;
      sscanf(n->GetTitle(), "written %lld, not written: no particles %lld", &a, &b);
      nWritten += a; nEmpty += b;
    } else
      title[k->GetName()] = n->GetTitle();
    cycles.push_back(Form("%s;%d", k->GetName(), k->GetCycle()));
  }
  for (const TString &c : cycles) f.Delete(c);
  for (auto &t : title) TNamed(t.first, t.second).Write();
  TNamed("source", sources.Data()).Write();
  TNamed("tEvents", Form("written %lld, not written: no particles %lld", nWritten, nEmpty)).Write();
  f.Close();
}
EOF
root -l -b -q "work_shards/merge_meta.C(\"${FINAL_OUT}\")" > work_shards/merge_meta.log 2>&1 \
    || echo "WARNING: metadata not merged, see work_shards/merge_meta.log" >&2
fi

[[ "${KEEP_SHARDS:-0}" == "1" ]] || rm -rf work_shards

RATE=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_EVENTS} / ${ELAPSED}}")
echo "=========================================================="
echo " Finished successfully in ${ELAPSED} seconds (${RATE} events/sec)!"
[[ "${MERGE:-1}" == "0" ]] && echo " Output ROOT Files: ${FINAL_OUT%.root}_MB*.root" || echo " Output ROOT File: ${FINAL_OUT}"
echo "=========================================================="
