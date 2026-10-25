#!/bin/bash
# Every experiment of the paper's evaluation, all seven datasets, all systems,
# on a quiet machine: each step waits until the machine is idle, records the
# load while it runs (other users' CPU sampled every 10 s), and is redone once
# if someone else used more than DISTURB % CPU meanwhile.
# Output: experiments/results/paper_clean/<dataset>/ (one CSV per step) and
# experiments/results/paper_clean/load.csv; then paper_comparison_clean.md.
# Started once by an @reboot crontab entry, which it removes.
set -u
( crontab -l 2>/dev/null | grep -v run_paper_clean ) | crontab - 2>/dev/null
G=/home/cxl6029/projects/GIR; B=$G/experiments/baselines; P=$G/data/protocol
OUT=$G/experiments/results/paper_clean; mkdir -p $OUT
LOAD=$OUT/load.csv
[ -f $LOAD ] || echo "dataset,step,attempt,start,end,seconds,load1_start,load1_end,max_other_cpu_pct,disturbed" > $LOAD
DISTURB=${DISTURB:-300}   # other users' CPU, in % of one core
ME=$(id -un)
other_cpu() { ps -eo user=,pcpu= | awk -v me=$ME '$1 != me && $1 != "root" {s += $2} END {printf "%d", s}'; }
load1() { cut -d' ' -f1 /proc/loadavg; }
log() { echo "[$(date -u +%F' '%T)] $*"; }

wait_quiet() {  # 3 checks in a row 20 s apart: load < 4 and others < DISTURB/2 (gives up waiting after 6 h)
    local ok=0 waited=0
    while [ $ok -lt 3 ]; do
        l=$(load1); o=$(other_cpu)
        if awk -v l=$l -v o=$o -v d=$DISTURB 'BEGIN { exit !(l < 4 && o < d / 2) }'; then ok=$((ok + 1)); else ok=0; fi
        [ $ok -lt 3 ] && sleep 20 && waited=$((waited + 20))
        if [ $waited -ge 21600 ]; then log "machine still busy after 6 h (load $l, others ${o}%): going on"; return; fi
    done
}

# step <dataset> <name> <log file> <command...>: run gated, sampled, redone once if disturbed.
step() {
    local D=$1 name=$2 out=$3; shift 3
    for attempt in 1 2; do
        wait_quiet
        local s=$(date +%s) l0=$(load1) samp=$OUT/.samp.$$
        ( while true; do other_cpu; echo; sleep 10; done ) > $samp 2>/dev/null &
        local sp=$!
        "$@" > "$out" 2>&1
        local rc=$?
        kill $sp 2>/dev/null; wait $sp 2>/dev/null
        local e=$(date +%s) l1=$(load1) mx=$(sort -n $samp 2>/dev/null | tail -1); rm -f $samp
        mx=${mx:-0}
        local dist=0; [ "$mx" -gt "$DISTURB" ] && dist=1
        echo "$D,$name,$attempt,$(date -u -d @$s +%FT%T),$(date -u -d @$e +%FT%T),$((e - s)),$l0,$l1,$mx,$dist" >> $LOAD
        log "$D $name: rc $rc, $((e - s)) s, others max ${mx}%$([ $dist = 1 ] && echo ' (disturbed)')"
        [ $dist = 0 ] && return $rc
        [ $attempt = 1 ] && log "$D $name: redoing (disturbed)"
    done
    return $rc
}

# failure rows for a baseline that crashed or timed out
fail_rows() {  # <csv> <dataset> <system> <reason> <metrics...>
    local csv=$1 ds=$2 sys=$3 why=$4; shift 4
    for m in "$@"; do echo "$ds,$sys,$m,FAILED:$why,-1,-,0,64" >> $csv; done
}

run_dataset() {  # <dataset> <sweep> <big: 0/1>
    local D=$1 SWEEP=$2 BIG=$3 MAXT=64 T=$([ $3 = 1 ] && echo 14400 || echo 3600)
    local O=$OUT/$D; mkdir -p $O
    log "== $D (sweep $SWEEP)"
    cd $G/gir_constructor
    step $D gir $O/gir.log env OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 \
        --threads $SWEEP --sources 16 --bc-sources 64 --csv $O/gir.csv
    step $D gir_bod $O/gir_bod.log env OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 \
        --skip-queries --blocks-on-demand --csv $O/gir_bod.csv
    if [ $BIG = 0 ]; then
        step $D gir_sym $O/gir_sym.log env OMP_NUM_THREADS=$MAXT ./bench_gir --datasets ${D}_sym --ops 50000 --batch 10000 --budget 300 \
            --sources 16 --bc-sources 64 --csv $O/gir_sym.csv
    fi
    for m in fmu lazy inc lazy-inc; do
        rm -f $O/b_gir_stream_$m.csv
        step $D gir_stream_$m $O/gir_stream_$m.log env OMP_NUM_THREADS=$MAXT $B/gir/gir_stream $P/$D $D $O/b_gir_stream_$m.csv 10000 10 \
            $([ $m = fmu ] || echo $m)
    done
    cd $B
    rm -f $O/b_terrace.csv
    step $D terrace $O/terrace.log env CILK_NWORKERS=$MAXT timeout $T $B/terrace/terrace_driver $P/$D $D $O/b_terrace.csv 10000 \
        || fail_rows $O/b_terrace.csv $D Terrace "rc$?" insert_MEPS delete_MEPS mixed_MEPS bfs_ms
    rm -f $O/b_terrace_sym.csv
    step $D terrace_sym $O/terrace_sym.log env CILK_NWORKERS=$MAXT timeout $T $B/terrace/terrace_driver $P/${D}_sym ${D}_sym $O/b_terrace_sym.csv 10000 \
        || fail_rows $O/b_terrace_sym.csv ${D}_sym Terrace "rc$?" insert_MEPS delete_MEPS mixed_MEPS bfs_ms bc_ms_total
    rm -f $O/b_aspen.csv
    step $D aspen $O/aspen_sym.log env NUM_THREADS=$MAXT timeout $T $B/aspen/aspen_driver $P/${D}_sym ${D}_sym $O/b_aspen.csv 10000 \
        || fail_rows $O/b_aspen.csv ${D}_sym Aspen "rc$?" insert_MEPS delete_MEPS mixed_MEPS bfs_ms bc_ms_total
    # Spruce's mixed stream hangs (insert at a vertex whose last edge was deleted; it hung on
    # every graph of the first full run), so it is recorded as such rather than waited for.
    rm -f $O/b_spruce.csv
    step $D spruce $O/spruce.log env OMP_NUM_THREADS=$MAXT SPRUCE_SKIP_MIXED=1 timeout $T $B/spruce/spruce_driver $P/$D $D $O/b_spruce.csv 10000 \
        || fail_rows $O/b_spruce.csv $D Spruce "rc$?" insert_MEPS delete_MEPS bfs_ms sssp_ms bc_ms_total
    rm -f $O/b_dzig.csv
    step $D dzig $O/dzig.log timeout $T $B/dzig/run_dzig.sh $P/$D $D $O/b_dzig.csv $MAXT 10000 10
    for t in ${SWEEP//,/ }; do
        [ $t -eq $MAXT ] && continue
        rm -f $O/b_sweep_terrace_$t.csv $O/b_sweep_aspen_$t.csv $O/b_sweep_spruce_$t.csv
        step $D sweep_terrace_$t $O/sweep_terrace_$t.log env CILK_NWORKERS=$t QUERIES_ONLY=1 NO_BC=1 timeout $T $B/terrace/terrace_driver $P/$D $D $O/b_sweep_terrace_$t.csv 10000
        step $D sweep_aspen_$t $O/sweep_aspen_$t.log env NUM_THREADS=$t QUERIES_ONLY=1 NO_BC=1 timeout $T $B/aspen/aspen_driver $P/${D}_sym ${D}_sym $O/b_sweep_aspen_$t.csv 10000
        step $D sweep_spruce_$t $O/sweep_spruce_$t.log env OMP_NUM_THREADS=$t QUERIES_ONLY=1 NO_BC=1 NO_SSSP=1 timeout $T $B/spruce/spruce_driver $P/$D $D $O/b_sweep_spruce_$t.csv 10000
    done
    cat $O/b_*.csv > $O/baselines.csv 2>/dev/null
    log "== $D done"
}

log "start; waiting for a quiet machine"
for D in twitter_combined web-Google wiki-topcats cit-Patents livejournal; do run_dataset $D 1,2,4,8,16,32,64 0; done
run_dataset twitter 1,2,8,32,64 1
run_dataset friendster 1,2,8,32,64 1
cd $B && PAPER_DIR=$OUT PAPER_OUT=$G/experiments/results/paper_comparison_clean.md python3 make_paper_comparison.py
log "ALL DONE"
