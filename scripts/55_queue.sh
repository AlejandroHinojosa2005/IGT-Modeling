#!/bin/bash
# Runs the hierarchical fits (scripts/55_stan_fit.R) on this computer, chained so they keep going unattended (nights,
# weekends) until every fit is done, and can be paused whenever the computer is needed for something else.
# Usage (from the project folder):
#   bash scripts/55_queue.sh start [JOBS] [DATA ...]  start in the background. JOBS: fits at once, each on 4 cores
#                                                     (default: performance cores / 4). DATA: default replication, then
#                                                     horstmann; e.g. "start 2 horstmann" for Horstmann only.
#   bash scripts/55_queue.sh status                   what is done, what is running, the end of the log
#   bash scripts/55_queue.sh pause                    freeze the running fits at once (they keep their place) and start
#                                                     no new ones; the computer is free for other work
#   bash scripts/55_queue.sh resume                   carry on from where they were
#   bash scripts/55_queue.sh stop                     start no new fits; the running ones finish and are saved
#   bash scripts/55_queue.sh kill                     stop now; the running fits are lost and redone next time
#   bash scripts/55_queue.sh test                     one quick fit to simulated data, to check the setup
# Order: for each data set, the cross-validation fits of NOVA and ORL + alpha, then the other models' folds, then the
# fits to all trials; after each round, the reruns of fits that did not pass (scripts/55e_stan_jobs.R lists what is
# left). Saved fits are never redone, so start can be run again at any time, and Kodiak can take over whatever is
# left (kodiak/stan_cv.sh skips saved fits too). The log is logs/stan_queue.log; each fit's own output is in
# logs/stan_<data>_<model>_<fold>.log. On a Mac the queue keeps the computer awake (caffeinate) while it runs; keep it
# plugged in, with the lid open.
# The registration says the code is committed before any fit is run, so start refuses to run with uncommitted changes
# to the scripts, the Stan files or predictions/hierarchical.md.

cd "$(dirname "$0")/.." || exit 1
LOG=logs/stan_queue.log
PAUSE=logs/stan_queue.PAUSE
STOP=logs/stan_queue.STOP
mkdir -p logs
FIT_PROCS="file=scripts/55_stan_fit.R"                 # the R process of each fit
STAN_PROCS="stan/hier/igt_.*method=sample"             # its four Stan chains
CODE_FILES="scripts/55_stan_common.R scripts/55_stan_fit.R scripts/55c_stan_cv_score.R scripts/55d_stan_predictions.R
            scripts/55e_stan_jobs.R scripts/55_queue.sh stan/hier predictions/hierarchical.md"

log() { echo "$(date '+%a %b %d %H:%M')  $*" >> "$LOG"; }
queue_pids() { pgrep -f "^bash scripts/55_queue.sh _run" ; }
default_jobs() {
  local cores
  cores=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)
  [ "$cores" -ge 8 ] && echo $((cores / 4)) || echo 1
}
check_committed() {
  git rev-parse --git-dir > /dev/null 2>&1 || return 0          # not a git checkout: nothing to check
  local f untracked=""
  for f in $CODE_FILES; do git ls-files --error-unmatch "$f" > /dev/null 2>&1 || untracked="$untracked $f"; done
  if [ -n "$untracked" ] || ! git diff --quiet HEAD -- $CODE_FILES; then
    echo "Commit the scripts, the Stan files and predictions/hierarchical.md first: the registration says the code is"
    echo "committed before any fit is run. Uncommitted:"
    git status --short --untracked-files=all -- $CODE_FILES
    exit 1
  fi
}

case "$1" in
start)
  if [ -n "$(queue_pids)" ]; then echo "The queue is already running (bash scripts/55_queue.sh status)."; exit 1; fi
  check_committed
  JOBS=${2:-$(default_jobs)}
  if [ $# -ge 2 ]; then shift 2; else shift 1; fi
  DATA=${*:-replication horstmann}
  rm -f "$STOP" "$PAUSE"
  log "queue started: $DATA, $JOBS fits at a time, code at commit $(git rev-parse --short HEAD 2> /dev/null)"
  if command -v caffeinate > /dev/null; then
    nohup caffeinate -is bash scripts/55_queue.sh _run "$JOBS" $DATA >> "$LOG" 2>&1 &
  else
    nohup bash scripts/55_queue.sh _run "$JOBS" $DATA >> "$LOG" 2>&1 &
  fi
  sleep 1
  echo "Started: $DATA, $JOBS fits at a time. Check on it with: bash scripts/55_queue.sh status"
  ;;

_run)                                    # the queue itself (started by start, or by scripts/55_run_all.sh)
  JOBS=$2
  shift 2
  Rscript -e 'library(cmdstanr); for (f in Sys.glob("stan/hier/igt_*.stan")) invisible(cmdstan_model(f))' > /dev/null 2>&1 ||
    { log "the Stan models did not compile; run: Rscript -e 'cmdstanr::cmdstan_model(\"stan/hier/igt_nova.stan\")'"; exit 1; }
  for DATA in "$@"; do
    PREV=""
    while :; do
      [ -f "$STOP" ] && { log "stopped"; exit 0; }
      LIST=$(Rscript scripts/55e_stan_jobs.R "$DATA") || { log "could not list the fits left for $DATA"; exit 1; }
      if [ -z "$LIST" ]; then log "$DATA: every fit is done"; break; fi
      if [ "$LIST" = "$PREV" ]; then
        log "$DATA: these fits failed twice in a row, so the queue moves on (see their logs/stan_${DATA}_... files): $(echo $LIST)"
        break
      fi
      PREV=$LIST
      log "$DATA: $(printf '%s\n' "$LIST" | wc -l | tr -d ' ') fits in this round, $JOBS at a time"
      printf '%s\n' "$LIST" | xargs -L 1 -P "$JOBS" bash scripts/55_queue.sh _job "$DATA"
    done
  done
  log "queue finished"
  ;;

_job)                                    # one fit (called by _run through xargs)
  DATA=$2; MODEL=$3; FOLD=$4; EXTRA=$5
  while [ -f "$PAUSE" ]; do sleep 30; done
  [ -f "$STOP" ] && exit 0
  NAME="$DATA $MODEL $FOLD${EXTRA:+ $EXTRA}"
  FITLOG="logs/stan_${DATA}_${MODEL}_${FOLD}${EXTRA:+_$EXTRA}.log"
  log "start  $NAME"
  T0=$(date +%s)
  if Rscript scripts/55_stan_fit.R "$DATA" "$MODEL" "$FOLD" $EXTRA > "$FITLOG" 2>&1; then
    log "done   $NAME: $(( ($(date +%s) - T0) / 60 )) min, $(grep -h 'the convergence rule' "$FITLOG" | head -1)"
  else
    log "FAILED $NAME: see $FITLOG"
  fi
  exit 0
  ;;

status)
  if [ -n "$(queue_pids)" ]; then
    if [ -f "$PAUSE" ]; then echo "The queue is PAUSED (bash scripts/55_queue.sh resume)."
    elif [ -f "$STOP" ]; then echo "The queue is stopping: the running fits finish, no new ones start."
    else echo "The queue is running."; fi
  else
    echo "The queue is not running."
  fi
  Rscript scripts/55e_stan_jobs.R summary horstmann replication
  RUNNING=$(ps -axo etime=,command= | grep "file=scripts/55_stan_fit.R --args" | grep -v grep | sed 's/ *\([^ ]*\) .*--args \(.*\)/  \2  (running for \1)/')
  [ -n "$RUNNING" ] && { echo "Running now:"; echo "$RUNNING"; }
  echo "End of the log ($LOG):"
  tail -n 12 "$LOG" 2> /dev/null
  ;;

pause)
  touch "$PAUSE"
  pkill -STOP -f "$STAN_PROCS"
  pkill -STOP -f "$FIT_PROCS"
  log "paused"
  echo "Paused. The running fits are frozen in place and no new ones will start. To carry on: bash scripts/55_queue.sh resume"
  ;;

resume)
  rm -f "$PAUSE"
  pkill -CONT -f "$FIT_PROCS"
  pkill -CONT -f "$STAN_PROCS"
  log "resumed"
  echo "Resumed."
  ;;

stop)
  touch "$STOP"
  log "stop requested: the running fits will finish"
  echo "No new fits will start; the running ones finish and are saved. Start again later with: bash scripts/55_queue.sh start"
  ;;

kill)
  touch "$STOP"
  rm -f "$PAUSE"
  pkill -f "^caffeinate -is bash scripts/55_queue.sh _run"
  pkill -f "^bash scripts/55_queue.sh _run"
  pkill -f "^xargs .*scripts/55_queue.sh _job"
  pkill -f "^bash scripts/55_queue.sh _job"
  pkill -CONT -f "$FIT_PROCS"; pkill -CONT -f "$STAN_PROCS"
  pkill -f "$FIT_PROCS"
  pkill -f "$STAN_PROCS"
  log "killed: the fits that were running will be redone next time"
  echo "Stopped. The fits that were running are lost and will be redone next time."
  ;;

check)                                   # used by scripts/55_run_all.sh
  check_committed
  ;;

test)
  Rscript scripts/55_stan_fit.R simulated orl_alpha 1 quick force &&
    Rscript scripts/55_stan_fit.R simulated orl_alpha 1 quick long force
  ;;

*)
  sed -n 2,22p "$0"
  ;;
esac
