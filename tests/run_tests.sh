#!/bin/bash
# lassip regression suite.
#
#   ./tests/run_tests.sh                 # check the build in src/ against tests/expected/
#   ./tests/run_tests.sh --regen         # re-record tests/expected/ from the current build
#   ./tests/run_tests.sh --bin /path/to/lassip
#   ./tests/run_tests.sh --keep          # keep the scratch directory for inspection
#   ./tests/run_tests.sh --list          # list case names and exit
#   ./tests/run_tests.sh <case> [...]    # run only the named cases
#
# Fixtures are derived at run time from data already in the repository (example/
# and testing/), so no new test data is committed. Stage-1 outputs are compared
# by hash of the decompressed text; stage-2 outputs are compared column by
# column, because a hash would hide which statistic moved.
#
# Regenerating: goldens must be recorded from a build you trust. The intended
# workflow is to record them once from the reference build, then require every
# later commit to reproduce them.

set -u

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$TESTS_DIR/.." && pwd)"
BIN="$ROOT_DIR/src/lassip"
EXPECTED="$TESTS_DIR/expected"
REGEN=0
KEEP=0
LIST=0
SELECTED=()
THREADS=${LASSIP_TEST_THREADS:-2}

while [ $# -gt 0 ]; do
    case "$1" in
        --regen) REGEN=1; shift ;;
        --keep)  KEEP=1; shift ;;
        --list)  LIST=1; shift ;;
        --bin)   BIN="$2"; shift 2 ;;
        --help|-h)
            sed -n '2,20p' "$BASH_SOURCE"; exit 0 ;;
        -*) echo "unknown option: $1" >&2; exit 64 ;;
        *)  SELECTED+=("$1"); shift ;;
    esac
done

# ---------------------------------------------------------------- utilities

if command -v md5sum >/dev/null 2>&1; then
    hash_cmd() { md5sum | awk '{print $1}'; }
elif command -v md5 >/dev/null 2>&1; then
    hash_cmd() { md5 -q; }
else
    echo "ERROR: need md5sum or md5 on PATH." >&2; exit 70
fi

hash_gz() { gzip -dc "$1" | hash_cmd; }

NPASS=0; NFAIL=0; NSKIP=0
FAILED_CASES=()

pass() { printf '  \033[32mok\033[0m    %s\n' "$1"; NPASS=$((NPASS+1)); }
fail() { printf '  \033[31mFAIL\033[0m  %s: %s\n' "$1" "$2"; NFAIL=$((NFAIL+1)); FAILED_CASES+=("$1"); }
skip() { printf '  \033[33mskip\033[0m  %s: %s\n' "$1" "$2"; NSKIP=$((NSKIP+1)); }

# compare_hash <case> <actual.gz> <golden-name>
compare_hash() {
    local name="$1" actual="$2" golden="$EXPECTED/$3.md5"
    if [ ! -f "$actual" ]; then fail "$name" "no output at $actual"; return; fi
    local h; h="$(hash_gz "$actual")"
    if [ "$REGEN" = 1 ]; then
        echo "$h" > "$golden"; pass "$name (recorded $3.md5)"; return
    fi
    if [ ! -f "$golden" ]; then fail "$name" "no golden $3.md5 (run --regen)"; return; fi
    if [ "$h" = "$(cat "$golden")" ]; then pass "$name"
    else fail "$name" "hash mismatch vs $3.md5"; fi
}

# compare_table <case> <actual.gz> <golden-name> <rtol> [tolerant-col-regex] [max-frac-differing]
#
# Columns whose name matches the tolerant regex may differ in up to the given
# fraction of rows without failing; every other column must agree within rtol
# (integer-looking columns must agree exactly). Header must match exactly.
compare_table() {
    local name="$1" actual="$2" gname="$3" rtol="$4"
    local tolcols="${5:-}" maxfrac="${6:-0}"
    local golden="$EXPECTED/$gname.tsv.gz"
    if [ ! -f "$actual" ]; then fail "$name" "no output at $actual"; return; fi
    if [ "$REGEN" = 1 ]; then
        cp "$actual" "$golden"; pass "$name (recorded $gname.tsv.gz)"; return
    fi
    if [ ! -f "$golden" ]; then fail "$name" "no golden $gname.tsv.gz (run --regen)"; return; fi

    local report
    report=$(awk -v rtol="$rtol" -v tolcols="$tolcols" -v maxfrac="$maxfrac" '
        function isnum(x) { return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) }
        function isint(x) { return (x ~ /^[+-]?[0-9]+$/) }
        FNR==NR {
            if (FNR==1) { ehdr=$0; en=NF; for(i=1;i<=NF;i++) ename[i]=$i }
            else { for(i=1;i<=NF;i++) E[FNR,i]=$i; erows=FNR }
            next
        }
        {
            if (FNR==1) {
                if ($0 != ehdr) { print "header differs"; bad=1; exit }
                next
            }
            if (FNR > erows) { print "actual has more rows than golden"; bad=1; exit }
            for (i=1;i<=en;i++) {
                e=E[FNR,i]; a=$i
                if (e==a) continue
                if (!isnum(e) || !isnum(a)) { diffcol[i]++; hard[i]++; continue }
                if (isint(e) && isint(a)) { diffcol[i]++; hard[i]++; continue }
                d = (e==0) ? (a==0?0:1) : (e-a)/e; if (d<0) d=-d
                if (d > maxd[i]) maxd[i]=d
                if (d > rtol) { diffcol[i]++; if (tolcols=="" || ename[i] !~ tolcols) hard[i]++ }
            }
            arows=FNR
        }
        END {
            if (bad) exit
            if (arows != erows) { print "row count differs: golden " erows-1 " actual " arows-1; exit }
            nrow = erows - 1
            msg=""
            for (i=1;i<=en;i++) {
                if (!diffcol[i]) continue
                frac = diffcol[i]/nrow
                if (hard[i] > 0) {
                    msg = msg sprintf("%s: %d/%d rows differ (max rel %.3g); ", ename[i], diffcol[i], nrow, maxd[i])
                } else if (frac > maxfrac+1e-12) {
                    msg = msg sprintf("%s: %d/%d rows differ (%.1f%% > %.1f%% allowed); ", ename[i], diffcol[i], nrow, 100*frac, 100*maxfrac)
                }
            }
            if (msg != "") print msg
        }
    ' <(gzip -dc "$golden") <(gzip -dc "$actual"))

    if [ -z "$report" ]; then pass "$name"; else fail "$name" "$report"; fi
}

run_lassip() {
    local logfile="$1"; shift
    if ! "$BIN" "$@" > "$logfile" 2>&1; then
        echo "--- lassip failed: $BIN $* ---" >&2
        tail -20 "$logfile" >&2
        return 1
    fi
    return 0
}

selected() {
    [ ${#SELECTED[@]} -eq 0 ] && return 0
    local c; for c in "${SELECTED[@]}"; do [ "$c" = "$1" ] && return 0; done
    return 1
}

# ---------------------------------------------------------------- fixtures

CASES="spec_phased stats_only spec_unphased spec_twopop spec_filter1 spec_filter0
       spec_missing spec_missing_tol missing_determinism nullwin_threads
       cluster_bestcomp cluster_softem cluster_conserved
       spec_unphased_missing cluster_unphased salti_unphased
       spec_medium avg_spec lassi lassi_nullspec salti_bp salti_nw salti_cm
       cm_map_mismatch cm_max_gap"

if [ "$LIST" = 1 ]; then for c in $CASES; do echo "$c"; done; exit 0; fi

if [ ! -x "$BIN" ]; then echo "ERROR: no lassip binary at $BIN (run make first)." >&2; exit 70; fi

SMALL="$ROOT_DIR/testing/small.vcf.gz"
YRI="$ROOT_DIR/example/YRI.chr22.vcf.gz"
for f in "$SMALL" "$YRI"; do
    [ -f "$f" ] || { echo "ERROR: missing repository test input $f" >&2; exit 70; }
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/lassip-tests.XXXXXX")"
cleanup() { if [ "$KEEP" = 1 ]; then echo "scratch kept at $WORK"; else rm -rf "$WORK"; fi; }
trap cleanup EXIT

mkdir -p "$EXPECTED"

# one population containing every sample in small.vcf.gz
gzip -dc "$SMALL" | awk '/^#CHROM/ {for(i=10;i<=NF;i++) print $i"\tPOP1"; exit}' > "$WORK/small.pop1.txt"
# the same samples split into two populations, in header order
gzip -dc "$SMALL" | awk '/^#CHROM/ {n=NF-9; for(i=10;i<=NF;i++) print $i"\t"((i-9)<=int(n/2)?"POPA":"POPB"); exit}' > "$WORK/small.pop2.txt"
# every sample of the YRI example
gzip -dc "$YRI" | awk '/^#CHROM/ {for(i=10;i<=NF;i++) print $i"\tYRI"; exit}' > "$WORK/yri.pop.txt"

# deterministic missing-data fixture: on every 3rd record, blank every 7th sample
gzip -dc "$SMALL" | awk 'BEGIN{OFS="\t"}
    /^#/ {print; next}
    { r++; if (r%3==0) for(i=10;i<=NF;i++) if ((i-9)%7==0) $i="./."; print }' | gzip > "$WORK/small.missing.vcf.gz"

# genetic map for the --dist-type cm path: 1 cM per Mb over the small fixture
gzip -dc "$SMALL" | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, ($3=="."?"locus"NR:$3), $2/1000000.0, $2}' > "$WORK/small.map"

# the same map under a contig name the spectra do not use: getMapInfo used to
# index a map<> entry it had just default-created and dereference a null pointer
awk 'BEGIN{OFS="\t"} {$1="nosuchchr"; print}' "$WORK/small.map" > "$WORK/small.badchr.map"

# a 20k-SNP slice of the YRI example for a realistic multi-window stage-1 case
gzip -dc "$YRI" | awk '/^#/{print;next} {n++; if(n<=20000) print; else exit}' | gzip > "$WORK/yri.slice.vcf.gz"

# every fifth record fully missing: with --max-lmiss 1 --keep-monomorphic those
# loci survive filtering, and at --max-hmiss 0 every haplotype of every window
# covering one is dropped, so every window is null. Exercises the null-window
# counter, which used to be incremented from every thread unsynchronised.
gzip -dc "$SMALL" | awk 'BEGIN{OFS="\t"} /^#/{print;next} \
    { r++; if (r%5==0) for(i=10;i<=NF;i++) $i="./."; print }' | gzip > "$WORK/small.allmissing.vcf.gz"

echo "lassip regression suite"
echo "  binary : $BIN"
echo "  data   : testing/small.vcf.gz, example/YRI.chr22.vcf.gz (+ derived fixtures)"
if [ "$REGEN" = 1 ]; then echo "  mode   : RECORDING goldens into $EXPECTED"; fi
echo

# ---------------------------------------------------------------- stage 1

if selected spec_phased; then
    if run_lassip "$WORK/spec_phased.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 --out "$WORK/sp"; then
        compare_hash spec_phased "$WORK/sp.POP1.lassip.hap.spectra.gz" spec_phased
    else fail spec_phased "run failed"; fi
fi

if selected stats_only; then
    if run_lassip "$WORK/stats_only.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" \
        --hapstats --winsize 50 --winstep 10 --out "$WORK/so"; then
        compare_hash stats_only "$WORK/so.POP1.lassip.hap.stats.gz" stats_only
    else fail stats_only "run failed"; fi
fi

if selected spec_unphased; then
    if run_lassip "$WORK/spec_unphased.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" \
        --unphased --calc-spec --hapstats --k 5 --winsize 50 --winstep 10 --out "$WORK/su"; then
        compare_hash spec_unphased "$WORK/su.POP1.lassip.mlg.spectra.gz" spec_unphased
    else fail spec_unphased "run failed"; fi
fi

if selected spec_twopop; then
    if run_lassip "$WORK/spec_twopop.log" --vcf "$SMALL" --pop "$WORK/small.pop2.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 --out "$WORK/tp"; then
        compare_hash spec_twopop_a "$WORK/tp.POPA.lassip.hap.spectra.gz" spec_twopop_a
        compare_hash spec_twopop_b "$WORK/tp.POPB.lassip.hap.spectra.gz" spec_twopop_b
    else fail spec_twopop "run failed"; fi
fi

if selected spec_filter1; then
    # --filter-level 1 pools populations and writes a single combined file
    if run_lassip "$WORK/spec_filter1.log" --vcf "$SMALL" --pop "$WORK/small.pop2.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 --filter-level 1 --out "$WORK/f1"; then
        compare_hash spec_filter1 "$WORK/f1.lassip.hap.spectra.gz" spec_filter1
    else fail spec_filter1 "run failed"; fi
fi

if selected spec_filter0; then
    if run_lassip "$WORK/spec_filter0.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --filter-level 0 --keep-monomorphic --out "$WORK/f0"; then
        compare_hash spec_filter0 "$WORK/f0.lassip.hap.spectra.gz" spec_filter0
    else fail spec_filter0 "run failed"; fi
fi

if selected spec_missing; then
    # Pinned to --hap-cluster garud-shuffle: this golden was recorded under the
    # pre-1.3 rule and keeping it proves that path is still bit-identical.
    # Deterministic since --seed exists: the per-window shuffle used when clustering
    # haplotypes that carry missing data is seeded from --seed mixed with the
    # window's SNP boundaries, so output no longer depends on the clock or on
    # --threads. Before that, four identical runs of this case gave four files.
    if run_lassip "$WORK/spec_missing.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 0 --hap-cluster garud-shuffle \
        --out "$WORK/ms"; then
        compare_hash spec_missing "$WORK/ms.POP1.lassip.hap.spectra.gz" spec_missing
    else fail spec_missing "run failed"; fi
fi

if selected spec_missing_tol; then
    # --match-tol > 0 actively merges haplotypes that differ at missing sites.
    # Pinned to the pre-1.3 rule, as above.
    if run_lassip "$WORK/spec_missing_tol.log" --vcf "$WORK/small.missing.vcf.gz" \
        --pop "$WORK/small.pop1.txt" --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --hap-cluster garud-shuffle \
        --out "$WORK/mt"; then
        compare_hash spec_missing_tol "$WORK/mt.POP1.lassip.hap.spectra.gz" spec_missing_tol
    else fail spec_missing_tol "run failed"; fi
fi

if selected missing_determinism; then
    # The property itself, independent of any golden: the same seed must give the
    # same spectrum whatever --threads is set to.
    ok=1
    for t in 1 3; do
        run_lassip "$WORK/det$t.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
            --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
            --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --threads "$t" --out "$WORK/det$t" || ok=0
    done
    if [ "$ok" = 1 ]; then
        h1=$(hash_gz "$WORK/det1.POP1.lassip.hap.spectra.gz")
        h3=$(hash_gz "$WORK/det3.POP1.lassip.hap.spectra.gz")
        if [ "$h1" = "$h3" ]; then pass "missing_determinism"
        else fail missing_determinism "--threads 1 and --threads 3 disagree on data with missing genotypes"; fi
    else fail missing_determinism "run failed"; fi
fi

if selected spec_medium; then
    if run_lassip "$WORK/spec_medium.log" --vcf "$WORK/yri.slice.vcf.gz" --pop "$WORK/yri.pop.txt" \
        --calc-spec --hapstats --k 10 --winsize 117 --winstep 60 --threads "$THREADS" --out "$WORK/med"; then
        compare_hash spec_medium "$WORK/med.YRI.lassip.hap.spectra.gz" spec_medium
    else fail spec_medium "run failed"; fi
fi

if selected cluster_bestcomp; then
    # The default rule. Also checks that omitting --hap-cluster selects it, and
    # that it ignores --seed: three runs, two of them with different seeds, must
    # all produce the same file.
    ok=1
    run_lassip "$WORK/cbc.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --hap-cluster best-comp \
        --out "$WORK/cbc" || ok=0
    run_lassip "$WORK/cbc_default.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --out "$WORK/cbc_default" || ok=0
    run_lassip "$WORK/cbc_seed.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 --seed 99 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --hap-cluster best-comp \
        --out "$WORK/cbc_seed" || ok=0
    if [ "$ok" = 0 ]; then fail cluster_bestcomp "run failed"
    elif [ "$(hash_gz "$WORK/cbc.POP1.lassip.hap.spectra.gz")" \
         != "$(hash_gz "$WORK/cbc_default.POP1.lassip.hap.spectra.gz")" ]; then
        fail cluster_bestcomp "best-comp is not the default"
    elif [ "$(hash_gz "$WORK/cbc.POP1.lassip.hap.spectra.gz")" \
         != "$(hash_gz "$WORK/cbc_seed.POP1.lassip.hap.spectra.gz")" ]; then
        fail cluster_bestcomp "--seed changed a best-comp run"
    else
        compare_hash cluster_bestcomp "$WORK/cbc.POP1.lassip.hap.spectra.gz" cluster_bestcomp
    fi
fi

if selected cluster_softem; then
    # Experimental: class sizes are fractional, so the spectrum columns are no
    # longer integers. Compared with a tolerance rather than by hash.
    if run_lassip "$WORK/cse.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
        --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --hap-cluster soft-em \
        --out "$WORK/cse"; then
        compare_table cluster_softem "$WORK/cse.POP1.lassip.hap.spectra.gz" cluster_softem 1e-9 "" 0
    else fail cluster_softem "run failed"; fi
fi

if selected cluster_conserved; then
    # Clustering decides how haplotypes are GROUPED, never how many there are:
    # nhaps is the sum of the class sizes, so it must agree across the three
    # rules window by window, and no window may report more classes than
    # haplotypes. This catches a rule that drops or double-counts a haplotype.
    ok=1
    for m in garud-shuffle best-comp soft-em; do
        run_lassip "$WORK/cc_$m.log" --vcf "$WORK/small.missing.vcf.gz" --pop "$WORK/small.pop1.txt" \
            --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 \
            --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 --hap-cluster "$m" \
            --out "$WORK/cc_$m" || ok=0
    done
    if [ "$ok" = 0 ]; then fail cluster_conserved "run failed"
    else
        report=$(for m in garud-shuffle best-comp soft-em; do
                     gzip -dc "$WORK/cc_$m.POP1.lassip.hap.spectra.gz" | awk -v m="$m" 'NR>2{print m, NR, $6, $7}'
                 done | awk '{ key=$2; n[key]=n[key]+1; if (n[key]==1) ref[key]=$3;
                               else if ($3 != ref[key]) mismatch++;
                               if ($4+0 > $3+0) toomany++ }
                             END{ printf "%d %d", mismatch+0, toomany+0 }')
        set -- $report
        if [ "$1" = 0 ] && [ "$2" = 0 ]; then pass "cluster_conserved (3 rules agree on nhaps)"
        else fail cluster_conserved "nhaps mismatches: $1, windows with uhaps > nhaps: $2"; fi
    fi
fi

if selected spec_unphased_missing; then
    # --unphased was only covered on data with no missing genotypes, so none of
    # the clustering code was exercised on multilocus genotype strings. Those
    # are drawn from {0,1,2,-} rather than {0,1,-}, which is the only path where
    # the packed representation uses all four symbols.
    if run_lassip "$WORK/spec_unphased_missing.log" --vcf "$WORK/small.missing.vcf.gz" \
        --pop "$WORK/small.pop1.txt" --unphased --calc-spec --hapstats --k 5 \
        --winsize 50 --winstep 10 --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 \
        --hap-cluster garud-shuffle --seed 3 --out "$WORK/umg"; then
        compare_hash spec_unphased_missing "$WORK/umg.POP1.lassip.mlg.spectra.gz" spec_unphased_missing
    else fail spec_unphased_missing "run failed"; fi
fi

if selected cluster_unphased; then
    # best-comp on unphased data with missing genotypes, at --match-tol 2 where
    # the three rules actually differ (at tol 0 they agree here: a 50-locus
    # three-state string is rarely ambiguous). Also asserts the result does not
    # move with --seed, and that all three rules conserve the genotype count.
    if run_lassip "$WORK/cluster_unphased.log" --vcf "$WORK/small.missing.vcf.gz" \
        --pop "$WORK/small.pop1.txt" --unphased --calc-spec --hapstats --k 5 \
        --winsize 50 --winstep 10 --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 \
        --out "$WORK/umb"; then
        compare_hash cluster_unphased "$WORK/umb.POP1.lassip.mlg.spectra.gz" cluster_unphased
        run_lassip "$WORK/cluster_unphased_seed.log" --vcf "$WORK/small.missing.vcf.gz" \
            --pop "$WORK/small.pop1.txt" --unphased --calc-spec --hapstats --k 5 \
            --winsize 50 --winstep 10 --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 \
            --seed 99 --out "$WORK/umb99"
        if [ "$(hash_gz "$WORK/umb.POP1.lassip.mlg.spectra.gz")" = \
             "$(hash_gz "$WORK/umb99.POP1.lassip.mlg.spectra.gz")" ]; then
            pass "cluster_unphased --seed has no effect"
        else fail cluster_unphased "--seed changed a best-comp unphased run"; fi
        ok=1
        for m in garud-shuffle best-comp soft-em; do
            run_lassip "$WORK/ucc_$m.log" --vcf "$WORK/small.missing.vcf.gz" \
                --pop "$WORK/small.pop1.txt" --unphased --calc-spec --hapstats --k 5 \
                --winsize 50 --winstep 10 --max-lmiss 0.5 --max-hmiss 0.5 --match-tol 2 \
                --hap-cluster "$m" --out "$WORK/ucc_$m" || ok=0
        done
        if [ "$ok" = 0 ]; then fail cluster_unphased "conservation run failed"
        else
            report=$(for m in garud-shuffle best-comp soft-em; do
                         gzip -dc "$WORK/ucc_$m.POP1.lassip.mlg.spectra.gz" | awk -v m="$m" 'NR>2{print m, NR, $6, $7}'
                     done | awk '{ key=$2; n[key]=n[key]+1; if (n[key]==1) ref[key]=$3;
                                   else if ($3 != ref[key]) mismatch++;
                                   if ($4+0 > $3+0) toomany++ }
                                 END{ printf "%d %d", mismatch+0, toomany+0 }')
            set -- $report
            if [ "$1" = 0 ] && [ "$2" = 0 ]; then pass "cluster_unphased (3 rules agree on nhaps)"
            else fail cluster_unphased "nhaps mismatches: $1, windows with uhaps > nhaps: $2"; fi
        fi
    else fail cluster_unphased "run failed"; fi
fi

# ---------------------------------------------------------------- stage 2

# all stage-2 cases read this spectrum
SPEC="$WORK/sp.POP1.lassip.hap.spectra.gz"

# an unphased spectra file, so stage 2 is covered on the .mlg path too
run_lassip "$WORK/spu.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" --unphased \
    --calc-spec --k 5 --winsize 50 --winstep 10 --out "$WORK/spu" >/dev/null 2>&1
SPECU="$WORK/spu.POP1.lassip.mlg.spectra.gz"
if [ ! -f "$SPEC" ]; then
    run_lassip "$WORK/spec_for_stage2.log" --vcf "$SMALL" --pop "$WORK/small.pop1.txt" \
        --calc-spec --hapstats --k 10 --winsize 50 --winstep 10 --out "$WORK/sp" || true
fi

if selected avg_spec && [ -f "$SPEC" ]; then
    if run_lassip "$WORK/avg_spec.log" --spectra "$SPEC" --avg-spec --out "$WORK/av"; then
        compare_hash avg_spec "$WORK/av.lassip.null.spectra.gz" avg_spec
    else fail avg_spec "run failed"; fi
fi

if selected lassi && [ -f "$SPEC" ]; then
    if run_lassip "$WORK/lassi.log" --spectra "$SPEC" --lassi --threads "$THREADS" --out "$WORK/la"; then
        compare_table lassi "$WORK/la.lassip.hap.out.gz" lassi 1e-6
    else fail lassi "run failed"; fi
fi

if selected lassi_nullspec && [ -f "$SPEC" ]; then
    run_lassip "$WORK/nullspec_make.log" --spectra "$SPEC" --avg-spec --out "$WORK/nsrc" || true
    if run_lassip "$WORK/lassi_nullspec.log" --spectra "$SPEC" --lassi \
        --null-spec "$WORK/nsrc.lassip.null.spectra.gz" --threads "$THREADS" --out "$WORK/ln"; then
        compare_table lassi_nullspec "$WORK/ln.lassip.hap.out.gz" lassi_nullspec 1e-6
    else fail lassi_nullspec "run failed"; fi
fi

# The saltiLASSI A column is the argmax of a likelihood that is flat in A over
# part of its range, so which of two adjacent grid points wins is decided by
# rounding and moves under any change to summation order or compiler flags.
# m and L must match to 1e-6; A may move in up to 5% of windows.
SALTI_TOL_COLS='_A$'
SALTI_TOL_FRAC=0.05

if selected salti_bp && [ -f "$SPEC" ]; then
    if run_lassip "$WORK/salti_bp.log" --spectra "$SPEC" --salti --dist-type bp \
        --max-extend-bp 200000 --threads "$THREADS" --out "$WORK/sb"; then
        compare_table salti_bp "$WORK/sb.lassip.hap.out.gz" salti_bp 1e-6 "$SALTI_TOL_COLS" "$SALTI_TOL_FRAC"
    else fail salti_bp "run failed"; fi
fi

if selected salti_nw && [ -f "$SPEC" ]; then
    if run_lassip "$WORK/salti_nw.log" --spectra "$SPEC" --salti --dist-type nw \
        --max-extend-nw 5 --threads "$THREADS" --out "$WORK/sn"; then
        compare_table salti_nw "$WORK/sn.lassip.hap.out.gz" salti_nw 1e-6 "$SALTI_TOL_COLS" "$SALTI_TOL_FRAC"
    else fail salti_nw "run failed"; fi
fi

if selected salti_cm && [ -f "$SPEC" ]; then
    if run_lassip "$WORK/salti_cm.log" --spectra "$SPEC" --salti --dist-type cm \
        --map "$WORK/small.map" --max-extend-cm 0.2 --threads "$THREADS" --out "$WORK/sc"; then
        compare_table salti_cm "$WORK/sc.lassip.hap.out.gz" salti_cm 1e-6 "$SALTI_TOL_COLS" "$SALTI_TOL_FRAC"
    else fail salti_cm "run failed"; fi
fi

if selected nullwin_threads; then
    # The window count in the spectra header is windows minus null windows, and
    # the finalize stage reads it to size its arrays. Incrementing nullWins from
    # every thread without synchronisation made it wrong and different on every
    # run: four 8-thread runs on this fixture reported 88, 56, 64 and 95 windows
    # where the answer is 0.
    ok=1
    for t in 1 8 8; do
        run_lassip "$WORK/nullwin_$t.log" --vcf "$WORK/small.allmissing.vcf.gz" --pop "$WORK/small.pop1.txt" \
            --calc-spec --k 10 --winsize 8 --winstep 1 --max-lmiss 1 --max-hmiss 0 \
            --keep-monomorphic --threads "$t" --out "$WORK/nw$t" || ok=0
    done
    if [ "$ok" = 1 ]; then
        h1=$(gzip -dc "$WORK/nw1.POP1.lassip.hap.spectra.gz" | head -1)
        h8=$(gzip -dc "$WORK/nw8.POP1.lassip.hap.spectra.gz" | head -1)
        if [ "$h1" = "$h8" ]; then pass "nullwin_threads ($(echo "$h1" | grep -o 'wins [0-9]*'))"
        else fail nullwin_threads "window count differs between 1 and 8 threads: [$h1] vs [$h8]"; fi
    else fail nullwin_threads "run failed"; fi
fi

if selected cm_map_mismatch && [ -f "$SPEC" ]; then
    # A genetic map that does not cover the spectra must be reported, not
    # crashed on: before the fix this segfaulted (exit 139) because the contig
    # lookup default-created a null entry and then dereferenced it. Windows in
    # a gap wider than the map's MAXGAP take the same path; they used to be
    # silently assigned the last successfully placed window's genetic position.
    # run the binary directly: run_lassip flattens every failure to 1, and the
    # exit code is the thing under test
    "$BIN" --spectra "$SPEC" --salti --dist-type cm --map "$WORK/small.badchr.map" \
        --max-extend-cm 0.2 --threads "$THREADS" --out "$WORK/cmbad" > "$WORK/cm_mismatch.log" 2>&1
    rc=$?
    if [ "$rc" = 65 ] && grep -q "does not place" "$WORK/cm_mismatch.log"; then
        pass "cm_map_mismatch (exit 65, reported)"
    else
        fail cm_map_mismatch "expected exit 65 with a 'does not place' error, got exit $rc"
    fi
fi

if selected cm_max_gap && [ -f "$SPEC" ]; then
    # --max-gap bounds how far --dist-type cm will interpolate. A tiny value
    # puts every window in an over-wide gap, which must be reported rather than
    # guessed; 0 means no limit and must reproduce the default run, since this
    # map has no gap anywhere near the 3 Mb default.
    "$BIN" --spectra "$SPEC" --salti --dist-type cm --map "$WORK/small.map" \
        --max-extend-cm 0.2 --max-gap 100 --threads "$THREADS" --out "$WORK/gapsmall" \
        > "$WORK/cm_max_gap.log" 2>&1
    rc=$?
    if [ "$rc" != 65 ] || ! grep -q "does not place" "$WORK/cm_max_gap.log"; then
        fail cm_max_gap "--max-gap 100 should exit 65 with a 'does not place' error, got exit $rc"
    elif run_lassip "$WORK/cm_max_gap0.log" --spectra "$SPEC" --salti --dist-type cm \
            --map "$WORK/small.map" --max-extend-cm 0.2 --max-gap 0 --threads "$THREADS" \
            --out "$WORK/gapnone"; then
        compare_table cm_max_gap "$WORK/gapnone.lassip.hap.out.gz" salti_cm 1e-6 "$SALTI_TOL_COLS" "$SALTI_TOL_FRAC"
    else fail cm_max_gap "--max-gap 0 run failed"; fi
fi

if selected salti_unphased && [ -f "$SPECU" ]; then
    # stage 2 reading a '#phased 0' spectra file: writes .mlg.out.gz rather than
    # .hap.out.gz, and nothing downstream of the header should care otherwise
    if run_lassip "$WORK/salti_unphased.log" --spectra "$SPECU" --salti --dist-type bp \
        --max-extend-bp 200000 --threads "$THREADS" --out "$WORK/su2"; then
        compare_table salti_unphased "$WORK/su2.lassip.mlg.out.gz" salti_unphased 1e-6 "$SALTI_TOL_COLS" "$SALTI_TOL_FRAC"
    else fail salti_unphased "run failed"; fi
fi

# ---------------------------------------------------------------- summary

echo
printf 'passed %d, failed %d, skipped %d\n' "$NPASS" "$NFAIL" "$NSKIP"
if [ "$NFAIL" -gt 0 ]; then
    printf 'failed: %s\n' "${FAILED_CASES[*]}"
    exit 1
fi
exit 0
