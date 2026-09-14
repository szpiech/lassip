# Branch `perf-optimizations`

Performance, reproducibility and CLI work on top of v1.2.2 (`ba8ca63`), one
reviewable commit per change, `make check` green at every commit.

All measurements are on `example/YRI.chr22.vcf.gz` — 278,607 SNPs, 108
individuals / 216 phased haplotypes, `--winsize 117 --winstep 12` → 23,208
windows, K = 10 — on an 8-core arm64 machine, built with `make` (`-O3`, no
architecture flags).

## Benchmark

| | v1.2.2 | this branch |
|---|---|---|
| stage 1 (`--calc-spec --hapstats`), 1 thread | 29.9 s | **4.4 s** |
| stage 1, 8 threads | 12.4 s wall / 41.8 s CPU | **2.2 s wall / 5.5 s CPU** |
| stage 2 `--lassi`, 8 threads | 0.34 s | 0.34 s |
| stage 2 `--salti`, 8 threads | ≈48 min wall / ≈6.4 h CPU | **15.5 s wall / 118 s CPU** |
| stage 2 `--salti`, peak memory | 1,009 MB | **48 MB** |

Stage-1 spectra are byte-identical to v1.2.2 in every configuration tested.

The saltiLASSI figures for v1.2.2 are extrapolated rather than run to
completion: measured CPU per window is 0.78 s at 1,000 windows, 0.89 s at
2,000 and 0.99 s at 4,000, rising only because the fraction of windows with a
truncated flanking region shrinks (the median window has 143 neighbours within
the default `--max-extend-bp 100000`, and 140 by 4,000 windows). Taking the
saturated 1.0 s/window gives ≈6.4 h of CPU for the full contig, i.e. about
200× the 118 s measured here.

## Commits

| commit | change |
|---|---|
| `94158f0` | regression suite: 17 checks over both stages, `tests/run_tests.sh` |
| `8d8eaaf` | Makefile: header dependencies, `ARCHFLAGS`, `make check` |
| `6539db8` | `--seed`: make runs on data with missing genotypes reproducible |
| `2aa0f19` | remove 1,357 lines of unreachable code |
| `fecd9d9` | saltiLASSI: hoist `log`/`exp` out of the (A, m, eps) grid search |
| `fc93ed5` | saltiLASSI: one sweep-spectrum table per contig, not per window |
| `380d7f8` | `hfs_window`: skip the clustering pass when it cannot merge anything |
| `03cda7c` | VCF reader: no per-genotype map lookups |
| `dd18692` | filter loci in one pass instead of two matrix copies |
| `15bc273` | drop the leftover per-haplotype debug write |
| `faeb0c6` | CLI: help on stdout and exit 0, `--version`, usage, exit codes |
| `18da792` | CLI: group the help by topic |
| `f05ea0b` | CLI: `isFlagSet` instead of sentinel default strings |
| `50ecaa2` | CLI: validate `--dist-type` at both stages, fix its label, use the file's K |
| `1954260` | docs: refresh README, fix `example/do_lassip_YRI.bash` |
| `46a1e42` | `--match-tol`: group at `<= MATCH_TOL` differences, as documented |
| `edeae4f` | example: regenerate the committed outputs, fix the second command |

## Behavioural differences

Everything below is a deliberate change in output or behaviour; nothing else
changed.

1. **saltiLASSI `A`.** `m` and the likelihood ratio are unchanged in every
   window tested. `A` moves to an *adjacent* point on the 101-point log-A grid
   in about 1% of windows (1/94, 0/94 and 1/94 in the three `--salti` test
   cases). The likelihood is flat in `A` there, so which grid point attains
   the maximum is decided by rounding; v1.2.2 has the same instability and
   moves a different subset of windows under any change to summation order or
   compiler flags. This is noted in `calcMTA`, and the fix belongs at the model
   level — report the interval of `A` within some ΔL of the optimum, or adopt
   an explicit tie-breaking rule.
2. **Runs on data with missing genotypes now reproduce.** Previously four
   identical single-threaded runs gave four different spectra, at the default
   `--match-tol 0`. Output for such data therefore differs from any particular
   v1.2.2 run — there was no stable value to preserve. `--seed 0` restores the
   old clock-seeded behaviour.
3. **`--calc-spec --dist-type cm`** used to write a header whose 5th column had
   no name; that column is now labelled `ppos`, which is what it has always
   contained.
4. **The null-spectrum check** at the `--spectra` stage now uses the K recorded
   in the spectra file rather than the `--k` on the command line, so its
   threshold is no longer wrong when they differ.
5. **`--filter-level 2` stderr**: one `Filtering N loci in POP` line per
   population instead of two (missingness and monomorphism are now one pass).
6. **Flags belonging to the other stage now print a warning** instead of being
   silently ignored.
7. **Exit codes**: 64 for a command line or validation error, 65 for bad input
   data, 70 for an unexpected exception; `--help` and `--version` exit 0. An
   error raised inside a reader previously escaped `main` and aborted the
   process.
8. **`--match-tol` groups at `<=` the given number of differences**, as its
   help text says, instead of `<`. The new `--match-tol t` reproduces the old
   `t+1`; `--match-tol 0` on data without missing genotypes is unchanged.

## Left for you to decide

- **Order-dependence of haplotype clustering.** `--seed` makes the order
  reproducible; it does not remove the dependence. Merging into the *most
  frequent* compatible haplotype, or union-find over the ≤tol graph, would
  remove the arbitrariness, but both change the method.

Resolved since the first version of this file:

- **`--match-tol` semantics** (commit `46a1e42`). The merge test is now
  `d <= MATCH_TOL`, matching the flag's help text. Every setting previously
  behaved as the one below it; the new `--match-tol t` reproduces the old
  `t+1` exactly. Unchanged at `--match-tol 0` without missing genotypes.
  Note the README's 06FEB2025 entry states the old rule and so disagreed with
  the flag help; the help was taken as the intent.
- **Example outputs** (commit `edeae4f`). Regenerated, and the script's second
  command corrected to read the per-population filename the first command
  writes at the default `--filter-level 2`.

## Not done

- The VCF is still decompressed and scanned twice (once to count loci, once to
  parse); removing the counting pass needs a growable genotype matrix.
- Genotypes are still one `char` per haplotype per locus for a 2-bit alphabet;
  packing would cut memory 4× and make comparisons word-parallel.
- `hfs_window` still rebuilds each haplotype string from scratch, so with
  `--winsize 117 --winstep 12` every site is re-read into a new `std::string`
  about ten times; an incremental or rolling-hash HFS would remove that.
- Threading still uses one `pthread_t` per thread per stage with a static
  stride and a cast-to-`void*(*)(void*)` worker, and work orders are leaked.
- `main` is still one long function; population data is still threaded through
  half a dozen parallel `map<string, T*>`.
