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
| stage 1 (`--calc-spec --hapstats`), 1 thread | 30.2 s | **1.8 s** |
| stage 1, 8 threads | 12.6 s | **1.3 s** |
| stage 1, peak memory | 229 MB | **87 MB** |
| stage 2 `--lassi`, 8 threads | 0.34 s | 0.34 s |
| stage 2 `--salti`, 8 threads | ≈48 min wall / ≈6.4 h CPU | **18 s / 118 s** |
| stage 2 `--salti`, peak memory | 1,009 MB | **53 MB** |

Stage-1 spectra are byte-identical to v1.2.2 in every configuration tested.

The saltiLASSI figures for v1.2.2 are extrapolated rather than run to
completion: measured CPU per window is 0.78 s at 1,000 windows, 0.89 s at
2,000 and 0.99 s at 4,000, rising only because the fraction of windows with a
truncated flanking region shrinks (the median window has 143 neighbours within
the default `--max-extend-bp 100000`, and 140 by 4,000 windows). Taking the
saturated 1.0 s/window gives ≈6.4 h of CPU for the full contig, i.e. about
200× the 118 s measured here.

## Commits

Every commit on the branch, oldest first, except the docs commit that last
revised this file.

| commit | change |
|---|---|
| `94158f0` | tests: add golden-file regression suite covering both lassip stages |
| `8d8eaaf` | build: track header dependencies, split arch flags, add make check |
| `6539db8` | make runs on data with missing genotypes reproducible: add --seed |
| `2aa0f19` | remove unreachable code (1357 lines) |
| `fecd9d9` | salti: hoist log/exp out of the (A, m, eps) grid search |
| `fc93ed5` | salti: one sweep-spectrum table per contig instead of one per window |
| `380d7f8` | hfs_window: skip the clustering pass when it cannot merge anything |
| `03cda7c` | readHaplotypeDataVCF: parse genotypes without per-field map lookups |
| `dd18692` | filter loci in one pass instead of two matrix copies |
| `15bc273` | hfs_window: drop the leftover per-haplotype debug write |
| `faeb0c6` | cli: help on stdout and exit 0, --version, usage line, distinct exit codes |
| `18da792` | cli: group the help by topic instead of alphabetically |
| `f05ea0b` | cli: ask param_t whether a flag was supplied instead of comparing sentinels |
| `50ecaa2` | cli: validate --dist-type at both stages, fix its column label, use the file's K |
| `1954260` | docs: refresh README and fix the example script |
| `cae6da1` | docs: branch summary with the before/after benchmark |
| `46a1e42` | match-tol: group at <= MATCH_TOL differences, as the flag documents |
| `edeae4f` | example: regenerate the committed outputs and fix the second command |
| `f3a5988` | docs: fold the match-tol and example fixes into the branch summary |
| `d004ac4` | store genotypes two bits per locus instead of one char |
| `03e8d51` | read the VCF once instead of twice |
| `a2f7bed` | hfs_window: tally packed windows instead of building char strings |
| `92cb365` | docs: record the representation, reader and HFS work |
| `da5b22c` | threading: std::thread, dynamic work claiming, and fix the null-window race |
| `c02ccfc` | docs: record the threading work |
| `11e95aa` | split main into registration, config, validation and the two stages |
| `638aab8` | one struct per population instead of seven parallel maps |
| `2a3a5db` | docs: record the structural work |
| `9af93b2` | docs: list every commit in the branch summary table |
| `f041c2f` | GMapData: bracket by binary search, and stop ignoring the placement failure |
| `2c741e8` | docs: record the genetic-map fixes |
| `70c1c7b` | add --max-gap to control how far --dist-type cm interpolates |
| `aaeecfa` | docs: record --max-gap |
| `1f8a38e` | add --hap-cluster: best-comp (new default), garud-shuffle, soft-em |

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
7. **The spectra header's window count is correct under `--threads > 1`.**
   It is windows minus null windows, and the null count was incremented from
   every thread without synchronisation, so any run with a null window and
   more than one thread wrote a wrong, non-reproducible number -- which the
   finalize stage then reads. See `tests/run_tests.sh::nullwin_threads`.
8. **Exit codes**: 64 for a command line or validation error, 65 for bad input
   data, 70 for an unexpected exception; `--help` and `--version` exit 0. An
   error raised inside a reader previously escaped `main` and aborted the
   process.
9. **Haplotype clustering is selectable and the default changed** to
   `--hap-cluster best-comp`. `garud-shuffle` reproduces the old behaviour;
   `soft-em` is experimental and makes class sizes fractional. Output on data
   without missing genotypes is unchanged at `--match-tol 0` (verified
   md5-identical across all three on the chr22 example), and the two
   missing-data regression goldens are pinned to `garud-shuffle` so the old
   path stays covered. On chr22 with 2% of genotypes blanked, `best-comp` runs
   in 11.4 s against `garud-shuffle`'s 107.8 s at one thread.
10. **A genetic map that cannot place a window is an error.** Previously
   `--dist-type cm` wrote the previous window's genetic position for any
   window in a map gap (3,951 windows sharing two positions in a test with a
   6 Mb hole) and segfaulted outright when the map named a contig the spectra
   do not use. Both now exit 65 naming the count and the first position. The
   gap threshold, hardcoded at 3 Mb, is now `--max-gap` (same default; 0
   interpolates across any gap).
11. **A null spectrum too flat to grid-search exits 65, not 64.** It is a
   property of the input, like the `checkNull` failure beside it.
12. **`--match-tol` groups at `<=` the given number of differences**, as its
   help text says, instead of `<`. The new `--match-tol t` reproduces the old
   `t+1`; `--match-tol 0` on data without missing genotypes is unchanged.

## Left for you to decide

- **Whether `soft-em` should stop being experimental.** It is the most accurate
  rule tested at high missingness and the best behaved for H12, but class sizes
  come out fractional, which changes what a spectrum is throughout the code
  (`HaplotypeFrequencySpectrum::sortedCount` is now `double`, so the machinery
  allows it, but the null-spectrum averaging and the K truncation were designed
  around counts).

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

Nothing remains from the original review list. Items noticed along the way and
not acted on:

- I/O errors are still not distinguished from data errors (both exit 65)
  because roughly 30 sites in the data layer throw untyped ints. `EXIT_IOERR`
  is defined for when they are typed.

Done since the first version of this file: genotypes packed two bits per locus
(`d004ac4`), the VCF read once rather than twice (`03e8d51`), `hfs_window`
tallying packed windows (`a2f7bed`), the threading rework (`da5b22c`), `main`
split into named stages (`11e95aa`), and per-population results in one struct
rather than seven parallel maps (`638aab8`).

One tradeoff to know about: reading the file once costs about 22 MB more peak
memory than the two-pass reader did on this dataset, because the growable copy
and the locus metadata coexist briefly. That transient is bounded by the size
of the packed matrix rather than by a second decompression pass, and `03e8d51`
is the single commit to revert if peak memory ever matters more than read time.

One review claim did not survive measurement: the static stride schedule was
said to leave cores idle at the tail. It does not. saltiLASSI parallel
efficiency on the YRI contig is 99% at 2 and 4 threads and 85% at 8 under
either schedule, the 8-thread figure being this machine's efficiency cores.
The dynamic schedule is kept because it costs nothing measurable and bounds
the tail when window costs are skewed, but it is not why `da5b22c` exists --
the data race, the undefined-behaviour cast and the leak are.
