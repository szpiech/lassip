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
| `cbb058b` | docs: list the --hap-cluster commit in the branch summary |
| `59723c1` | portability: include <cstring>, use fabs, initialise MAX_EXTEND |
| `1d23854` | docs: record the portability fixes |
| `cd72a4b` | help: remove a method and citation that do not exist |
| `5372bc7` | docs: record the citation correction |
| `d00eb19` | bump the version to 1.3.0 |
| `07c3471` | docs: record the version bump |
| `fc3cf8b` | linux v1.3.0 |
| `4eb0255` | macos-arm v1.3.0 |
| `d6f8942` | tests: cover --unphased with missing genotypes, and stage 2 from .mlg spectra |
| `9d7d080` | docs: record the unphased coverage work |
| `8355dd8` | tests: make the phased/unphased coverage matrix symmetric |
| `d9d8090` | docs: record the symmetric coverage matrix |
| `2953ca4` | macos-arm binary |
| `b31494f` | readHaplotypeDataVCF: skip leading whitespace before the CHROM token |
| `962ed12` | docs: record the leading-whitespace fix |
| `add1534` | stage 2: delimit contigs by the chr column, not by the input file |

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
13. **Two reader errors say more.** The allele-coding error now names the
   genotype it read, the sample and the position, because that is where a
   malformed or shifted record surfaces and the bare message gave nothing to
   look at; a record with no CHROM field is reported instead of being parsed
   onward. Both exit 65 like the other input-data errors.

## Multiple contigs

`add1534` is Phase 0 of multi-contig support (see `multi-contig-proposal.md`):
stage 2 delimits contigs by the **chr column** rather than by the input file, so
`--spectra` takes one multi-contig file or N single-contig files
interchangeably, and the two give byte-identical results.

It is also a live bug fix. `readSpecData` built one `SpectrumData` per file, so
a file holding several contigs became one block, and `calcMTA` bounds its
flanking scan only by `SpectrumData::nwins`. Measured on two contigs sharing a
coordinate range (188 windows, boundary at row 94), merged file vs the same two
files separately:

| | effect before the fix |
|---|---|
| `--dist-type nw` | contig 2's entire `pos` column wrong (94 rows), since stage 2 writes `dist[w] = w`, the index within the block; plus `m`/`A`/`L` differing in 4/5/10 rows inside the 5-window reach of the boundary, `L` by up to 82% |
| `--dist-type bp` | `m`/`A`/`L` differ across rows 79-95, a contiguous run spanning the boundary's reach; `A` moves by a factor of 38 (6.9e6 -> 1.8e5) |
| `--dist-type cm` | likewise |

The null spectrum is byte-identical either way, so nothing global is involved --
but the affected windows are every contig end.

The boundary is now the edge of a block rather than a bounds check, which makes
`--dist-type nw` correct by construction since its distance is an index into
the block. A first pass reads only the `chr` and `start` columns to learn the
layout, and rejects a file that returns to a contig it has left, is not
ascending within a contig, disagrees with its own window count, or repeats a
contig another file supplied. The header may carry
`contigs <n> <name> <nwins> ...` after the population names (D3), validated
when present and appended rather than inserted so the file stays readable by
1.2.x. Nothing writes it yet -- stage 1 gains that in Phase 1.

Three cases added; suite is 48 checks. Against the pre-fix binary all three
`--salti` arms fail, as do both validation cases. The `--lassi` and
`--avg-spec` arms pass pre-fix and are coverage rather than bug-pinning:
LASSI uses no flanking windows, and the averaging is order-independent.

## Test coverage

`d3b1a61` closes a coverage hole I should have noticed earlier: only 1 of the
22 cases passed `--unphased`, and it used the fixture with no missing
genotypes, so the clustering code had never run on multilocus genotype
strings -- the one path whose alphabet is `{0,1,2,-}` and therefore the only
one where the packed representation uses all four symbols.

The covered part was verified unchanged first: with no missing data,
`--unphased` output is byte-identical to v1.2.2 on both `testing/small.vcf.gz`
and `example/YRI.chr22.vcf.gz` (4750 hom-alt genotypes, so all three symbols
appear). Three cases added -- `spec_unphased_missing` (legacy rule pinned on
the `.mlg` path), `cluster_unphased` (default rule, `--seed` invariance, and
count conservation across all three rules) and `salti_unphased` (stage 2 from
a `#phased 0` spectra file). Suite is 28 checks.

`8355dd8` then filled the rest of the unphased row, so the matrix is
symmetric: every combination of {1,2} populations x {clean, missing} x
{stage 1, stage 2} x {filter 0,1,2} tested for phased data is now tested for
unphased too, and vice versa -- `cluster_twopop` had to be added because the
unphased row reached a combination (two populations with missing genotypes)
that the phased row never had.

| combination | phased | unphased |
|---|---|---|
| 1 pop, clean, stage 1, filter 0 | yes | yes |
| 1 pop, clean, stage 1, filter 2 | yes | yes |
| 1 pop, clean, stage 2, filter 2 | yes | yes |
| 1 pop, missing, stage 1, filter 2 | yes | yes |
| 2 pop, clean, stage 1, filter 1 | yes | yes |
| 2 pop, clean, stage 1, filter 2 | yes | yes |
| 2 pop, missing, stage 1, filter 2 | yes | yes |

34 cases, 40 checks; 23 phased, 11 unphased. The useful result: all nine
checks of the clean unphased cases pass against the v1.2.2 binary, so
`--unphased` output is byte-identical to 1.2.2 at every filter level, for one
and two populations, and through `--lassi`, `--avg-spec`, `--null-spec` and
`--salti`. That is what validates the packed genotype representation, the
single-pass reader and the single-pass filter on the `{0,1,2,-}` path.

Against the v1.2.2 binary the four missing-data cluster cases fail --
`spec_unphased_missing`, `cluster_unphased`, `cluster_unphased_twopop` and
`cluster_twopop` -- since `--hap-cluster` does not exist there, so they pin
genuinely new behaviour. Every clean case passes against it, including
`salti_unphased`, `lassi_unphased`, `avg_spec_unphased` and
`lassi_nullspec_unphased`: they use only flags v1.2.2 has, and their input
spectra are byte-identical between the two builds. Those close coverage gaps
rather than pinning changes.

`b31494f` fixes a regression the suite could not have caught, and adds the
case that would have. Every VCF in `example/` and `testing/` is tab-delimited
with no leading whitespace, so no fixture exercised the one assumption the
rewritten reader in `03cda7c` added: that a record's first character begins
the CHROM field. Given a record written as `" 1<TAB>417<TAB>..."`, the reader
took the empty string before the space as the contig name and shifted every
later field by one, so the FORMAT column was read as the first sample's
genotype and the run died on "Alleles must be coded 0/1/. only" -- an error
naming the genotypes, which were fine. `vcf_leading_space` is
`testing/small.vcf.gz` with a space prepended to every data line, run with the
`spec_phased` flags and compared against the `spec_phased` golden, so it
asserts that leading whitespace changes nothing rather than only that the run
survives. Suite is 34 cases, 40 checks. The general lesson is that fixtures
derived from repository data share the repository's formatting, so a reader
rewrite wants at least one fixture that is deliberately formatted differently.

The exercise also produced a result worth knowing: the `--hap-cluster` choice
matters much less for unphased data. Same VCF, same window, `--match-tol 0` --
phased: 149 distinct strings, 72 of them compatible with more than one other,
rules return 83/82/82 classes. Unphased: 81 distinct strings, 3 ambiguous,
all three rules return 74.

## Version

`d00eb19` sets `VERSION` to `1.3.0` (it was still `1.2.2`, so `--version` and
the startup banner reported the version the branch was built from). The README
changelog entry is headed `15SEP2026 - v1.3.0.` following the convention of the
entries above it -- change the date when you tag, nothing reads it. No output
file carries a version stamp, so no golden changed.

Two things outside `src/` still refer to older versions and are yours to decide
on at release: `bin/` holds prebuilt binaries up to `lassip-v1.2.1` (and no
1.2.2), and the manual PDFs in the repository root are `v1.1.2`-era, so they
predate `--seed`, `--hap-cluster`, `--max-gap`, `--version` and the grouped
help.

## Correction

`cd72a4b` removes a method and a citation from the help text that I had
invented: the PREAMBLE added in `faeb0c6` listed "2-population LASSI (Harris
and DeGiorgio 2020, Genetics 210:1429)". No such method is implemented and no
such paper exists -- it came from misreading `calc_LASSI_stats2` as a
two-population variant (the `2` is a stage suffix) and garbling the G123
citation already in the README. The remaining citations now follow the README's
own forms, and the fourth slot names a real thing: `--unphased` writes the two
haplotype-statistic columns as `g123`/`g2g1`, which is Harris et al. 2018.

If you are reviewing anything else I added to user-facing text, the citations
are the part to check first -- everything else in the branch is verifiable
against the code or the test suite, but a reference is not.

## Portability

`59723c1` fixes a build break on Linux/libstdc++: `lassip-data.cpp` used
`memset`/`memcpy` without `<cstring>`, which libc++ supplies transitively and
libstdc++ does not. Introduced by `03e8d51` on this branch. The same commit
hardens two things that predate the branch and are latent rather than fatal:
four unqualified `abs()` calls on doubles became `fabs()` (an integer overload
would truncate `--dist-type cm` distances silently), and `MAX_EXTEND` in
`lassip-wintools.cpp` is now initialised with an explicit else rather than left
indeterminate if `--dist-type` ever admits a fourth value.

The tree compiles under GCC 16 with `-Wall -Wextra` as well as clang; the two
latent items were found by GCC's diagnostics. Remaining warnings are
pre-existing unused parameters (`PHASED` in `filterHaplotypeData`, `DIST_TYPE`
in `initResults`, `len` in the three clustering routines).

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
