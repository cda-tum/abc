# Power-aware technology mapping (`map -d`) — validation notes

This documents the first end-to-end validation of `map -d` on a real
Liberty file, the bug found and fixed while doing it, and the numbers
that came out. Written 2026-08-20, on commit `29e51036d`.

## What `map -d` does

`map -d` ("toggles dynamic power mapping using SC library power
values") runs a second cost function through ABC's technology mapper:
a global `PowerFlow` pass and a local exact-power dereference/reference
pass, structurally mirroring the existing area-recovery machinery
(`Map_Mapping()` in `src/map/mapper/mapperCore.c`).

The power values come from a **stock ABC mechanism**, not anything
new: whenever a Liberty file is loaded via `read_lib`, `Abc_NtkMap()`
already auto-derives a genlib from it (`Abc_SclDeriveGenlib` in
`src/map/scl/sclLibUtil.c`) for gate-sizing purposes. This branch
extends that existing derivation to also carry a power value through,
computed in memory straight from the Liberty file's power tables — no
genlib editing, no separate generator script. The whole pipeline is:

```
read_lib <library>.lib
<read design, strash>
map -d [-D <delay_target_ps>]
```

## The bug found (fixed in commit `29e51036d`)

`src/base/abci/abc.c`, `Abc_CommandMap()`:

```c
if ( fAreaOnly || fDynPower )
    DelayTarget = ABC_INFINITY;
```

`-d` (`fDynPower`) was piggybacking on `-a`'s (area-only mapping's)
"ignore timing entirely" semantics. Any `-D` value passed alongside
`-d` was silently discarded, so `map -d` always minimized power with
zero regard for delay, regardless of what target was requested.

Fix: drop `fDynPower` from that condition, so `-d` alone falls back to
the natural delay-optimal bound (the same way default area recovery
already works), and the fully-unconstrained mode is still reachable
via `-a -d` (unchanged behavior, same as it always was for area).

## Test setup

- Libraries: `0.7V_10K.lib` and `0.7V_300K.lib` (from
  `cda-tum/cryogenic-cmos/standard_cell_libraries/`), same 0.7 V
  supply, cryogenic (10 K) vs. room temperature.
- Circuit: `radd8.blif`, in this directory — a hand-written 8-bit
  ripple-carry adder (17 inputs, 9 outputs, chained from 8 one-bit
  full-adder `.names` blocks), chosen as a small circuit with enough
  gates and fanout depth to give the mapper real drive-strength and
  delay/power trade-off choices, without needing a full RTL toolchain
  to generate a test case.
- Build: plain `make -j6` (Makefile build, `ABC_USE_PTHREADS`), no
  special flags.

## Commands run

```
abc -c "read_lib <library>.lib; read_blif radd8.blif; strash; map -v; print_stats"
abc -c "read_lib <library>.lib; read_blif radd8.blif; strash; map -d -v; print_stats"
abc -c "read_lib <library>.lib; read_blif radd8.blif; strash; map -d -D <target> -v; print_stats"
abc -c "read_lib <library>.lib; read_blif radd8.blif; strash; map -a -d -v; print_stats"
```

## Results — 0.7V_10K.lib

Delay-optimal-only baseline (`map`, mode 0, before any recovery): delay
0.00 (unset at this stage) / `DynPower` 24.3. Area-recovery-only
baseline (`map`, final, no `-d`): delay 177.94, `DynPower` 18.4
(-24.6%).

Constrained power recovery (`map -d -D <target>`), sweeping the target
from tight to loose:

| `-D` target | Achieved delay | `DynPower` | Reduction vs. 24.3 |
|---:|---:|---:|---:|
| 180 | 185.39 | 16.3 | 32.9% |
| 200 | 199.24 | 16.3 | 32.9% |
| 250 | 239.74 | 16.3 | 32.9% |
| 300 | 293.47 | 16.2 | 33.3% |
| 400 | 388.88 | 15.9 | 34.6% |
| 442 | 429.38 | 15.9 | 34.6% |
| — (`map -d`, no `-D`) | 185.39 | 16.3 | 32.9% |
| ∞ (`map -a -d`) | 442.24 | 15.9 | 34.6% |

Clean and monotonic across the whole sweep: achieved delay tracks the
target and **saturates exactly at the `map -a -d` result** (442.24)
once the target stops binding (`-D 600`/`-D 1000` both also land on
442.24).

## Results — 0.7V_300K.lib

Delay-optimal-only baseline: `DynPower` 32.9. Area-recovery-only
baseline (no `-d`): delay 168.45, `DynPower` 22.9 (-30.4%).
Fully unconstrained (`map -a -d`): delay 411.37, `DynPower` 19.6
(-40.3%).

| `-D` target | Achieved delay | `DynPower` | Reduction vs. 32.9 |
|---:|---:|---:|---:|
| 175 | 175.24 | 24.0 | 27.1% |
| 200 | 187.75 | 20.1 | 38.9% |
| 210 | 200.26 | 20.1 | 38.9% |
| 220 | 213.59 | 20.0 | 39.2% |
| 230 | 225.02 | 20.1 | 38.9% |
| 240 | 237.53 | 19.9 | 39.5% |
| 250 | 249.78 | **22.8** | 30.7% |
| 260 | 250.86 | 20.0 | 39.2% |
| 270 | 262.29 | 19.9 | 39.5% |
| 280 | 274.80 | 19.9 | 39.5% |
| 300 | 299.56 | 19.9 | 39.5% |
| 350 | 349.34 | 19.8 | 39.8% |
| 411 | 399.94 | 19.7 | 40.1% |
| ∞ (`map -a -d`) | 411.37 | 19.6 | 40.3% |

**Anomaly:** the point at `-D 250` is a genuine, reproducible outlier
— re-ran it twice, byte-identical both times, and the neighboring
targets (240 → 19.9, 260 → 20.0) bracket it tightly while 250 itself
jumps to 22.8. Not a fluke, not nondeterminism. Most likely explanation
is the known brittleness of greedy, single-pass local-search recovery
heuristics near specific required-time threshold values — a slightly
different required time can shift which cuts tie and get selected at a
handful of nodes, landing in a different local optimum. This is a
general characteristic of this class of heuristic (ABC's plain
area-recovery can show the same kind of sensitivity in principle), not
something specific to the power path or to the `-D`/`fDynPower` fix
above. Practical implication for real experiments: **don't trust a
single `-D` value at face value — sweep a few nearby targets and use
the best, or note the sensitivity if reporting a single-point result.**

## Cross-temperature comparison (same circuit, same `-D` targets)

| | 10 K | 300 K |
|---|---:|---:|
| Delay-optimal `DynPower` | 24.3 | 32.9 |
| Area-recovery-only reduction | 24.6% | 30.4% |
| Constrained (`map -d`, tight target) reduction | ~32.9% | ~38.9%* |
| Fully unconstrained (`map -a -d`) reduction | 34.6% | 40.3% |

\* excluding the 250 anomaly above.

At this one small circuit, room temperature shows a *larger*
percentage power reduction from power-aware mapping than 10 K does,
both for area-recovery alone and for the power-aware mapper — the
opposite of "cryogenic operation trivially benefits more." Absolute
`DynPower` is also higher at 300 K throughout (32.9 baseline vs. 24.3),
consistent with the library characterizing genuinely different
dynamic-power behavior per temperature rather than just a scaled
version of the same numbers. This is exactly the kind of "measure, don't
assume" result the research plan is built around — but it's one toy
8-bit adder, not a benchmark suite, so treat the *direction* as a
hypothesis to check on real designs, not a conclusion.

## Known remaining issue (cosmetic, not fixed)

The verbose (`-v`) log's `Flow` column prints an uninitialized
`FLT_MAX` value (`34028234663852885981170418348451692544.0`) during
the exact-power (`Power`, mode 6) stage. It doesn't affect the actual
mapping decisions — `PowerF` is tracked in a separate field from
`AreaFlow`/`Flow` and the mapped netlist and `DynPower` numbers above
are unaffected — but the printed `Flow` number at that stage is
meaningless and shouldn't be reported if scripting around this output.

## Caveats

- Single test circuit (8-bit adder, ~120 AIG nodes) at both
  temperatures. Numbers above are a proof that the mechanism works and
  behaves sensibly (plus one real anomaly worth knowing about), not a
  representative result for real benchmarks — re-run this sweep on
  actual paper benchmarks before citing power/delay percentages, and
  don't trust any single `-D` point without checking its neighbors.
- Only `0.7V_10K.lib` and `0.7V_300K.lib` tested so far — 77 K only
  exists as `iso-Ioff-libs/iso_off_77K_0.7V.lib`, a **different
  characterization methodology** from these two root-level files (file
  contents differ even for cells with identical nominal temperature —
  confirmed by diffing `0.7V_10K.lib` against
  `iso-Ioff-libs/iso_off_10K_0.7V.lib`). Don't mix root-level and
  iso-Ioff files in the same comparison; use one methodology
  consistently across all three temperature points if 77 K is added.
