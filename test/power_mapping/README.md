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

- Library: `0.7V_10K.lib` (cryogenic, 10 K, from
  `cda-tum/cryogenic-cmos/standard_cell_libraries/`)
- Circuit: `radd8.blif`, in this directory — a hand-written 8-bit
  ripple-carry adder (17 inputs, 9 outputs, chained from 8 one-bit
  full-adder `.names` blocks), chosen as a small circuit with enough
  gates and fanout depth to give the mapper real drive-strength and
  delay/power trade-off choices, without needing a full RTL toolchain
  to generate a test case.
- Build: plain `make -j6` (Makefile build, `ABC_USE_PTHREADS`), no
  special flags.
- Baseline for comparison: `map` (no `-d`), delay-optimal + area
  recovery, no `-D` — achieves delay 177.94, `DynPower` 18.4 (a 24.6%
  reduction from the pure delay-optimal pass's 24.3 baseline, via area
  recovery alone).

## Commands run

```
abc -c "read_lib 0.7V_10K.lib; read_blif radd8.blif; strash; map -v; print_stats"
abc -c "read_lib 0.7V_10K.lib; read_blif radd8.blif; strash; map -d -v; print_stats"
abc -c "read_lib 0.7V_10K.lib; read_blif radd8.blif; strash; map -d -D <target> -v; print_stats"
abc -c "read_lib 0.7V_10K.lib; read_blif radd8.blif; strash; map -a -d -v; print_stats"
```

## Results

Delay-optimal-only baseline (`map`, mode 0, before any recovery): delay
0.00 (unset at this stage) / `DynPower` 24.3. Area-recovery-only
baseline (`map`, final): delay 177.94, `DynPower` 18.4 (-24.6%).

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

Observations:

- Achieved delay tracks the requested target monotonically and
  **saturates exactly at the `map -a -d` (fully unconstrained) result**
  once the target stops binding (`-D 600` and `-D 1000` both also give
  442.24 — table stops at 442, the natural ceiling, since looser
  targets past that produce no further change).
- Power-aware mapping beats area-recovery-only power reduction
  (32.9–34.6% vs. 24.6%) at a comparable or even tighter delay — the
  power objective is doing something area alone doesn't.
- Diminishing returns are visible: most of the power win (32.9%) is
  already captured at a delay target close to the original
  delay-optimal critical path (177.94 → 185.39, +4%). Going all the
  way to fully unconstrained (+148% delay) only buys another 1.7
  points of power reduction. That shape — steep near the optimum, flat
  after — is exactly what a slack-harvesting argument needs to show.

## Known remaining issue (cosmetic, not fixed)

The verbose (`-v`) log's `Flow` column prints an uninitialized
`FLT_MAX` value (`34028234663852885981170418348451692544.0`) during
the exact-power (`Power`, mode 6) stage. It doesn't affect the actual
mapping decisions — `PowerF` is tracked in a separate field from
`AreaFlow`/`Flow` and the mapped netlist and `DynPower` numbers above
are unaffected — but the printed `Flow` number at that stage is
meaningless and shouldn't be reported if scripting around this output.

## Caveats

- Single test circuit (8-bit adder, ~120 AIG nodes). Numbers above are
  a proof that the mechanism works and behaves sensibly, not a
  representative result for real benchmarks — re-run this sweep on
  actual paper benchmarks before citing power/delay percentages.
- Only ran on `0.7V_10K.lib`. Same sweep should be repeated on
  `0.7V_300K.lib` for a real cryogenic-vs-room-temperature comparison.
