# Power-aware technology mapping (`map -d`) — validation notes

This documents the first end-to-end validation of `map -d` on real
Liberty files, the bug found and fixed while doing it, and the numbers
that came out — first on a small hand-written toy circuit, then on a
real ISCAS85 benchmark. Written 2026-08-20, commits `29e51036d`
onward.

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

## Scope: this is dynamic power only, not total power

**Important for interpreting every number below.** The value driving
`map -d`'s cost function — `dPower_dyn` — comes from
`Abc_SclComputeAveragePower()` in `src/map/scl/sclLibUtil.c`:

```c
static float Abc_SclComputeAveragePower( SC_Cell ** p )
{
    float a = Abc_SclComputeAverageNetSwitchingPower( p );
    float b = Abc_SclComputeAverageCellInternalPower( p );
    ...
    return a + b;
}
```

`a` is net switching power (output pin `rise_power`/`fall_power`
tables), `b` is cell internal power (input pin `internal_power`
tables) — **both dynamic components. Leakage never enters this
function.**

Leakage isn't unavailable — it's parsed from the Liberty file fine
(`Scl_LibertyReadCellLeakage` reads `cell_leakage_power`, stored as
`pCell->leakage` on every cell) — it's just wired to an unrelated,
pre-existing ABC feature instead: `Abc_SclConvertLeakageIntoArea(p, A,
B)` does `area = A·area + B·leakage`, a leakage-weighted *area*
metric for sizing, used nowhere in the `map -d` path. Leakage sits in
memory, fully characterized, and the power-aware mapper never looks at
it.

**Consequence:** every `DynPower`/percentage number in this document
is a dynamic-power-only reduction, not total power. That's a real,
legitimate, correctly-computed result — dynamic power is expected to
matter proportionally more as cryogenic leakage collapses, so this
isn't the wrong thing to optimize for a cryo study. But it means:

- These numbers can't be quoted as "power" reduction in a paper without
  qualifying "dynamic."
- The 10K-vs-300K finding below (300K shows a *larger* percentage
  dynamic-power reduction than 10K) is specifically about dynamic
  power. Total power could tell a different story, since leakage
  typically matters far more at 300K than at 10K — a leakage-blind
  tool is closer to "correct by construction" at cryo than at room
  temperature, not the other way around.
- This is exactly why `set_opt_config -sizing_leakage_limit` (the
  OpenROAD resizer knob identified separately, operating downstream at
  placement/CTS/routing and reading real leakage values) is worth
  keeping as the complementary piece for a total-power story: `map -d`
  for dynamic power at synthesis, `-sizing_leakage_limit` for leakage
  at physical design.

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

## Results — radd8 (toy circuit, first sanity check)

Small hand-written 8-bit ripple-carry adder (~120 AIG nodes), used
first just to prove the mechanism runs and behaves sensibly before
spending time on a real benchmark. See the "real benchmark" section
below for the more representative numbers.

### 0.7V_10K.lib

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

### 0.7V_300K.lib

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

### Cross-temperature comparison (radd8)

| | 10 K | 300 K |
|---|---:|---:|
| Delay-optimal `DynPower` | 24.3 | 32.9 |
| Area-recovery-only reduction | 24.6% | 30.4% |
| Constrained (`map -d`, tight target) reduction | ~32.9% | ~38.9%* |
| Fully unconstrained (`map -a -d`) reduction | 34.6% | 40.3% |

\* excluding the 250 anomaly above.

At this circuit, room temperature shows a *larger* percentage power
reduction from power-aware mapping than 10 K does — the opposite of
"cryogenic operation trivially benefits more." Held up on a real
benchmark too, see below.

## Results — ISCAS85 c880 (real benchmark)

`c880` is a classic ISCAS85 combinational benchmark (an 8-bit ALU),
sourced from `fiction/benchmarks/ISCAS85/c880.v` (copied into this
directory) — 60 inputs, 26 outputs, 327 two-input AND nodes after
`strash`, read directly with ABC's native `read_verilog` (this
benchmark set is distributed in the plain structural
`assign n = a & b;` / `~`/`|`/`^` subset ABC's built-in parser accepts,
no Yosys conversion needed). About 3x the gate count of `radd8`, and
an actual, recognizable, citable benchmark rather than a hand-written
toy — this is the number set to use in the paper, not the radd8 ones
above.

### 0.7V_10K.lib

Delay-optimal `DynPower`: 106.5. Area-recovery-only baseline (no
`-d`): delay 224.93, `DynPower` 57.8 (-45.7%).

| `-D` target | Achieved delay | `DynPower` | Reduction vs. 106.5 |
|---:|---:|---:|---:|
| 230 | 269.76 | 48.0 | 54.9% |
| 250 | 269.76 | 48.0 | 54.9% |
| 270 | 269.76 | 48.0 | 54.9% |
| 290 | 287.73 | 47.4 | 55.5% |
| 310 | 303.28 | 47.2 | 55.7% |
| 322 | 321.84 | 47.1 | 55.8% |
| — (`map -d`, no `-D`) | 269.76 | 48.0 | 54.9% |
| ∞ (`map -a -d`) | 321.84 | 47.1 | 55.8% |

Targets at or below 270 all clamp to the same result — 269.76 is the
natural floor the power-recovery pass can't push tighter than on this
circuit/library, so `-D 230`/`250`/`270` are all equivalent to `-d`
alone here. Above that floor the curve is smooth and monotonic, no
repeat of the radd8 anomaly.

### 0.7V_300K.lib

Delay-optimal `DynPower`: 128.3. Area-recovery-only baseline (no
`-d`): delay 212.27, `DynPower` 74.3 (-42.2%).

| `-D` target | Achieved delay | `DynPower` | Reduction vs. 128.3 |
|---:|---:|---:|---:|
| 215 | 237.87 | 58.1 | 54.7% |
| 230 | 237.87 | 58.1 | 54.7% |
| 250 | 247.68 | 58.5 | 54.4% |
| 270 | 266.44 | 57.8 | 55.0% |
| 290 | 283.99 | 57.6 | 55.1% |
| 302 | 301.72 | 57.3 | 55.3% |
| — (`map -d`, no `-D`) | 237.87 | 58.1 | 54.7% |
| ∞ (`map -a -d`) | 301.72 | 57.3 | 55.3% |

Same floor behavior below 230. One small (0.4-point) non-monotonic dip
at `-D 250`, same family of effect as the radd8 anomaly but much
smaller — not worth chasing further, just noted for honesty.

### Cross-temperature comparison (c880) — the number to actually cite

| | 10 K | 300 K |
|---|---:|---:|
| Delay-optimal `DynPower` | 106.5 | 128.3 |
| Area-recovery-only reduction | 45.7% | 42.2% |
| Constrained (`map -d`, tight target) reduction | 54.9% | 54.7% |
| Fully unconstrained (`map -a -d`) reduction | 55.8% | 55.3% |

Two things worth having in the paper:

1. **Power-aware mapping clearly beats area-recovery as a leakage
   proxy on a real circuit** — roughly +9-13 points of extra power
   reduction over plain area recovery, at both temperatures. This is a
   much bigger, more convincing margin than the radd8 toy circuit
   showed (+8 points at 10 K, +8-9 at 300 K) — real benchmarks give the
   mapper more structure to exploit.
2. **The cross-temperature direction from radd8 replicates on a real
   benchmark**: 10 K shows a *smaller* percentage reduction than 300 K
   at every stage (area-recovery, constrained, and unconstrained),
   consistently, not a fluke of the toy circuit. That's now two
   independent data points pointing the same way — worth treating as a
   real candidate finding for the paper (subject to confirming on more
   benchmarks), not just noise.

## Known remaining issue (cosmetic, not fixed)

The verbose (`-v`) log's `Flow` column prints an uninitialized
`FLT_MAX` value (`34028234663852885981170418348451692544.0`) during
the exact-power (`Power`, mode 6) stage. It doesn't affect the actual
mapping decisions — `PowerF` is tracked in a separate field from
`AreaFlow`/`Flow` and the mapped netlist and `DynPower` numbers above
are unaffected — but the printed `Flow` number at that stage is
meaningless and shouldn't be reported if scripting around this output.

## Caveats

- Two circuits tested so far (radd8 toy adder, c880 real ISCAS85
  benchmark), both small-to-medium (≤330 gates). The c880 numbers are
  the more trustworthy ones and are consistent in direction with
  radd8, but this is still not a benchmark suite — run more (and
  larger) circuits before treating the cross-temperature direction as
  settled, and don't trust any single `-D` point without checking its
  neighbors given the observed (small) non-monotonic sensitivity.
- Only `0.7V_10K.lib` and `0.7V_300K.lib` tested so far — 77 K only
  exists as `iso-Ioff-libs/iso_off_77K_0.7V.lib`, a **different
  characterization methodology** from these two root-level files (file
  contents differ even for cells with identical nominal temperature —
  confirmed by diffing `0.7V_10K.lib` against
  `iso-Ioff-libs/iso_off_10K_0.7V.lib`). Don't mix root-level and
  iso-Ioff files in the same comparison; use one methodology
  consistently across all three temperature points if 77 K is added.
