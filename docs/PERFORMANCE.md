# Performance evidence and provenance

These measurements belong to the frozen native candidate before the namespace
was changed to `rivin_jones`. Packaging is verified by fresh coefficient and
invariant tests; it is not a new performance sweep.

The 25 retained and ten unseen complete-polynomial cases all beat Sage 10.9
on macOS arm64 / Python 3.12.16, with speedups from 5.61x to 3,364x. All 29,005
calibration/timed computations matched full integer coefficient maps. Imports
were excluded, braid inputs were fresh, and engine order was randomized in
five batches. There was no cross-call answer cache.

The separate T(2,7), color-7 observation took 67.361 seconds with the preceding
cable engine, 0.001701 seconds with the new dispatcher, and 0.5915 seconds with
Sage. This was a single call per engine and is not part of the repeated sweep.

Five completed direct-root comparisons beat Sage's polynomial-plus-evaluation
route. The sixth (weaving [1,-2]*4, color 6) exhausted its 180-second comparison
worker budget. A separate native call completed in 0.0116 seconds and agreed
with native polynomial evaluation; a fresh Sage call missed a 30-second
process deadline. That case remains unranked and lacks a completed independent
Sage comparison at that color.

Figure-eight root checks at color 200 took 0.352–0.377 seconds at 864/992 bits.
These were validation observations, not a repeated speed comparison. The
intervals enclose finite-N values; they do not bound error to the conjectured
limiting hyperbolic volume. Trace rows still grow as N^(strands-1).

The [retained](evidence/benchmark-retained.json),
[unseen](evidence/benchmark-heldout.json), [single-call](evidence/benchmark-hard.json)
and [root](evidence/benchmark-roots.json) reports preserve inputs, timings,
work counts and failures. Machine paths and environment inventories were
removed from this export. [Provenance](evidence/provenance.json) records the
source manifest and fixture hashes. The [405-group direct audit](evidence/audit-release.json)
and [330-group colored regression](evidence/colored-regression.json) are retained
as pre-packaging evidence, including all 83 ordinary outputs.

These are finite-corpus results. A new platform, compiler, dependency stack,
namespace/interface change or algorithm revision needs its own appropriate
validation; the old timing record is not a blanket portability or speed claim.
