# Corridor buffer scaling — spec

Date: 2026-07-04
Scope: `pipelines/channel_attrs/` corridor construction (consumed by the
fine-raster extraction tasks; see the ddrs program spec
`ddrs docs/superpowers/specs/2026-07-04-leakance-gate-program-design.md` §3).
Status: agreed design; supersedes the flat two-tier (100 m / 200 m) corridor
widths for reaches where bankfull width makes them too narrow.

## The question

How should the corridor buffer half-width scale with river size (stream
order)? Intuition says "smaller buffer for smaller streams" — the literature
says something more specific: the buffer serves two different masters, and
which one dominates flips with river size.

## Two regimes, one hinge

The corridor must cover (a) the flowline's POSITIONAL ERROR and (b) the
CHANNEL + BANKS themselves.

- MERIT flowlines carry a typical lateral positional error of 100–300 m
  (Amatulli et al. 2022 benchmarked 90 m-DEM networks against NHDPlus HR:
  <50% of stream cells within 100 m of the reference; flat valleys are the
  worst case per Yamazaki et al. 2019). This error is ORDER-INDEPENDENT.
- Channel width spans 1–2 m (order 1) to >1 km (order 9–10). It grows
  geometrically with order.

For small streams, width ≪ positional error: shrinking the buffer below
~100 m does not make extraction "more precise" — it makes the corridor MISS
the channel entirely, because the vector is off by more than the buffer.
For big rivers, width ≥ positional error: a fixed 100 m corridor samples a
strip of the channel, not the channel.

**The rule (hinge, not power law):**

```
half_width(reach) = max( E_pos , alpha * w_bankfull(reach) / 2 * 3 )
                  = max( 100 m , 1.5 * w_bankfull(reach) )

E_pos = 100 m   positional-error floor (Amatulli 2022; StreamCat precedent
                for the fixed value — Hill et al. 2016 use 100 m riparian
                buffers on NHDPlus)
alpha ~ 3       channel + banks/hyporheic margin (half-width = 1.5 * w)
```

## Literature chain for the width–order relationship

If expressed in stream order, the chain is:

1. **Leopold & Maddock (1953)** — downstream hydraulic geometry:
   `w ∝ Q^0.5`.
2. **Moody & Troutman (2002)** — the coefficient: `w ≈ 7.2 · Q^0.5`
   (w in m, Q in m³/s).
3. **Horton/Strahler laws** — drainage area (and hence Q) grows
   geometrically with order, area ratio `R_A ≈ 4–5` per order.
4. Composing 1–3: width grows ≈ **1.9–2× per order** — the ratio
   **Downing et al. (2012, Inland Waters)** confirm empirically in their
   global width-by-order distributions (order 1 ≈ 1.6–2 m, roughly doubling
   per order).

The order-form of the rule and where the hinge engages:

```
w(omega) ~ 2 m * 1.9^(omega - 1)
half_width(omega) = max(100 m, 1.5 * w(omega))

omega:        1     2     3     4     5     6     7      8
w(omega) m:   2     4     7    14    26    50    95    180
half-width:  100   100   100   100   100   100  ~140   ~270
                                                 ^ hinge engages ~order 6-7
```

The buffer is FLAT at the error floor through order ~5–6. This is why the
fixed-100 m practice in the riparian literature (StreamCat) works as well as
it does: most reaches never reach the hinge.

## Implementation rule for this pipeline

DO NOT use order as the primary scale — order is a proxy for width with
±1-order scatter (Downing's order-width distributions overlap heavily). The
pipeline already ingests two per-reach width estimates:

- `bankfull_width` — Zarrabi et al. 2025 (all NHDPlus-crosswalked reaches)
- `channel_width_obs` — SWORD/GRWL Landsat-observed (rivers ≥ 30 m)

Priority: `channel_width_obs` where present (observed beats modeled), else
`bankfull_width`, else the order fallback `w(omega) = 2 * 1.9^(omega-1)` m
using the MERIT shapefile's `order` field (crosswalk-failure reaches only).

```
half_width = max(100.0, 1.5 * width_est)     # meters, EPSG:5070
```

Deliverable: a third corridor set `corridors_scaled.parquet` built with this
rule, used by the FINE-raster extraction tasks (imperviousness-class data if
ever rasterized; alluvium overlay). The coarse-grid WTD sampling
(nearest-channel-cell, ~250 m–1 km rasters) is unaffected — no fine buffer
there by design (RiverATLAS-style native-grid association).

Consequences for the existing sets:
- `corridors_100m.parquet` — remains valid for all reaches below the hinge
  (the overwhelming majority) and as the sensitivity column.
- `corridors_wide.parquet` (flat 200 m) — superseded by `corridors_scaled`;
  the flat 200 m misses most of the channel on >1 km-wide rivers (lower
  Mississippi class) where the scaled rule gives ≥1.5 km half-widths.

## Caveats

- The order fallback uses global-average coefficients; arid-region channels
  run wider-and-shallower than the L&M exponents suggest — acceptable for a
  fallback that only fires on crosswalk-failure reaches.
- `alpha = 3` (i.e. 1.5×w half-width) is a judgment call bounding banks +
  hyporheic margin; sensitivity to alpha ∈ {2, 4} is cheap to add if a
  downstream result hinges on it.
- Where the corridor is width-scaled AND the reach is braided, the GRWL
  width already includes the braid belt (it measures the water surface at
  mean discharge), so no extra braiding factor is applied.

## References

- Leopold, L.B. & Maddock, T. (1953). The hydraulic geometry of stream
  channels and some physiographic implications. USGS Prof. Paper 252.
- Moody, J.A. & Troutman, B.M. (2002). Characterization of the spatial
  variability of channel morphology. Earth Surf. Process. Landforms 27(12).
- Downing, J.A. et al. (2012). Global abundance and size distribution of
  streams and rivers. Inland Waters 2(4), 229–236.
- Amatulli, G. et al. (2022). Hydrography90m. ESSD 14, 4525–4550.
  DOI:10.5194/essd-14-4525-2022.
- Hill, R.A. et al. (2016). The StreamCat Dataset. JAWRA 52(1), 120–128.
- Yamazaki, D. et al. (2019). MERIT Hydro. WRR 55(6), 5053–5073.
- Zarrabi, M. et al. (2025). Bankfull and mean-flow channel geometry via ML
  (CONUS). WRR 61, e2024WR037997.
- Altenau, E.H. et al. (2021). SWORD. WRR 57, e2021WR030054.
- Full verified-citation list: ddrs
  `docs/2026-07-04-leakance-literature-review.md` §6.
