# Revision analyses A1–A11 — results summary

Branch: `revision-response`  
Generated from scripts in `scripts/revision/`.  
Re-run everything with:

```bash
Rscript scripts/revision/00_run_all.R
```

Editable side tables (fill/verify before finalising the response letter):

- `data/revision/tephra_events.csv` (curated lake×event inventory + refs; regional rows excluded; `include_sensitivity` flags A6 predictor ages)
- `data/revision/fish_stocking.csv` (Azul = 1792 CE; other lakes still NA)
- `data/revision/lake_trophic_state.csv` (Caldeirão seeded; others from Table S3)

---

## A1 — Rarefied richness (R1)

**Rarefaction depths** (highest n retaining ≥70% of 30-yr bins):

| Group | n | bins kept |
|---|---:|---:|
| Producers (diatoms) | 500 | 209 / 267 (78%) |
| Consumers (chironomids) | 50 | 180 / 237 (76%) |

Chironomid counts remain low (lake medians often <<100); n = 50 is the honest common depth.

**Period means of rarefied richness** (no within-lake min–max scaling):

| Group | pre-1600 | post-1600 |
|---|---:|---:|
| Producers | 31.7 | 29.3 |
| Consumers | 6.84 | 6.92 |

**Implication for claims:** A1 answers R1’s sedimentation / time-standardisation concern about **temporal richness curves** (Fig. 4 style). After rarefaction, regional producer richness is only mildly lower post-1600 and consumers are flat — so do **not** hang the title on a long-term decline in rarefied regional richness.

That is **not** the Mendoza & Araújo “simplification” claim used in this paper. Simplification here means increasing prevalence of **species-poor CTS** (yellow CTS1). That claim is tested by Fig. 3d (CTS occupancy over time) + Fig. 3e / **A10** (CTS1 remains the lowest-richness structure with and without rarefaction). A1 does not overturn it.

Figures: `outputs/revision/figures/A1_fig4_rarefied_richness.png`, `A1_fig3e_lake_rarefied_richness.png`, `A1_sample_counts.png`

---

## A2 — Resolution-matched turnover (R2c, R3#8)

| Resolution | Producers earlier | Consumers earlier | Mean mag. diff (P−C) |
|---|---:|---:|---:|
| species | 6 | 2 | 0.51 |
| genus (diatoms) | 6 | 1 | 0.01 |
| guild | 5 | 1 | 0.07 |

**Implication:** “Earlier” largely survives (≈6 of 9 at species/genus). “Stronger” collapses once diatoms are genus-aggregated or compared in guild space — magnitude differences near zero. Prefer Option B wording from the revision plan (“producer turnover preceded consumer turnover in 6 of 9 lakes”; drop “stronger”).

---

## A3 — Coverage / common-period (R3#10)

Coverage diagnostics and varpart under:

- full record
- 1400–2000 CE only
- bins with `n_lakes ≥ 6`
- island-balanced regional means
- lake-conditioned partial RDA

Full-record phase pattern (adj. R²): vegetation dominates Phase 3 (~0.56 pure), Phase 4 (~0.25), with NAO rising in Phase 5 (~0.11). Common-period run only retains Phases 4–5 (as expected). Early-phase (1–2) R² is effectively zero / unsupported — qualify “pristine phase” as Pico/Corvo-dominated.

---

## A4 — NAO–vegetation correlation (R3#2)

| Slice | Pearson r | p | Spearman ρ | p |
|---|---:|---:|---:|---:|
| overall | −0.37 | 0.002 | −0.41 | <0.001 |
| pre-1300 | 0.01 | 0.95 | −0.03 | 0.84 |
| post-1300 | 0.23 | 0.30 | 0.26 | 0.24 |

Overall weak negative association; **no** clear coupling within pre- or post-1300 slices. Supports treating NAO and vegetation as partially independent predictors.

---

## A5 — Human vs climatic vegetation (R3#3)

Median first presence: **Plantago ~1445 CE**, **Cereals ~1522 CE** (site-level first positives; check record-start artefacts).

Arboreal % vs NAO is weak within eras (pre/post 1300 ns) but positive overall (r ≈ 0.38) — largely a shared long-term trend, not evidence that NAO drives clearance. Soften “catchment vegetation” to **regional vegetation**.

---

## A6 — Volcanism sensitivity (R3#5–6)

Tephra inventory: **27 lake-specific** rows in `data/revision/tephra_events.csv` (regional-only rows removed; raw + refs under `data/revision/tephra_*`). A6 uses **18 SUPPORTED/TENTATIVE** ages (7 lakes: Azul, Caveiro, Empadadas Norte, Ginjal, Peixinho, Prata, Santiago) as a **predictor**, not a filter.

**Predictor coding:** on the regional joined 30-yr series, `tephra` = 1 if the bin centre is within ±1 thirty-year bin of any retained lake-specific tephra age, else 0 (27/66 bins = 1). Optional intensity = count of lakes with overlapping influence windows. Full n retained in all phases for the main test.

**Primary result — varpart with vs without tephra** (adj. R²; Pure estimate = vegetation):

| Phase | Pure NAO (base → +tep) | Pure veg (base → +tep) | Pure tephra | Shared (+tep) | Total AB → ABC | n |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 0 → 0 | 0 → 0 | 0 | 0 | 0 → 0 | 23 |
| 2 | 0 → 0 | 0 → 0 | 0 | 0 | 0 → 0 | 10 |
| 3 | 0.088 → 0.104 | 0.559 → 0.246 | **0.054** | 0.184 | 0.533 → 0.587 | 14 |
| 4 | 0 → 0.073 | 0.254 → 0.096 | **0.031** | 0.086 | 0.256 → 0.288 | 10 |
| 5 | 0.114 (base) | 0.140 (base) | — | — | — | 8 |

Phase 5 has no tephra-influenced bins in this window, so the 3-predictor model is undefined there (baseline unchanged).

**Implication:** Tephra explains little unique variance (≤ ~0.05–0.07 adj. R²). Adding it reallocates some Phase 3–4 vegetation pure-R² into shared fractions (tephra windows overlap clearance-era bins), but vegetation still dominates Phase 3 and NAO still rises in Phase 5. **Do not** treat volcanism as overturning the NAO/vegetation conclusions. An exclusion diagnostic (±1 bin drop) still destroys n (Phase 3: 14→5) and is **not** the main A6 result.

Figures/tables: `outputs/revision/figures/A6_tephra_varpart.png`, `A6_varpart_tephra_sensitivity.csv`, `A6_varpart_delta.csv`

---

## A7 — CTS5/CTS6 split (R3#4)

Before 750 CE, k = 6 recovers a distinct **CTS5** (~15% of early samples) that k = 5 / lumping hides. Early composition is still dominated by CTS2 (~67–69%), so “relative early stability” can be kept if qualified, but lumping CTS5→CTS5/6 should be justified or reversed.

---

## A8 — Morphometry / CTS1 (H4)

CTS1 occupancy vs Zmax is weak (Spearman ≈ 0.10, n = 8). Deep lakes (Funda, Santiago, Azul) do not uniquely monopolise CTS1 under the k = 6 ranking used here. Fill Table S3 TP / trophic state for a proper test; Prata still lacks morphometry.

---

## A9 — Onset vs stocking (H3)

Species-level onsets: producers earlier in 6 lakes; consumers earlier in **Caveiro** and **Santiago**; Caldeirão similar.

Azul stocking = 1792 CE; both producer (1486) and consumer (1740) onsets **precede** stocking, so Azul is not a post-stocking consumer-first case under the ≥800 CE onset definition. Need stocking dates for Caveiro and Santiago before testing the top-down mechanism.

---

## A10 — Fig. 3e observed vs rarefied richness by CTS

Paired boxplots for each CTS (CTS5→CTS1): **observed** total species richness (diatoms + chironomids; uses stored `total_nspp_by_lake_core` when available) and **rarefied** richness (diatoms to n = 400; chironomids to n = 20; sample-level, ≥70% retention).

| CTS | Observed median | Rarefied median |
|---|---:|---:|
| CTS5 | 44 | 34.8 |
| CTS4 | 43 | 31.6 |
| CTS3 | 41 | 32.0 |
| CTS2 | 28.5 | 26.9 |
| CTS1 | 25 | 20.8 |

ANOVA p ≪ 0.001 for both series. The richness gradient across CTS persists after rarefaction (CTS5 > CTS1), but absolute differences shrink.

**Implication for the simplification claim (Mendoza & Araújo):** CTS1 (yellow) remains the **species-poorest** structure under rarefaction (median 25 → 20.8). So if Fig. 3d still shows rising CTS1 prevalence toward the present, the structural “simplification” argument holds with or without rarefaction. A10 supports keeping that framing; it does not require retitling away from simplification on richness-gradient grounds.

Figures: `outputs/revision/figures/A10_fig3e_observed_vs_rarefied_cts.png`, `A10_fig3e_observed_vs_rarefied_cts_dodged.png`

---

## A11 — Alternative Figure 5 with tephra overlay

Alternative variance-partitioning figure (Fig. 5 a–c) with the **A6 tephra signal** overlaid on every panel (annotation only; no bins dropped).

**Panels:**
- **(a)** Historical phases — stacked Pure Climate / Pure Vegetation / Shared (adj. R²)
- **(b)** Moving-window varpart (300-yr window, 30-yr step)
- **(c)** Effect size (Vegetation − Climate) by window

**Tephra visual:** lake-specific SUPPORTED/TENTATIVE ages (`include_sensitivity`; same coding as A6), ±1×30-yr merged windows. Panels **(b–c)**: vertical taupe year bands. Panel **(a)**: full-column taupe shading for phases that overlap any tephra window (categorical x-axis). Legend/caption: “Tephra (±1×30-yr)”.

Figures: `outputs/revision/figures/A11_alt_figure5_with_tephra.png`, `A11_alt_figure5_with_tephra.pdf`  
Script: `scripts/revision/A11_alt_figure5_tephra.R`  
Tables: `A11_tephra_bands.csv`, `A11_fig5_phase_varpart.csv`, `A11_fig5_window_varpart.csv`, `A11_fig5_effect_diff.csv`

---

## Suggested claim updates (for response letter)

1. Drop or heavily qualify “stronger” producer response (A2).
2. Keep **simplification** as rising prevalence of species-poor CTS1 (Fig. 3d + A10/Fig. 3e); do not equate it with a regional rarefied-richness decline (A1). Soften any Abstract/Discussion lines that say “long-term declines in species richness” unless restricted to consumers or to CTS composition.
3. Restrict early-phase / pristine statements to lakes actually covering that interval (A3).
4. Soften catchment → regional vegetation (A5).
5. Keep “earlier in 6 of 9 lakes”; name Caveiro and Santiago as consumer-first (A2/A9).
