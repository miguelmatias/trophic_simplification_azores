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

## A6 — Volcanism sensitivity / alternative Figure 5 (R3#5–6)

Tephra inventory: **27 lake-specific** rows in `data/revision/tephra_events.csv` (no regional rows). A6 uses **18 SUPPORTED/TENTATIVE** ages (`include_sensitivity`; 7 lakes) as a **lake-specific predictor**, not a regional curve and not a data filter.

**Predictor coding (lake × 30-yr bin matrix):** for each lake×bin, `tephra = 1` if that bin centre is within ±1×30-yr of **that lake’s** retained tephra ages, else 0. This flags **28/283** lake×bins (9.9%; 6 lakes with any tephra flag in the abundance table — Santiago has inventory ages but no overlapping bio bins). No regional “any-lake tephra” indicator.

**Design (prefer apples-to-apples):** both columns are **lake-level with `Condition(lake)`** on `bio_lake_30` joined to regional NAO + vegetation (± lake-specific tephra). This is a **lake-level alternative Fig. 5**, not a replot of the published regional Fig. 5. Significance: residual permutation tests (`n_perm = 999`; published Fig. 5 used 9999) with faded NS fractions (Shared always opaque).

**Layout — `A6_alt_figure5_standard_vs_tephra` (2 columns × 3 rows):**
- **Left:** NAO + vegetation only (+ Condition(lake))
- **Right:** NAO + vegetation + lake-specific tephra (+ Condition(lake))
- **(a)** Historical phases — stacked Pure Climate / Pure Vegetation / Pure tephra (right only) / Shared
- **(b)** Moving window (300-yr window, 30-yr step) — same components
- **(c)** Effect size (Vegetation − Climate) by window (tephra partialled on the right)

**Primary phase result** (adj. R² after Condition(lake); Pure estimate = vegetation):

| Phase | Pure NAO (base → +tep) | Pure veg (base → +tep) | Pure tephra (p) | Shared (+tep) | n bins / lakes |
|---|---:|---:|---:|---:|---:|
| 1 | 0 → 0 | 0 → 0 | **0.014 (p = 0.066, NS)** | 0.004 | 43 / 4 |
| 2 | 0 → 0 | 0 → 0 | 0 (p = 0.966) | 0 | 28 / 5 |
| 3 | ~0 → ~0 | **0.062 → 0.062 (p = 0.001)** | 0 (p = 0.831) | 0 | 71 / 8 |
| 4 | 0 → 0 | **0.013 → 0.014 (p = 0.023)** | 0 (p = 0.959) | 0 | 73 / 9 |
| 5 | 0.003 (NS) | **0.014 (p = 0.032)** | — (not estimable) | — | 63 / 9 |

**Phase 5 right blank is not “zero tephra effect”:** Phase 5 (≥1750 CE; n = 63 bins / 9 lakes) has `tephra_prevalence = 0` (last lake×bin tephra flag is 1620 CE). The 3-predictor model is **skipped** when `length(unique(tephra)) < 2` (`tephra_varies = FALSE` → all adj. R² / p = NA). Left Phase 5 is the standard 2-predictor fit only; do not read the empty right bar as a null result. Pure-vegetation significance with tephra is therefore claimable for **Phases 3–4 only** on the right (Phase 3: 0.062, p = 0.001; Phase 4: 0.014, p = 0.023); Phase 5 vegetation (0.014, p = 0.032) is left-column only.

**Early moving-window pure tephra (selected midpoints; n_perm = 999):** several ~470–660 CE windows are significant at α = 0.05 and plot at full opacity:

| Midpoint (CE) | Window | Pure tephra adj. R² | p (pure tephra) | Sig? |
|---:|---|---:|---:|:---:|
| 480 | 330–630 | 0.049 | 0.031 | yes |
| 510 | 360–660 | 0.045 | 0.039 | yes |
| 540 | 390–690 | 0.043 | 0.050 | borderline |
| 600 | 450–750 | 0.047 | 0.033 | yes |
| 660 | 510–810 | 0.037 | 0.045 | yes |

(Pure veg/NAO in these windows are NS.) This is a plausible **local volcanic signal before the main vegetation-change interval** in the Pico/Corvo/Prata–dominated early record; it does not overturn the later vegetation-dominated pattern (~1000–1400 CE).

**Post-~1780 right-side gap:** windows with midpoints 1800, 1830, 1860 (n = 78–80 bins) have `tephra_prevalence = 0` / `tephra_varies = FALSE` and are **dropped** (NA, not a failed RDA fit). Midpoint 1770 still has one tephra-flagged bin (prev ≈ 0.013) and is estimable. Left-column late windows remain.

**Interpretation (response letter / SI caption):**  
This lake-level alternative Figure 5 places the standard NAO + vegetation partitioning (left) beside the same design with lake-specific tephra added (right), both after Condition(lake), so the tephra contrast is apples-to-apples rather than a regional-vs-lake comparison. Across historical phases where tephra varies (Phases 1–4), pure tephra is small: Phase 1 adj. R² ≈ 0.014 (p ≈ 0.066, NS); Phases 2–4 ≈ 0 (all NS). Pure vegetation remains the only significant unique fraction in Phases 3–4 and is essentially unchanged when tephra is included; Phase 5 (+tephra) is not estimable (no tephra variation after 1620 CE), so the empty right Phase 5 bar is “n/a”, not a null effect. Moving-window panels keep the vegetation-dominated interval (~1000–1400 CE) on both sides, but early windows (~480–660 CE) show significant pure-tephra fractions (adj. R² ≈ 0.037–0.049; p ≈ 0.031–0.045) consistent with a local volcanic signal before major vegetation change. Right-side windows after ~1780 are missing because tephra is invariant (dropped), not because fits failed. Magnitudes are smaller than published regional Figure 5 because Condition(lake) absorbs among-lake structure; compare left vs right within this figure. Overall, lake-specific volcanism does not reallocate the main NAO/vegetation chronology, but early tephra is not uniformly non-significant.

Figures/tables: `outputs/revision/figures/A6_alt_figure5_standard_vs_tephra.png` (+ PDF), `A6_tephra_varpart.png` (phase-only extract of row a), `A6_varpart_tephra_sensitivity.csv`, `A6_varpart_delta.csv`, `A6_varpart_results_with_p.csv`, `A6_window_varpart.csv`, `A6_window_sig.csv`, `A6_alt_fig5_*.csv`, `A6_tephra_lake_bin_matrix.csv`, `A6_phase_sig_*.csv`  
Script: `scripts/revision/A06_volcanism_sensitivity.R`

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

## A11 — Fig. 5 rebuild (no tephra bars)

Rebuild of variance-partitioning panels (Fig. 5 a–c) with significance fading. **No taupe tephra bands** — volcanism sensitivity is handled in **A6** (lake-level alternative Fig. 5: standard vs +tephra with Condition(lake)), not as shading on Fig. 5.

**Panels:**
- **(a)** Historical phases — stacked Pure Climate / Pure Vegetation / Shared (adj. R²; faded if NS)
- **(b)** Moving-window varpart (300-yr window, 30-yr step)
- **(c)** Effect size (Vegetation − Climate) by window

Figures: `outputs/revision/figures/A11_alt_figure5_with_tephra.png`, `A11_alt_figure5_with_tephra.pdf`  
Script: `scripts/revision/A11_alt_figure5_tephra.R`  
Tables: `A11_fig5_phase_varpart.csv`, `A11_fig5_window_varpart.csv`, `A11_fig5_effect_diff.csv`

---

## Suggested claim updates (for response letter)

1. Drop or heavily qualify “stronger” producer response (A2).
2. Keep **simplification** as rising prevalence of species-poor CTS1 (Fig. 3d + A10/Fig. 3e); do not equate it with a regional rarefied-richness decline (A1). Soften any Abstract/Discussion lines that say “long-term declines in species richness” unless restricted to consumers or to CTS composition.
3. Restrict early-phase / pristine statements to lakes actually covering that interval (A3).
4. Soften catchment → regional vegetation (A5).
5. Keep “earlier in 6 of 9 lakes”; name Caveiro and Santiago as consumer-first (A2/A9).
