# Revision analyses A1–A9 + A1x — results summary

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

That is **not** the Mendoza & Araújo “simplification” claim used in this paper. Simplification here means increasing prevalence of **species-poor CTS** (yellow CTS1). That claim is tested by Fig. 3d (CTS occupancy over time) + Fig. 3e / **A1x** (CTS1 remains the lowest-richness structure with and without rarefaction). A1 does not overturn it.

Figures: `outputs/revision/figures/A1_fig4_rarefied_richness.png`, `A1_fig3e_lake_rarefied_richness.png`, `A1_sample_counts.png`

See also **A1x** (observed vs rarefied richness by CTS), which follows A1 and addresses the same rarefaction / simplification question at the CTS level.

---

## A1x — Fig. 3e observed vs rarefied richness by CTS

Companion to **A1**: paired boxplots for each CTS (CTS5→CTS1) addressing whether the CTS richness gradient (and thus the structural simplification claim) holds after rarefaction. **Observed** total species richness (diatoms + chironomids; uses stored `total_nspp_by_lake_core` when available) vs **rarefied** richness (diatoms to n = 400; chironomids to n = 20; sample-level, ≥70% retention).

| CTS | Observed median | Rarefied median |
|---|---:|---:|
| CTS5 | 44 | 34.8 |
| CTS4 | 43 | 31.6 |
| CTS3 | 41 | 32.0 |
| CTS2 | 28.5 | 26.9 |
| CTS1 | 25 | 20.8 |

ANOVA p ≪ 0.001 for both series. The richness gradient across CTS persists after rarefaction (CTS5 > CTS1), but absolute differences shrink.

**Implication for the simplification claim (Mendoza & Araújo):** CTS1 (yellow) remains the **species-poorest** structure under rarefaction (median 25 → 20.8). So if Fig. 3d still shows rising CTS1 prevalence toward the present, the structural “simplification” argument holds with or without rarefaction. A1x supports keeping that framing; it does not require retitling away from simplification on richness-gradient grounds.

Figure: `outputs/revision/figures/A1x_fig3e_observed_vs_rarefied_cts.png`  
Script: `scripts/revision/A1x_fig3e_rarefied_cts.R`  
Tables: `A1x_cts_observed_vs_rarefied_richness.csv`, `A1x_cts_richness_summary.csv`, `A1x_anova_pvalues.csv`

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

Tephra inventory: **27 lake-specific** rows in `data/revision/tephra_events.csv` (no regional rows). A6 uses **11 SUPPORTED-only** ages (`include_sensitivity`; 6 lakes) as a **lake-specific predictor**, not a regional curve and not a data filter. **Seven rows dropped** from the previous SUPPORTED/TENTATIVE set: four **Caveiro proxy-inferred** layers (100, 500, 1350, 1615 CE; Björck Table 3, no visible tephra) plus three **TENTATIVE** Empadadas/Ginjal entries (Furnas I / Fogo 1563 / Ginjal basal ash).

**Predictor coding (lake × 30-yr bin matrix):** for each lake×bin, `tephra = 1` if that bin centre is within ±1×30-yr of **that lake’s** retained SUPPORTED ages, else 0. This flags **13/283** lake×bins (4.6%; 4 lakes with any tephra flag — Santiago/Peixinho inventory ages no longer overlap bins after proxy removal). No regional “any-lake tephra” indicator.

**Design (prefer apples-to-apples):** both columns are **lake-level with `Condition(lake)`** on `bio_lake_30` joined to regional NAO + vegetation (± lake-specific tephra). This is a **lake-level alternative Fig. 5**, not a replot of the published regional Fig. 5. Significance: residual permutation tests (`n_perm = 999`; published Fig. 5 used 9999) with faded NS fractions (Shared always opaque).

**Layout — `A6_alt_figure5_standard_vs_tephra` (2 columns × 3 rows):**
- **Left:** NAO + vegetation only (+ Condition(lake))
- **Right:** NAO + vegetation + lake-specific tephra where `tephra` varies; **NAO + vegetation-only fallback** (same partial-RDA framework) where tephra is invariant, marked with **\***
- **(a)** Historical phases — stacked Pure Climate / Pure Vegetation / Pure tephra (right, full model only) / Shared
- **(b)** Moving window (300-yr window, 30-yr step) — same components
- **(c)** Effect size (Vegetation − Climate) by window (tephra partialled on the right when estimable)

**Primary phase result** (adj. R² after Condition(lake); Pure estimate = vegetation):

| Phase | Pure NAO (base → +tep) | Pure veg (base → +tep) | Pure tephra (p) | Right column |
|---|---:|---:|---:|---|
| 1 | 0 → 0 | 0 → 0 | 0.005 (p = 0.204, NS) | full +tephra |
| 2 | 0 → 0 | 0 → 0 | 0 (p = 0.966) | full +tephra |
| 3 | ~0 → ~0 | **0.062 → 0.059 (p = 0.001)** | 0 (p = 0.663) | full +tephra |
| 4 | 0 → 0 | **0.013 (p = 0.037)** | — | **\* fallback** (no tephra flags in phase) |
| 5 | 0.003 (NS) | **0.014 (p = 0.032)** | — | **\* fallback** (no tephra variation) |

Phases 4–5 on the right repeat the left-column NAO + vegetation partition (asterisk); they are not empty and do not imply a null tephra effect.

**Early moving-window pure tephra (470–660 CE midpoints; SUPPORTED-only; n_perm = 999):** after removing proxy-inferred Caveiro layers, **no window remains significant** at α = 0.05:

| Midpoint (CE) | Window | Pure tephra adj. R² | p (pure tephra) | Sig? |
|---:|---|---:|---:|:---:|
| 480 | 330–630 | 0.013 | 0.224 | no |
| 510 | 360–660 | 0.023 | 0.120 | no |
| 540 | 390–690 | 0.035 | 0.068 | no |
| 600 | 450–750 | 0.006 | 0.302 | no |
| 660 | 510–810 | 0.013 | 0.222 | no |

The earlier significant early-tephra pattern was driven largely by **proxy-only Caveiro ages** in the old inventory; with SUPPORTED layers only, early pure-tephra fractions are small and NS. The later vegetation-dominated interval (~1000–1400 CE) is unchanged.

**Post-~1780 right column:** windows with no tephra variation (e.g. midpoints 1800, 1830, 1860) now show the **NAO + vegetation fallback** with \*, matching the left column for those windows rather than blank panels. Tephra **presence** ticks on right panels (b)–(c) mark unique 30-yr bin centres where lake-specific `tephra = 1`.

**Caption / encoding:** figure caption uses UTF-8 “±1 × 30-yr” and explains full +tephra vs \* fallback bars.

**Interpretation (response letter / SI caption):**  
Lake-level alternative Fig. 5 contrasts NAO + vegetation (left) with the same design plus lake-specific **SUPPORTED** tephra where it varies (right), both with Condition(lake). Pure tephra is small in Phases 1–3 and never significant at α = 0.05 after the stricter inventory; pure vegetation in Phase 3 remains the main significant unique fraction and is barely altered by tephra. Phases 4–5 and late moving windows on the right use an explicit NAO + vegetation fallback (\*) when tephra is invariant, so no information is hidden as empty cells. Early-window pure-tephra peaks that appeared with proxy-inferred Caveiro layers **do not survive** SUPPORTED-only filtering. Overall, defensible lake-specific volcanism does not reallocate the main NAO/vegetation chronology.

Figures/tables: `outputs/revision/figures/A6_alt_figure5_standard_vs_tephra.png`, `A6_varpart_tephra_sensitivity.csv`, `A6_varpart_delta.csv`, `A6_varpart_results_with_p.csv`, `A6_window_varpart.csv`, `A6_window_sig.csv`, `A6_alt_fig5_*.csv`, `A6_tephra_lake_bin_matrix.csv`, `A6_tephra_presence_ages.csv`, `A6_phase_sig_*.csv`  
Script: `scripts/revision/A06_volcanism_sensitivity.R`; tephra build: `scripts/revision/_build_tephra_events.R`

---

## A7 — CTS5/CTS6 split (R3#4)

**What was wrong (initial revision):** A7 **re-ran AMD separately at k = 5 and k = 6** (300 iterations) and compared those partitions in the temporal/guild figures. That is **not** the published workflow in `main_script.Rmd`, which runs **one k = 6** solution (asymptote-selected), **relabels raw clusters by mean `total_nspp_by_lake_core`**, then **merges richness ranks 6 into 5** (`amd_clusts[amd_clusts == 6] <- 5`). The old script only wrote a naive `CTS6 → CTS5` lump on euplanctonic labels to CSV (`k6_lumped_to_5`) and **did not** use it in the figures—so guild profiles for “k = 5” in the plot were from an **independent k = 5 AMD**, not the manuscript merge.

**Corrected comparison (three schemes):**

| Scheme | Definition |
|--------|------------|
| **k = 5 (independent AMD)** | Separate fuzzy AMD at k = 5, euplanctonic CTS labels |
| **k = 6** | Single k = 6 AMD partition (5000 iterations, same seed as merge path) |
| **k = 5 (merged from k = 6, manuscript)** | Same k = 6 partition → richness rank → merge ranks 6+5 → euplanctonic CTS1–5 |

**CTS display labels:** Euplanctonic reordering (highest → **CTS1**, viridis **yellow**) for all A7 panels. Manuscript **merge criterion** uses **`total_nspp_by_lake_core`** on raw k = 6 clusters (`A7_k6_manuscript_merge_map.csv`, `A7_cts_diversity_by_scheme.csv`).

**Manuscript vs revision CTS names:** In `main_script.Rmd`, CTS1–CTS5 are labelled by **ascending richness** after the k = 6 → k = 5 merge (CTS1 = lowest richness / euplanctonic-dominated raw cluster 5; the merged pair = **CTS5** = highest richness). Revision **A7** and **A1x** use **euplanctonic reorder** (`revision_assign_cts` / `revision_manuscript_k6_merged_to_k5`) for **colours and comparability with A8 k = 6** (yellow CTS1). That can swap **CTS3 vs CTS4** names relative to manuscript richness order; **sample membership** of merged k = 5 still matches the manuscript merge.

**Downstream:** `A7_sample_cts5_assignments.csv` is now the **manuscript merged k = 5** (for A1x / Fig. 3e). **A8** still uses **k = 6** assignments for euplanctonic **CTS1** vs morphometry.

Temporal occupancy is normalised **within each 30-yr Age(CE) bin** (proportions sum to 1; `A7_cts_temporal_bin_sums.csv`).

**Interpretation:** Independent k = 5 AMD can split guild structure differently from k = 6 (e.g. CTS3 mean profiles need not match). The **merged** k = 5 reuses the k = 6 sample partition and only collapses **richness ranks 5 and 6**—the two **highest** mean-`total_nspp_by_lake_core` k = 6 groups (raw AMD clusters **2** and **6** after richness relabelling; manuscript `amd_clusts[amd_clusts == 6] <- 5`). Guild means for merged euplanctonic **CTS5** are **mixtures** of k = 6 **CTS4** + **CTS6**, not equal to naive euplanctonic `CTS5`+`CTS6` relabelling without the richness merge.

Figures: `outputs/revision/figures/A7_cts_guild_profiles_k5_k6_merged.png` (alias of updated guild plot), `A7_cts_temporal.png` (three facets).  
CSVs: `A7_cts_guild_profiles.csv`, `A7_k6_manuscript_merge_map.csv`, `A7_cts_diversity_by_scheme.csv`, assignment files `A7_sample_cts5_*` / `A7_sample_cts6_assignments.csv`.  
Script: `scripts/revision/A07_split_cts5_cts6.R`

---

## A8 — Morphometry / CTS1 (H4; reviewer L415 / L484)

With **corrected CTS1** (euplanctonic-dominated, k = 6), CTS1 occupancy vs **Zmax** is **strong and positive** (Spearman **ρ ≈ 0.81, p ≈ 0.016, n = 8**). Deep lakes carry more CTS1 (Santiago ≈ 0.95, Funda ≈ 0.41; shallow lakes ≈ 0). **Prata is excluded** from depth tests: it is now a peatland (Zmax = 0 / unavailable in metadata).

Same Spearman tests (n limited by missing Table S3):
- CTS1 vs **lake area**: ρ ≈ 0.07, p ≈ 0.88 (n = 7; Ginjal area missing)
- CTS1 vs **trophic state / TP**: not estimable (only Caldeirão has Table S3 values; n = 1)
- **Turnover** vs Zmax: ρ ≈ −0.10, p ≈ 0.84 (n = 8); **Funda** is a high-turnover outlier. Area/trophic likewise NS or missing.

**Implication:** after fixing CTS labelling, morphometry (depth) supports H4 for euplanctonic CTS1; fill Table S3 before claiming trophic-state contingency.

---

## A9 — Onset vs stocking (H3)

Species-level onsets: producers earlier in 6 lakes; consumers earlier in **Caveiro** and **Santiago**; Caldeirão similar.

Azul stocking = 1792 CE; both producer (1486) and consumer (1740) onsets **precede** stocking, so Azul is not a post-stocking consumer-first case under the ≥800 CE onset definition. Need stocking dates for Caveiro and Santiago before testing the top-down mechanism.

---

## Suggested claim updates (for response letter)

1. Drop or heavily qualify “stronger” producer response (A2).
2. Keep **simplification** as rising prevalence of species-poor CTS1 (Fig. 3d + A1x/Fig. 3e); do not equate it with a regional rarefied-richness decline (A1). Soften any Abstract/Discussion lines that say “long-term declines in species richness” unless restricted to consumers or to CTS composition.
3. Restrict early-phase / pristine statements to lakes actually covering that interval (A3).
4. Soften catchment → regional vegetation (A5).
5. Keep “earlier in 6 of 9 lakes”; name Caveiro and Santiago as consumer-first (A2/A9).
