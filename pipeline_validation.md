# Compound Verification-Tier Pipeline — Validation Report

**Status:** pipeline built and validated end-to-end against a **DEVELOPMENT stand-in set**.
The real HMDB predicted-metabolite export has **not** been processed. No number in this
document is a result about HMDB, and the development-set tier fractions must never be
reported as findings.

**Date:** 2026-08-04 · **Dev set:** 227 unique real metabolites from PubChem

---

## 1. Headline: retrosynthesis feasibility

**The full set is feasible. Sampling is not required.**

| Quantity | Measured value |
|---|---|
| Engine | AiZynthFinder **4.4.1** (local, MIT, no credentials, no rate limit) |
| Benchmark | 227 dev compounds, 64 worker processes |
| Wall time | **89.1 s** |
| **Throughput @64 workers** | **9,172 compounds/hour** |
| Per-compound latency | mean 12.64 s, median 12.45 s, p90 22.74 s, max 42.52 s |
| Engine errors | **0 / 227** |
| Compounds hitting the 120 s cap | **0 / 227** |
| Projected @200 workers (85% scaling eff.) | ~24,400 compounds/hour |

**All three measurements taken (nothing omitted):**

| Run | n | Workers | Wall | Throughput |
|---|---|---|---|---|
| Dev sample | 40 | 20 | 34.3 s | 4,198 /hr |
| Scale replicate | 160 | 64 | 57.2 s | **10,070 /hr** |
| Full dev set | 227 | 64 | 89.1 s | **9,172 /hr** |

**Worker scaling is real:** 20 → 64 workers raised throughput ~2.2× (4,198 → 9,172–10,070/hr),
i.e. ~68% parallel efficiency for a 3.2× worker increase.

**Set-size scaling is NOT favourable — throughput fell ~9% from n=160 to n=227 at fixed
worker count.** Two caveats on that comparison, both pointing the same way:

- The n=160 run was the 40-compound sample **repeated 4×**, so it contained only 40 unique
  structures. Repeated targets do not benefit from any cache in AiZynthFinder, but they do
  make the per-compound difficulty mix narrower and easier than the full 227.
- The full 227-compound set contains the hard bile-acid and conjugate tail (max latency
  42.5 s vs 18.4 s in the sample), which drags the mean up.

So the honest reading is that the n=160 figure is **optimistic** and the **9,172/hr
full-set number is the one to plan against** — it is the only measurement over the real
difficulty distribution. All projections below use 9,172/hr. Whether throughput rises or
falls with set size beyond n=227 has **not** been measured; the model-load amortization
argument predicts improvement, but the data here do not demonstrate it, and the observed
direction is the opposite.

**Projected full-set runtimes**

| Set size | @64 workers | @200 workers |
|---|---|---|
| 5,000 | 0.5 h | 0.2 h |
| 20,000 | 2.2 h | 0.8 h |
| 50,000 | 5.5 h | 2.1 h |
| 217,000 (whole HMDB, worst case) | 23.7 h | 8.9 h |

### Recommendation

**Run the complete set.** Even the extreme case — all ~217k HMDB entries — finishes
overnight on this machine, and any realistic predicted-metabolite subset (typically
10³–10⁴ compounds) finishes in under three hours. Sampling would buy nothing and would
add a sampling-error term to the paper's headline number for no reason.

*If* a sample is ever needed (e.g. a much larger set, or a slower engine), the Wilson
sample sizes are: ±5% → n=385; ±3% → n=1,068; ±2% → n=2,401; ±1% → n=9,604. At the dev
set's tier-4 rate the interval on n=227 is 4.8% [2.7, 8.5] — too wide for a headline
claim, which is the quantitative reason to prefer the census.

**Practical note:** the dominant cost is *not* search — it is the ~30 s model + 663 MB
ZINC stock load per worker. The runner amortizes this by loading once per worker process
and mapping the whole chunk. **Note this predicts improving throughput with set size, but
that is not what we measured** (§1: throughput fell ~9% from n=160 to n=227). At these set
sizes the harder-compound tail outweighs load amortization. Plan against the measured
9,172/hr, not against an amortization argument.

---

## 2. IBM RXN — not used, as instructed

Per your correction, the configured `IBM_TOKEN` / `IBM_INSTANCE` credentials are **IBM
Quantum** credentials, not IBM RXN for Chemistry. They were **not read, not printed, not
written anywhere, and no authentication attempt was made against `rxn.res.ibm.com`.**
(An early read attempt in this session was denied by the credential broker before your
correction arrived; nothing was retrieved.) The `rxn4chemistry` package is in any case
uninstallable in this environment — every published version pins `requests==2.23.0` or
`beautifulsoup4==4.9.0`, which conflict with the platform's global pip constraints.
AiZynthFinder is the primary and only retrosynthesis engine used.

---

## 3. What was built

`tier_pipeline.py` — importable, documented, with **every threshold in a single `CONFIG`
dict at the top**. No threshold is hard-coded in any function body.

**Pipeline stages**

1. **Standardization** — parse → `Cleanup` → largest-fragment (salt/solvate stripping) →
   neutralization → canonical tautomer → InChIKey. Failures are counted and typed
   (`parse_fail` / `standardize_fail` / `empty_after_std`), never silently dropped.
   *Dev-set result: 232 inputs, **0 parse failures**, 5 synonym duplicates collapsed
   (e.g. phloretic acid = 3-(4-hydroxyphenyl)propionic acid), 227 unique.*
2. **Purchasability** — PubChem CID resolution (InChIKey first, canonical SMILES
   fallback) → PUG-View categorized sources → "Chemical Vendors".
3. **Retrosynthesis** — AiZynthFinder, parallel, with in-stock targets detected separately.
4. **Crystallizability proxy** — RDKit descriptors + explicit threshold rule.
5. **Tier assignment** — four tiers, explicit precedence, per-compound audit trail.

**Environment note:** AiZynthFinder requires `networkx<3.0`, which conflicts with the
platform's global pip constraint (`networkx==3.4.2`). It is installed in an isolated venv
built with `PIP_CONSTRAINT` cleared, and requires
`LD_LIBRARY_PATH=/data/.claude-science/conda/envs/python/lib` at runtime so RDKit's
drawing module can find `libXrender.so.1`. Both are recorded in the runner script.

---

## 4. What "purchasable" operationally means — and does not

**What is actually measured:** whether PubChem's PUG-View categorized-sources endpoint
lists ≥ `min_distinct_vendors` distinct suppliers under the **"Chemical Vendors"**
category for the CID that the standardized structure resolves to.

**What that is NOT:**

- **Not evidence that an authentic analytical standard exists.** Most "Chemical Vendors"
  entries are catalog listings for research-grade material with no certificate of
  analysis, no stated purity, and no quantitative reference value. LC-MS/MS verification
  requires a certified reference standard — a far smaller supplier set.
- **Not evidence of current stock.** Many listings are make-to-order, or advertise
  compounds the vendor has never synthesized. PubChem carries no stock or price data.
- **Not stereochemically exact.** A CID match can return a listing for a different
  stereoisomer or a racemate — unacceptable as a standard for a stereochemically defined
  metabolite. This is a live risk for bile acids, where epimers at C-3/C-7/C-12 are
  distinct metabolites with distinct biology.
- **A lagging snapshot**, in both directions.

**Mitigation.** The pipeline reports **two** columns:

| Column | Definition | Dev-set count |
|---|---|---|
| `purchase_status` (permissive) | any vendor annotation | 207 / 227 (91.2%) |
| `has_reference_standard_vendor` (strict) | ≥1 supplier from a curated list of analytical-standard houses (Sigma-Aldrich, Cerilliant, TRC, Steraloids, Cayman, LGC, USP, …) | 179 / 227 (78.9%) |

The permissive definition is the **default** because it inflates tier 1 and therefore
*shrinks* the headline tier-4 residual — the conservative direction for a paper arguing a
verification gap exists. `require_reference_standard_vendor_for_tier1: True` switches to
the strict definition; on the dev set this moves 28 compounds out of tier 1.

**Honest caveat for the paper:** even the strict column is a vendor-name heuristic, not a
catalog check. A defensible tier-1 claim ultimately needs a real catalog query (Sigma,
TRC, Cayman APIs) or manual spot-checking of a random subsample. Recommend manually
verifying ~50 sampled tier-1 assignments before publication and reporting that agreement rate.

**Caching & rate limits.** Every response (including 404s) is cached to
`cache/pubchem/` keyed by URL hash; re-runs cost zero requests. Requests are throttled to
≤4/s (under PubChem's documented 5/s cap) with exponential backoff on 503/429/500.
Observed dev-set cost: ~3 s/compound, dominated by PUG-View payload size (~190 KB/compound).
13/227 structures did not resolve to any CID.

---

## 5. Retrosynthesis: engine provenance and the TOOL-LIMITED distinction

**Exact configuration** (record this in the paper's methods):

| Component | File |
|---|---|
| Engine | aizynthfinder 4.4.1 (pip) |
| Expansion policy | `uspto_model.onnx` + `uspto_templates.csv.gz` |
| Ringbreaker policy | `uspto_ringbreaker_model.onnx` + templates |
| Filter policy | `uspto_filter_model.onnx` |
| **Stock** | `zinc_stock.hdf5` — **17,422,831 InChIKeys** |
| Search budget | 120 s / 100 iterations per compound |

All obtained via AiZynthFinder's own `download_public_data` mechanism.

**"In stock" means exactly:** membership in that 17.4M-InChIKey ZINC snapshot. It is not a
live vendor check and carries no price, purity, or lead-time. A building block "in stock"
by this definition may be unobtainable in practice; a genuinely purchasable reagent absent
from the snapshot is treated as not-in-stock.

### "No route found" is TOOL-LIMITED — with direct evidence

`RouteStatus.NO_ROUTE_FOUND_TOOL_LIMITED` is a distinct enum value, never conflated with
"hard to synthesize". The dev set produced decisive evidence for why this matters:

- **AiZynthFinder found no route to trimethylamine** — a bulk commodity chemical made
  industrially at megatonne scale from methanol and ammonia.
- **It failed on essentially every bile acid** — chenodeoxycholic, hyodeoxycholic,
  allocholic, isodeoxycholic, β-muricholic, murideoxycholic, and the taurine/glycine
  conjugates.

These are not difficult molecules. Their real routes — C1 amination chemistry, steroid
semi-synthesis from natural feedstocks — are absent or under-represented in the USPTO
template library, which encodes medicinal-chemistry reaction space and systematically
under-covers the steroid, sugar, and natural-product chemistry that dominates gut
microbiome metabolites. **This is precisely the chemical space of the paper.**

Dev-set outcome distribution: `already_in_stock` 116, `no_route_found_TOOL_LIMITED` 75,
`solved_short_route` 31, `solved_long_route` 5.

**Implication for the paper:** the tier-2 → tier-3/4 boundary is the *softest* boundary in
the analysis, and it is biased in a known direction — toward over-populating tiers 3 and 4.
The residual tier is therefore an **upper bound** on the true verification gap, and should
be presented as such. Reporting the 75 TOOL_LIMITED compounds as a separate, explicitly
labelled column is strongly recommended. Consider a template library trained on a broader
corpus (e.g. Reaxys/Pistachio-derived, if licensing permits) as a robustness check.

---

## 6. Crystallizability proxy — every threshold, with provenance

**This is a negative screen (a proxy), not a prediction of crystallizability.**

| Parameter | Value | Grounding |
|---|---|---|
| `cryst_mw_max` | 350 | **Literature-derived — citation not yet verified.** From the crystallization-propensity literature (Wicker & Cooper line of work), where low MW + low flexibility is reported as a crude predictor of ready crystallization. **See caveat below.** |
| `cryst_rotb_max` | 4 | **Literature-derived — citation not yet verified**, same source. |
| `cryst_frac_rotb_max` | 0.30 | **Judgment call.** Size-normalizes flexibility so large rigid polycyclics (steroid nuclei) aren't penalized for size alone. No literature anchor. |
| `cryst_tpsa_max` | 140 | **Judgment call.** Borrowed from Veber's oral-bioavailability rule, which concerns membrane permeation, *not* crystallization. Used as a loose proxy for heavy solvation / hydrate formation. Weakly justified — sweep it. |
| `cryst_hbd_max` | 5 | **Judgment call**, Lipinski-derived. Note this term arguably has the *wrong sign*: strong directional H-bonding often *promotes* good crystal packing. |
| `cryst_logp_min` | −1.0 | **Judgment call.** Very hydrophilic species are hard to crystallize from aqueous media at trace abundance. |
| `cryst_abs_formal_charge_max` | 1 | **Judgment call.** Permanently charged species crystallize as salts whose counterion isn't determined by the predicted structure. |
| `cryst_n_violations_for_hard_flag` | 2 | **Judgment call — the most consequential knob in the pipeline.** |

### ⚠ Citation-verification caveat — action required before submission

The literature groundings above were assembled from search-result summaries, **not from
the primary articles**, and the persisted record of those searches does not contain the
underlying page text. Specifically:

- The **MW < 350 / < 4 rotatable bonds** pairing is the commonly quoted form of a result
  from the crystallization-propensity literature, but **the exact author, year, and
  numeric thresholds have not been checked against the source article.** An earlier draft
  attributed these to "Kleppe (2015)"; that attribution is **not verified** and has been
  removed from the code comments.
- Numeric ML-classifier accuracy figures quoted in an earlier draft (a held-out error rate
  and a blind-screen accuracy) have been **removed entirely** as unverifiable.
- A specific nanometre MicroED crystal-size range has likewise been **removed**.

**Before this pipeline's thresholds are cited in the manuscript, obtain and read:**

1. Wicker & Cooper, *Will it crystallise? Predicting crystallinity of molecular materials*,
   CrystEngComm (2015).
2. *Predicting the Crystallization Propensity of Drug-Like Molecules* (ScienceDirect, ~2016).
3. Wicker & Cooper, *Beyond Rotatable Bond Counts: Capturing 3D Conformational Flexibility
   in a Single Descriptor*, J. Chem. Inf. Model. (2017).
4. Cooper et al., *Can molecular flexibility control crystallization?*, Chem. Sci. (2021)
   — PMC8179050.

Confirm each number against the article, or reclassify the threshold as a judgment call.
This matters more than usual here: a position paper arguing that others overstate their
evidence cannot itself carry an unchecked citation.

### The literature caveat that must be stated in the paper

Cooper and co-workers (*Chem. Sci.* 2021, para-substituted benzoic acids) showed that
conformational flexibility alone does **not** limit crystallizability — a compound with
*more* rotatable bonds nucleated *faster* than less flexible analogues. Their conclusion:
it is the **nature** of the rotatable bonds (the conformational energy barrier relative to
the ~70 kJ/mol nucleation activation energy), not the **count**. A rotatable-bond count is
therefore a coarse surrogate for a quantity it does not actually measure.

Better descriptors exist — `nConf20` (Wicker & Cooper, *JCIM* 2017) and the Kier
flexibility index Φ — and published ML classifiers substantially outperform threshold
rules of this kind. (Specific accuracy figures are deliberately omitted here pending
verification against the primary articles; see the caveat above.) **We compute Kier Φ and
store it, but do not threshold on it.** Thresholds were chosen for transparency and
sweepability, not accuracy. The paper must not claim this is a crystallizability prediction.

**MicroED-specific note:** MicroED works from sub-micron crystals — orders of magnitude
smaller than single-crystal X-ray requires (a specific nanometre range is deliberately not
quoted here pending citation check) — which is why tier 3 is plausible at trace abundance at all — the technique routinely solves
structures from material that looks amorphous by eye, and compounds often crystallize
adequately during ordinary purification. Throughput and instrument access, not crystal
size, are the practical limits. "MicroED-plausible" is a claim about physical possibility,
not routine feasibility.

### Sensitivity of the headline number

Measured on the dev set (purchasability disabled so the proxy is actually exercised):

| `cryst_n_violations_for_hard_flag` | tier 3 | tier 4 | tier 4 % |
|---|---|---|---|
| 1 | 3 | 17 | 7.5% |
| **2 (default)** | **9** | **11** | **4.8%** |
| 3 | 11 | 9 | 4.0% |
| 4 | 14 | 6 | 2.6% |

The headline residual fraction moves by ~3× across a defensible range of a single
judgment-call parameter. **Report the sweep, not a point estimate.**

Additional ionization caveat: `estimate_charge_at_ph` is a SMARTS rule of thumb with no
pKa model behind it. It identifies obvious zwitterions and permanent charges and will be
wrong for anything with an unusual electronic environment.

---

## 7. Tier assignment

Precedence — first matching rule wins:

1. **`1_purchasable`** — vendor annotation satisfies the configured definition.
2. **`2_short_route_synthesizable`** — solved route ≤ `retro_max_steps_for_short_route`
   (default 3), **or** the target is itself an in-stock ZINC building block.
3. **`3_microed_plausible`** — not flagged likely-hard-to-crystallize.
4. **`4_residual`** — everything else.

Every compound carries `tier_rule_fired` and `tier_evidence`, the latter recording each
rule's verdict *and the evidence it fired on* (vendor counts, route status and step count,
the exact list of crystallizability violations). Tier assignment is fully auditable and
reproducible from the CSV alone.

**`is_purchasable()` re-derives its verdict from raw evidence columns**
(`n_distinct_vendors`, `has_reference_standard_vendor`) rather than reading the frozen
`purchase_status` string. This is deliberate: it makes both vendor thresholds genuinely
sweepable post-hoc **without re-hitting the PubChem API**. Reading the cached status
string would have silently pinned those thresholds to their lookup-time values — a bug
caught during validation.

### Dev-set end-to-end run (n=227) — NOT A RESULT

| Configuration | T1 | T2 | T3 | T4 |
|---|---|---|---|---|
| Default (permissive vendors) | 207 | 0 | 9 | 11 |
| Strict (reference-standard vendors only) | 179 | 14 | 12 | 22 |
| `min_distinct_vendors = 5` | 206 | 1 | 9 | 11 |
| Purchasability disabled (branch test) | 0 | 147 | 35 | 45 |

Tier 2 is empty under the default config because 91% of these well-known metabolites are
commercially listed and rule 1 takes precedence — the expected behaviour for a dev set of
*known* compounds, and exactly the opposite of what a set of *predicted* structures will
look like. The last two rows confirm the tier-2 branch fires correctly on real evidence
when tier 1 does not absorb everything: **all four tiers are exercised and verified.**

---

## 8. Known limitations

1. **The dev set is not representative.** 227 well-characterized metabolites chosen for
   coverage of chemical classes. Real HMDB predicted structures will be less
   commercially available, more likely to fail CID resolution, and more likely to hit the
   TOOL_LIMITED path. Dev-set tier fractions carry **no** information about the real set.
2. **Vendor annotation ≠ analytical standard** (§4). The single largest threat to the
   tier-1 claim.
3. **USPTO template bias** (§5) systematically inflates tiers 3 and 4 for exactly the
   natural-product chemistry this paper studies. The residual is an **upper bound**.
4. **The crystallizability proxy is not validated** against any crystallization outcome
   dataset. Thresholds are transparent, not accurate; the headline number moves ~3× across
   a defensible parameter range (§6).
5. **No stereochemistry check in purchasability** — a CID match may return a different
   stereoisomer.
6. **`estimate_charge_at_ph` has no pKa model.**
7. **Tautomer canonicalization can fail** on large systems; the code degrades gracefully
   and records `tautomer_skipped(...)`, but the standardized form may then differ.
8. **ZINC stock is a snapshot**, not a live availability check.
9. **13/227 dev compounds failed CID resolution** — for the real set this rate will be
   higher and must be reported as a distinct category, not folded into tier 4.
10. **Literature citations are unverified** (§6 caveat). The two "literature-grounded"
    thresholds rest on search-summary material, not on the primary articles. Until those
    are read and checked, treat **all seven** crystallizability thresholds as judgment
    calls for the purposes of the paper's claims.

---

## 9. Files

- `tier_pipeline.py` — complete importable pipeline; `CONFIG` at the top.
- `DEV_SET_dev_set_tiers.csv` — 227 rows × 60 columns, dev-set output. Every row carries
  `DEV_SET_WARNING = DEVELOPMENT_SET_NOT_A_RESULT`.
- `DEV_SET_sensitivity_sweep.csv` — threshold sweep.
- `DEV_SET_pipeline_validation.png` — latency, route outcomes, sensitivity.
- `run_retro.py` — parallel AiZynthFinder runner.

## 10. To run on the real HMDB data

1. Set `CONFIG["dataset_label"] = "HMDB_PREDICTED"`.
2. Load the export; map its structure column to `smiles_input`.
3. `standardize_frame()` → `deduplicate()` — **report the parse-failure count.**
4. `assess_purchasability()` per compound (cached; ~3 s/compound, parallelizable).
5. `run_retro.py` with 64–200 workers (~9,000–24,000 compounds/hour).
6. `compute_descriptors()` → `crystallizability_flag()` → `classify_route_result()` →
   `assign_tier()`.
7. Sweep `cryst_n_violations_for_hard_flag` ∈ {1,2,3,4} and
   `require_reference_standard_vendor_for_tier1` ∈ {False, True}; report the range.
