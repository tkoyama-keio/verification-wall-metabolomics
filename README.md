# The Verification Wall in Microbial Metabolomics

Analysis code and results for the paper, presented at the AI for Science
Workshop, NeurIPS 2026.

Koyama, T. The Verification Wall in Microbial Metabolomics.
AI for Science Workshop (NeurIPS 2026).

## Contents

- `tier_pipeline.py` sorts candidate metabolite structures into verification
  tiers. Every threshold is a named entry in one configuration dictionary at
  the top of the file.
- `pipeline_validation.md` records tool versions, the provenance of each
  threshold, measured throughput, and known limitations.
- `data/hmdb_verification_tiers.csv` has one row per compound,
  10,645 rows and 38 columns. Each carries its tier, the rule
  that fired, and the evidence behind it, so assignments are auditable from
  this table alone.
- `data/hmdb_sensitivity_sweep.csv` has the 72 threshold combinations
  behind the range reported in the paper.

## Input data

The input is the complete `hmdb_metabolites.xml` export from HMDB 5.0, which
is not redistributed here. Download it from hmdb.ca. The appendix of the
paper records the file timestamp used, 2021-10-23, and HMDB 5.0 remains the
current release.

## Reproducing the tier assignment

Retrosynthesis uses AiZynthFinder with the USPTO models and a ZINC stock
file. Versions are listed in `pipeline_validation.md`. Purchasability is
resolved against PubChem categorized sources. Both steps need network access
and take a few hours over the full set.

## Caveats that matter when reading the table

Crystallizability is a property-based screen, not a measurement, and was
never tested against crystallization outcomes. A `route_status` of
`no_route_found_TOOL_LIMITED` means the planner found no route, not that none
exists, and 97.0% of the residual tier is of this kind. Vendor annotation does
not establish that a certified analytical standard exists; the separate
`has_reference_standard_vendor` column is the stricter test.
