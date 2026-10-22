"""
tier_pipeline.py — Compound verification-tier assignment pipeline.

Assigns each computationally predicted metabolite structure to one of four
verification tiers, for the "verification wall" analysis in gut microbiome
metabolomics:

    TIER 1  purchasable              an authentic standard can simply be bought
    TIER 2  short_route_synthesizable  retrosynthetically reachable in <= N steps
                                     from in-stock building blocks
    TIER 3  microed_plausible        not cheaply synthesizable, but the molecule
                                     does not trip the "hard to crystallize" proxy
    TIER 4  residual                 none of the above -- the irreducible
                                     verification gap

Design commitments (these matter for the paper's argument):

  * EVERY threshold is a named key in CONFIG. Nothing is hard-coded in the
    body of a function, so a sensitivity analysis can sweep any parameter.
  * "No route found" by the retrosynthesis engine is recorded as a DISTINCT
    value (RouteStatus.NO_ROUTE_FOUND_TOOL_LIMITED) and is never conflated
    with "hard to synthesize". Search failure is a property of the tool and
    its template library, not of the molecule.
  * Tier 3 is a NEGATIVE screen (a proxy), not a positive prediction of
    crystallizability. See CRYSTALLIZABILITY NOTES below.
  * Every compound carries an audit trail: which rule fired, on what evidence.

Author: pipeline built against a DEVELOPMENT stand-in set; see
pipeline_validation.md for validation results and limitations.
"""

from __future__ import annotations

import json
import os
import time
import hashlib
import urllib.parse
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Optional

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, GraphDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")


# ===========================================================================
# CONFIG -- every tunable threshold lives here and nowhere else.
# ===========================================================================

CONFIG: dict[str, Any] = {
    # ---- provenance -------------------------------------------------------
    "dataset_label": "DEVELOPMENT",   # set to "HMDB_PREDICTED" for the real run

    # ---- PubChem / purchasability ----------------------------------------
    # Rate limit: PubChem permits max 5 requests/sec. We stay under it.
    "pubchem_min_interval_s": 0.25,        # -> <= 4 req/s
    "pubchem_timeout_s": 30,
    "pubchem_max_retries": 5,
    "pubchem_backoff_initial_s": 1.0,
    "pubchem_backoff_max_s": 30.0,
    "pubchem_cache_dir": "cache/pubchem",
    # The PUG-View source category that denotes a commercial supplier.
    "vendor_category_name": "Chemical Vendors",
    # Minimum number of DISTINCT vendor SourceNames for tier-1 "purchasable".
    # JUDGMENT CALL: 1 vendor is a weak signal (may be a catalog listing with
    # no stock); requiring >=2 distinct vendors is a modest robustness gain.
    "min_distinct_vendors": 1,
    # Vendors that specifically supply CERTIFIED ANALYTICAL / REFERENCE
    # STANDARDS as opposed to bulk research chemicals. Presence of one of
    # these is recorded separately as `has_reference_standard_vendor` and is
    # the operationally stricter definition of "purchasable" for MS
    # verification purposes. See pipeline_validation.md.
    "reference_standard_vendors": [
        "Sigma-Aldrich", "Merck", "Supelco", "Cerilliant",
        "Toronto Research Chemicals", "TRC", "LGC Standards",
        "Cayman Chemical", "Steraloids", "Avanti Polar Lipids",
        "USP", "U.S. Pharmacopeia", "Sigma", "MilliporeSigma",
        "Santa Cruz Biotechnology", "Tokyo Chemical Industry", "TCI",
        "Alfa Aesar", "Thermo Fisher Scientific", "Fisher Scientific",
        "CDN Isotopes", "Isosciences", "Carbosynth", "Biosynth",
    ],
    # If True, tier 1 requires a reference-standard vendor, not merely any
    # vendor. Default False (permissive) -- the permissive setting makes the
    # paper's residual-gap claim CONSERVATIVE, which is the honest direction.
    "require_reference_standard_vendor_for_tier1": False,

    # ---- retrosynthesis ---------------------------------------------------
    "retro_engine": "aizynthfinder",
    "retro_max_steps_for_short_route": 3,   # "few steps" == <= 3
    "retro_time_limit_s": 120,              # per-compound search cap
    "retro_iteration_limit": 100,
    "retro_n_workers": 32,

    # ---- crystallizability proxy -----------------------------------------
    # See CRYSTALLIZABILITY NOTES for provenance of each number.
    # NOTE: the two thresholds below are LITERATURE-DERIVED but their
    # citation is NOT YET VERIFIED against the primary article. Treat them as
    # judgment calls until checked. See CRYSTALLIZABILITY_NOTES.
    "cryst_mw_max": 350.0,            # literature-derived, citation UNVERIFIED
    "cryst_rotb_max": 4,              # literature-derived, citation UNVERIFIED
    "cryst_frac_rotb_max": 0.30,      # JUDGMENT CALL
    "cryst_tpsa_max": 140.0,          # JUDGMENT CALL (borrowed from Veber)
    "cryst_hbd_max": 5,               # JUDGMENT CALL (Lipinski-derived)
    "cryst_logp_min": -1.0,           # JUDGMENT CALL
    "cryst_abs_formal_charge_max": 1, # JUDGMENT CALL
    # A compound is flagged "likely hard to crystallize" when the number of
    # violated criteria reaches this count. A single violation is common among
    # perfectly crystallizable compounds, so we require a conjunction.
    # JUDGMENT CALL -- this is the single most important knob to sweep.
    "cryst_n_violations_for_hard_flag": 2,
    # Permanently-charged / zwitterionic species at physiological pH are
    # harder to obtain as neutral single crystals from solution.
    "physiological_ph": 7.4,

    # ---- misc -------------------------------------------------------------
    "random_seed": 20260804,
}


# ===========================================================================
# CRYSTALLIZABILITY NOTES  (verbatim provenance for each threshold)
# ===========================================================================
CRYSTALLIZABILITY_NOTES = """
cryst_mw_max = 350, cryst_rotb_max = 4
    LITERATURE-DERIVED, CITATION NOT YET VERIFIED AGAINST THE PRIMARY SOURCE.
    These two numbers come from the crystallization-propensity literature,
    where a combination of low molecular weight and low flexibility is
    reported as a crude predictor of ready crystallization from solution;
    the specific MW < 350 / fewer-than-4-rotatable-bonds pairing is the
    commonly quoted form. The relevant primary literature is the Wicker &
    Cooper line of work on crystallization propensity of drug-like molecules
    (see refs below).

    >>> ACTION REQUIRED BEFORE PUBLICATION: obtain the primary paper and
    >>> confirm the exact numbers, author, and year before citing them. This
    >>> attribution was assembled from secondary/search-summary material and
    >>> the specific author-year has NOT been checked against the article
    >>> itself. Do not put it in a manuscript uncited or unchecked.

    Candidate primary sources to verify against:
      - Wicker & Cooper, "Will it crystallise? Predicting crystallinity of
        molecular materials", CrystEngComm (2015).
      - "Predicting the Crystallization Propensity of Drug-Like Molecules"
        (J. Pharm. Sci. / ScienceDirect, ~2016).

    IMPORTANT CAVEAT, stated because the paper depends on honesty here:
    Cooper and co-workers (Chem. Sci. 2021, para-substituted benzoic acids)
    showed that conformational flexibility alone does NOT limit
    crystallizability -- a compound with MORE rotatable bonds nucleated
    FASTER than ones with fewer. Their conclusion is that it is the NATURE of
    the rotatable bonds (the conformational energy barrier relative to the
    ~70 kJ/mol activation energy for nucleation), not the COUNT, that
    matters. A rotatable-bond count is therefore a coarse surrogate for a
    quantity it does not actually measure. Better descriptors exist
    (nConf20, Wicker & Cooper, J. Chem. Inf. Model. 2017, "Beyond Rotatable
    Bond Counts: Capturing 3D Conformational Flexibility in a Single
    Descriptor"; and the Kier flexibility index), and published ML
    classifiers substantially outperform any threshold rule.

    >>> The specific accuracy figures previously quoted here (a held-out
    >>> error rate and a blind-screen accuracy) have been REMOVED: they came
    >>> from search-summary material and could not be verified against the
    >>> primary articles. Do not reintroduce them without checking the
    >>> papers directly.

    We use thresholds anyway because they are transparent and sweepable, but
    the paper should not claim this is a crystallizability PREDICTION.

cryst_frac_rotb_max = 0.30
    JUDGMENT CALL. Normalizes flexibility by molecular size so that large
    rigid polycyclics (steroid nuclei, bile acids) are not penalized purely
    for being large. No literature value anchors this number.

cryst_tpsa_max = 140
    JUDGMENT CALL, borrowed from Veber et al.'s oral-bioavailability rule,
    which is about membrane permeation and NOT about crystallization. Used
    here only as a proxy for "heavily polar, strongly solvated, prone to
    hydrate/solvate formation and oiling out". Weakly justified; sweep it.

cryst_hbd_max = 5
    JUDGMENT CALL, Lipinski-derived. Rationale is that many donors favour
    hydrate/solvate formation and conformational polymorphism. Note this cuts
    against itself: strong, directional H-bonding also frequently PROMOTES
    good crystal packing. The sign of this term is genuinely arguable.

cryst_logp_min = -1.0
    JUDGMENT CALL. Very hydrophilic species (sugars, polyamines, zwitterions)
    are hard to crystallize from aqueous media at trace abundance. No
    literature threshold.

cryst_abs_formal_charge_max = 1
    JUDGMENT CALL. Permanently charged species crystallize as salts whose
    counterion is not determined by the predicted structure.

cryst_n_violations_for_hard_flag = 2
    JUDGMENT CALL and the most consequential parameter in the whole proxy.
    Requiring 2+ violations rather than 1 prevents almost every polar
    metabolite from being flagged. THE HEADLINE TIER-4 FRACTION IS SENSITIVE
    TO THIS NUMBER -- report the sweep, not a single value.

MicroED-specific note:
    MicroED works from sub-micron crystals -- orders of magnitude smaller than
    single-crystal X-ray requires -- which is why tier 3 is plausible at all at
    trace abundance: the technique routinely obtains structures from powders
    that appear amorphous by eye, and often eliminates crystallization
    screening because compounds crystallize adequately during ordinary
    purification. That said, throughput, instrument access and cost remain the
    practical limits at discovery scale, so "MicroED-plausible" is a statement
    about physical possibility, not routine feasibility.

    >>> A specific nanometre crystal-size range quoted here earlier has been
    >>> removed as unverified. If the paper needs a numeric size threshold,
    >>> take it from a named MicroED methods paper and cite that paper.
"""


# ===========================================================================
# Enums / data model
# ===========================================================================

class Tier(str, Enum):
    PURCHASABLE = "1_purchasable"
    SHORT_ROUTE = "2_short_route_synthesizable"
    MICROED = "3_microed_plausible"
    RESIDUAL = "4_residual"
    UNASSIGNED = "0_unassigned"


class RouteStatus(str, Enum):
    """Retrosynthesis outcome.

    NO_ROUTE_FOUND_TOOL_LIMITED is deliberately verbose. It means: this
    engine, with this template library and this stock file, did not find a
    route within the search budget. It does NOT mean the molecule is hard to
    synthesize, and it must never be reported as such.
    """
    SOLVED_SHORT = "solved_short_route"
    SOLVED_LONG = "solved_long_route"
    NO_ROUTE_FOUND_TOOL_LIMITED = "no_route_found_TOOL_LIMITED"
    ALREADY_IN_STOCK = "already_in_stock"
    ERROR = "engine_error"
    NOT_RUN = "not_run"


class PurchaseStatus(str, Enum):
    VENDOR_LISTED = "vendor_listed"
    REFERENCE_STANDARD_VENDOR = "reference_standard_vendor_listed"
    NO_VENDOR = "no_vendor_annotation"
    NO_CID = "structure_not_in_pubchem"
    ERROR = "lookup_error"
    NOT_RUN = "not_run"


# ===========================================================================
# STEP 2 -- Standardization
# ===========================================================================

_UNCHARGER = rdMolStandardize.Uncharger()
_TAUTOMER = rdMolStandardize.TautomerEnumerator()


def standardize(smiles: str) -> tuple[Optional[Chem.Mol], dict]:
    """Parse -> sanitize -> largest-fragment -> neutralize -> canonical tautomer.

    Returns (mol, record). `mol` is None on failure; record['std_status'] is
    one of ok / parse_fail / standardize_fail / empty_after_std. The count of
    non-ok rows is a reported quantity -- for the real HMDB run it tells us
    what fraction of predicted structures we could not even read.
    """
    rec = {"std_status": "ok", "std_note": ""}
    if not smiles or not isinstance(smiles, str):
        rec["std_status"] = "parse_fail"
        rec["std_note"] = "empty or non-string input"
        return None, rec
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        rec["std_status"] = "parse_fail"
        rec["std_note"] = "RDKit MolFromSmiles returned None"
        return None, rec
    try:
        mol = rdMolStandardize.Cleanup(mol)
        n_frag = len(Chem.GetMolFrags(mol))
        if n_frag > 1:
            mol = rdMolStandardize.FragmentParent(mol)
            rec["std_note"] += f"stripped {n_frag - 1} salt/solvate fragment(s);"
        mol = _UNCHARGER.uncharge(mol)
        try:
            canon = _TAUTOMER.Canonicalize(mol)
            if canon is not None:
                mol = canon
            else:
                rec["std_note"] += "tautomer_canonicalize_returned_none;"
        except Exception as exc:                      # noqa: BLE001
            # Tautomer enumeration can blow up on large/unusual systems.
            # Degrade gracefully rather than dropping the compound.
            rec["std_note"] += f"tautomer_skipped({type(exc).__name__});"
        Chem.SanitizeMol(mol)
    except Exception as exc:                          # noqa: BLE001
        rec["std_status"] = "standardize_fail"
        rec["std_note"] += f"{type(exc).__name__}: {exc}"
        return None, rec
    if mol.GetNumHeavyAtoms() == 0:
        rec["std_status"] = "empty_after_std"
        rec["std_note"] += "no heavy atoms remain"
        return None, rec
    return mol, rec


def standardize_frame(records: list[dict], smiles_key: str = "smiles_input") -> list[dict]:
    """Standardize a list of dicts in place-ish; returns new list with std fields.

    Adds: smiles_std, inchikey_std, inchikey14 (skeleton key), std_status,
    std_note. Does NOT deduplicate -- call deduplicate() for that, so the
    parse-failure count is preserved separately from the dedup count.
    """
    out = []
    for rec in records:
        mol, srec = standardize(rec.get(smiles_key))
        row = dict(rec)
        row.update(srec)
        if mol is not None:
            row["smiles_std"] = Chem.MolToSmiles(mol)
            row["inchikey_std"] = Chem.MolToInchiKey(mol)
            row["inchikey14"] = row["inchikey_std"].split("-")[0]
        else:
            row["smiles_std"] = None
            row["inchikey_std"] = None
            row["inchikey14"] = None
        out.append(row)
    return out


def deduplicate(records: list[dict], key: str = "inchikey_std") -> tuple[list[dict], dict]:
    """Collapse rows sharing a standardized InChIKey. Returns (rows, stats).

    Merged names are retained in `merged_query_names` so the audit trail can
    show which input identifiers collapsed together.
    """
    seen: dict[str, dict] = {}
    order: list[str] = []
    n_dropped = 0
    unparsed = [r for r in records if not r.get(key)]
    for rec in records:
        k = rec.get(key)
        if not k:
            continue
        if k not in seen:
            seen[k] = dict(rec)
            seen[k]["merged_query_names"] = [rec.get("query_name")]
            order.append(k)
        else:
            seen[k]["merged_query_names"].append(rec.get("query_name"))
            n_dropped += 1
    rows = [seen[k] for k in order]
    stats = {
        "n_input": len(records),
        "n_unparsed": len(unparsed),
        "n_unique": len(rows),
        "n_duplicates_collapsed": n_dropped,
    }
    return rows, stats


# ===========================================================================
# STEP 5 -- Crystallizability proxy descriptors
# ===========================================================================

def compute_descriptors(mol: Chem.Mol) -> dict:
    """RDKit physicochemical descriptors relevant to crystallization propensity."""
    n_rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    # Fraction of rotatable bonds among all non-ring single bonds between
    # heavy atoms -- a size-normalized flexibility measure.
    n_bonds = mol.GetNumBonds()
    desc = {
        "mw": Descriptors.MolWt(mol),
        "exact_mass": Descriptors.ExactMolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "n_rot_bonds": n_rotb,
        "frac_rot_bonds": (n_rotb / n_bonds) if n_bonds else 0.0,
        "n_heavy_atoms": n_heavy,
        "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "n_rings": rdMolDescriptors.CalcNumRings(mol),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "n_stereocenters": len(Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False)),
        "formal_charge_neutral_form": Chem.GetFormalCharge(mol),
        # Kier molecular flexibility index (Phi) = kappa1 * kappa2 / n_heavy.
        # An established 2D flexibility descriptor, cited alongside rotatable
        # bond count in the crystallization-propensity literature.
        "kier_phi": (GraphDescriptors.Kappa1(mol) * GraphDescriptors.Kappa2(mol)
                     / n_heavy) if n_heavy else 0.0,
        "n_amide_bonds": rdMolDescriptors.CalcNumAmideBonds(mol),
        "n_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
    }
    desc["charge_at_ph74"] = estimate_charge_at_ph(mol)
    return desc


# Coarse ionization model. Substructure -> charge contribution at pH 7.4.
# This is a RULE OF THUMB, not a pKa calculation. Documented as such.
_ACID_PATTERNS = [
    ("carboxylic_acid", Chem.MolFromSmarts("[CX3](=O)[OX2H1]"), -1),
    ("sulfonic_acid", Chem.MolFromSmarts("[SX4](=O)(=O)[OX2H1]"), -1),
    ("sulfate_ester", Chem.MolFromSmarts("[OX2][SX4](=O)(=O)[OX2H1]"), -1),
    ("phosphate", Chem.MolFromSmarts("[PX4](=O)([OX2H1])[OX2H1]"), -2),
]
_BASE_PATTERNS = [
    ("aliphatic_primary_amine", Chem.MolFromSmarts("[NX3;H2;!$(NC=O);!$(Nc)]"), +1),
    ("aliphatic_secondary_amine", Chem.MolFromSmarts("[NX3;H1;!$(NC=O);!$(Nc);!R0]"), +1),
    ("guanidine", Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3]"), +1),
    ("quaternary_N", Chem.MolFromSmarts("[NX4+]"), +1),
]


def estimate_charge_at_ph(mol: Chem.Mol, ph: float | None = None) -> int:
    """Very coarse net-charge estimate at physiological pH.

    Counts strongly acidic groups (assumed deprotonated at pH 7.4) and
    strongly basic groups (assumed protonated). This is a SMARTS rule of
    thumb -- it has no pKa model behind it and will be wrong for anything
    with an unusual electronic environment. Used only to identify obviously
    zwitterionic / permanently charged species.
    """
    charge = 0
    for _name, patt, delta in _ACID_PATTERNS:
        if patt is None:
            continue
        charge += delta * len(mol.GetSubstructMatches(patt))
    for _name, patt, delta in _BASE_PATTERNS:
        if patt is None:
            continue
        charge += delta * len(mol.GetSubstructMatches(patt))
    return charge


def crystallizability_flag(desc: dict, config: dict = CONFIG) -> dict:
    """Threshold-based 'likely hard to crystallize' proxy.

    Returns a dict with the boolean flag, the count of violated criteria, and
    the explicit list of which criteria were violated (the audit trail).

    THIS IS A PROXY, NOT A PREDICTION. See CRYSTALLIZABILITY_NOTES.
    """
    violations: list[str] = []

    if desc["mw"] > config["cryst_mw_max"]:
        violations.append(f"mw>{config['cryst_mw_max']}({desc['mw']:.1f})")
    if desc["n_rot_bonds"] > config["cryst_rotb_max"]:
        violations.append(
            f"n_rot_bonds>{config['cryst_rotb_max']}({desc['n_rot_bonds']})")
    if desc["frac_rot_bonds"] > config["cryst_frac_rotb_max"]:
        violations.append(
            f"frac_rot_bonds>{config['cryst_frac_rotb_max']}"
            f"({desc['frac_rot_bonds']:.2f})")
    if desc["tpsa"] > config["cryst_tpsa_max"]:
        violations.append(f"tpsa>{config['cryst_tpsa_max']}({desc['tpsa']:.1f})")
    if desc["hbd"] > config["cryst_hbd_max"]:
        violations.append(f"hbd>{config['cryst_hbd_max']}({desc['hbd']})")
    if desc["clogp"] < config["cryst_logp_min"]:
        violations.append(
            f"clogp<{config['cryst_logp_min']}({desc['clogp']:.2f})")
    if abs(desc["charge_at_ph74"]) > config["cryst_abs_formal_charge_max"]:
        violations.append(
            f"|charge_pH7.4|>{config['cryst_abs_formal_charge_max']}"
            f"({desc['charge_at_ph74']})")

    n_viol = len(violations)
    hard = n_viol >= config["cryst_n_violations_for_hard_flag"]
    return {
        "cryst_n_violations": n_viol,
        "cryst_violations": ";".join(violations) if violations else "",
        "likely_hard_to_crystallize": hard,
    }


# ===========================================================================
# STEP 3 -- Purchasability via PubChem
# ===========================================================================
#
# OPERATIONAL DEFINITION OF "PURCHASABLE" -- read this before using tier 1.
#
# What we actually measure: whether PubChem's PUG-View "categorized sources"
# endpoint lists >= min_distinct_vendors distinct suppliers under the
# "Chemical Vendors" source category for the CID that our standardized
# structure resolves to.
#
# What that is NOT:
#   (a) It is not evidence that an AUTHENTIC ANALYTICAL STANDARD exists. Most
#       "Chemical Vendors" entries are catalog listings for research-grade
#       material with no certificate of analysis, no isotopic-purity
#       statement, and no quantitative reference value. For LC-MS/MS
#       verification you need a certified standard, which is a much smaller
#       set of suppliers.
#   (b) It is not evidence of CURRENT STOCK. Many catalog entries are
#       make-to-order or list compounds the vendor has never synthesized.
#       PubChem does not carry stock or price.
#   (c) It is not STEREOCHEMICALLY EXACT. A CID match on a flat structure can
#       return a vendor listing for a different stereoisomer or a racemate,
#       which is not an acceptable standard for a stereochemically defined
#       metabolite.
#   (d) Vendor annotations are a snapshot; PubChem source records lag reality
#       in both directions.
#
# Therefore this pipeline reports TWO purchasability columns:
#   purchase_status              -- any vendor at all (permissive)
#   has_reference_standard_vendor -- at least one supplier from a curated list
#                                    of analytical-standard houses (strict)
# The permissive definition inflates tier 1 and therefore SHRINKS the
# headline tier-4 residual. That is the conservative direction for a paper
# claiming a verification gap exists, which is why it is the default.

_LAST_REQUEST = [0.0]


def _cache_path(url: str, config: dict) -> str:
    os.makedirs(config["pubchem_cache_dir"], exist_ok=True)
    key = hashlib.sha256(url.encode()).hexdigest()[:32]
    return os.path.join(config["pubchem_cache_dir"], key + ".json")


def pubchem_get(url: str, config: dict = CONFIG) -> tuple[int, str]:
    """Rate-limited, disk-cached GET against PubChem.

    Honours the documented 5 req/s cap (we use 4), caches every terminal
    response (including 404, which is a real answer), and backs off
    exponentially on 503/429/500.
    """
    path = _cache_path(url, config)
    if os.path.exists(path):
        try:
            data = json.load(open(path))
            return data["status"], data["text"]
        except Exception:                              # noqa: BLE001
            os.remove(path)

    backoff = config["pubchem_backoff_initial_s"]
    for _attempt in range(config["pubchem_max_retries"]):
        elapsed = time.time() - _LAST_REQUEST[0]
        if elapsed < config["pubchem_min_interval_s"]:
            time.sleep(config["pubchem_min_interval_s"] - elapsed)
        try:
            resp = requests.get(url, timeout=config["pubchem_timeout_s"])
            _LAST_REQUEST[0] = time.time()
        except Exception:                              # noqa: BLE001
            time.sleep(backoff)
            backoff = min(backoff * 2, config["pubchem_backoff_max_s"])
            continue
        if resp.status_code in (429, 500, 503):
            # Server asked us to slow down -- throttle harder.
            time.sleep(backoff)
            backoff = min(backoff * 2, config["pubchem_backoff_max_s"])
            continue
        json.dump({"status": resp.status_code, "text": resp.text, "url": url},
                  open(path, "w"))
        return resp.status_code, resp.text
    return 0, ""


def resolve_cid(inchikey: str | None, smiles: str | None,
                config: dict = CONFIG) -> tuple[Optional[int], str]:
    """Resolve a structure to a PubChem CID. InChIKey first, then SMILES.

    Returns (cid, method). InChIKey is tried first because it is an exact
    structure match including stereochemistry; SMILES fallback is looser.
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    if inchikey:
        url = f"{base}/inchikey/{urllib.parse.quote(inchikey)}/cids/JSON"
        status, text = pubchem_get(url, config)
        if status == 200:
            try:
                cids = json.loads(text)["IdentifierList"]["CID"]
                if cids:
                    return int(cids[0]), "inchikey"
            except Exception:                          # noqa: BLE001
                pass
    if smiles:
        url = f"{base}/smiles/{urllib.parse.quote(smiles)}/cids/JSON"
        status, text = pubchem_get(url, config)
        if status == 200:
            try:
                cids = json.loads(text)["IdentifierList"]["CID"]
                if cids and int(cids[0]) != 0:
                    return int(cids[0]), "smiles"
            except Exception:                          # noqa: BLE001
                pass
    return None, "unresolved"


def get_vendor_info(cid: int, config: dict = CONFIG) -> dict:
    """Fetch categorized sources for a CID and extract vendor annotation."""
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/"
           f"compound/{cid}/JSON")
    status, text = pubchem_get(url, config)
    info = {
        "n_vendor_records": 0,
        "n_distinct_vendors": 0,
        "vendor_names_sample": "",
        "has_reference_standard_vendor": False,
        "reference_standard_vendors_found": "",
        "vendor_lookup_http_status": status,
    }
    if status != 200:
        return info
    try:
        cats = json.loads(text)["SourceCategories"]["Categories"]
    except Exception:                                  # noqa: BLE001
        return info

    vendor_names: list[str] = []
    for cat in cats:
        if cat.get("Category") != config["vendor_category_name"]:
            continue
        for src in cat.get("Sources", []):
            name = src.get("SourceName")
            if name:
                vendor_names.append(name)

    distinct = sorted(set(vendor_names))
    ref_lower = [v.lower() for v in config["reference_standard_vendors"]]
    found_ref = [v for v in distinct
                 if any(r in v.lower() or v.lower() in r for r in ref_lower)]

    info["n_vendor_records"] = len(vendor_names)
    info["n_distinct_vendors"] = len(distinct)
    info["vendor_names_sample"] = "|".join(distinct[:8])
    info["has_reference_standard_vendor"] = bool(found_ref)
    info["reference_standard_vendors_found"] = "|".join(found_ref[:6])
    return info


def assess_purchasability(inchikey: str | None, smiles: str | None,
                          config: dict = CONFIG) -> dict:
    """Full purchasability assessment for one structure."""
    result = {
        "cid": None,
        "cid_resolution_method": "unresolved",
        "purchase_status": PurchaseStatus.NOT_RUN.value,
        "n_vendor_records": 0,
        "n_distinct_vendors": 0,
        "vendor_names_sample": "",
        "has_reference_standard_vendor": False,
        "reference_standard_vendors_found": "",
    }
    try:
        cid, method = resolve_cid(inchikey, smiles, config)
    except Exception as exc:                           # noqa: BLE001
        result["purchase_status"] = PurchaseStatus.ERROR.value
        result["cid_resolution_method"] = f"error:{type(exc).__name__}"
        return result

    result["cid"] = cid
    result["cid_resolution_method"] = method
    if cid is None:
        result["purchase_status"] = PurchaseStatus.NO_CID.value
        return result

    info = get_vendor_info(cid, config)
    result.update({k: v for k, v in info.items()
                   if k != "vendor_lookup_http_status"})

    if info["n_distinct_vendors"] >= config["min_distinct_vendors"]:
        result["purchase_status"] = (
            PurchaseStatus.REFERENCE_STANDARD_VENDOR.value
            if info["has_reference_standard_vendor"]
            else PurchaseStatus.VENDOR_LISTED.value)
    else:
        result["purchase_status"] = PurchaseStatus.NO_VENDOR.value
    return result


def is_purchasable(row: dict, config: dict = CONFIG) -> bool:
    """Apply the configured purchasability definition to an assessed row.

    IMPORTANT: this re-derives the decision from the RAW evidence columns
    (n_distinct_vendors, has_reference_standard_vendor) rather than reading
    the frozen `purchase_status` string. That is deliberate -- it makes both
    `min_distinct_vendors` and `require_reference_standard_vendor_for_tier1`
    genuinely sweepable in a post-hoc sensitivity analysis WITHOUT re-hitting
    the PubChem API. Reading the cached status string instead would silently
    pin those two thresholds to whatever they were at lookup time.
    """
    status = row.get("purchase_status")
    # Terminal states carrying no vendor evidence can never be purchasable.
    if status in (PurchaseStatus.NO_CID.value, PurchaseStatus.ERROR.value,
                  PurchaseStatus.NOT_RUN.value):
        return False
    n_vendors = row.get("n_distinct_vendors") or 0
    if n_vendors < config["min_distinct_vendors"]:
        return False
    if config["require_reference_standard_vendor_for_tier1"]:
        return bool(row.get("has_reference_standard_vendor"))
    return True


# ===========================================================================
# STEP 4 -- Retrosynthesis (AiZynthFinder, local)
# ===========================================================================
#
# ENGINE PROVENANCE (record this in any results table):
#   aizynthfinder      4.4.1        (pip, MIT licence)
#   expansion policy   uspto_model.onnx + uspto_templates.csv.gz
#   ringbreaker policy uspto_ringbreaker_model.onnx (+ templates)
#   filter policy      uspto_filter_model.onnx
#   stock              zinc_stock.hdf5 -- 17,422,831 InChIKeys
#   obtained via AiZynthFinder's own `download_public_data` mechanism.
#
# WHAT "IN STOCK" MEANS HERE: membership in the bundled ZINC stock file,
# which is a list of 17.4M InChIKeys of ZINC-catalogued compounds. It is NOT
# a real-time vendor stock check and carries no price, purity, or lead-time
# information. A building block "in stock" by this definition may be
# unobtainable in practice, and conversely a genuinely purchasable reagent
# absent from this snapshot will be treated as not-in-stock.
#
# WHY "no route found" IS TOOL-LIMITED, EMPIRICALLY:
#   In our development-set benchmark, AiZynthFinder failed to find a route to
#   TRIMETHYLAMINE -- a bulk commodity chemical produced industrially at
#   megatonne scale. It also failed on essentially every bile acid
#   (chenodeoxycholic, hyodeoxycholic, allocholic, isodeoxycholic, ...).
#   These are not hard molecules; they are molecules whose real routes
#   (steroid semi-synthesis from natural feedstocks, C1 amination chemistry)
#   are absent or under-represented in the USPTO template library. A
#   USPTO-trained template policy encodes medicinal-chemistry reaction space,
#   which systematically under-covers the natural-product / steroid / sugar
#   chemistry that dominates gut microbiome metabolites.
#   => NO_ROUTE_FOUND_TOOL_LIMITED is a statement about the template library,
#      never about the molecule. Do not report it as "not synthesizable".


def classify_route_result(rec: dict, config: dict = CONFIG) -> tuple[str, str]:
    """Map a raw AiZynthFinder result record to (RouteStatus, evidence string).

    `rec` is expected to carry: retro_run_status, target_in_stock, is_solved,
    n_steps_top, n_routes -- as produced by the parallel runner.
    """
    if rec.get("retro_run_status") == "error":
        return (RouteStatus.ERROR.value,
                f"engine error: {rec.get('retro_error', '')[:120]}")
    if rec.get("retro_run_status") in (None, "not_run"):
        return RouteStatus.NOT_RUN.value, "retrosynthesis not executed"

    if rec.get("target_in_stock"):
        return (RouteStatus.ALREADY_IN_STOCK.value,
                "target itself is present in the ZINC stock file; "
                "retrosynthesis is vacuous for this compound")

    if not rec.get("is_solved"):
        return (RouteStatus.NO_ROUTE_FOUND_TOOL_LIMITED.value,
                f"no solved route within {config['retro_time_limit_s']}s / "
                f"{config['retro_iteration_limit']} iterations using the USPTO "
                f"template policy; TOOL-LIMITED, not proven infeasible")

    n_steps = rec.get("n_steps_top")
    if n_steps is None:
        return (RouteStatus.SOLVED_LONG.value,
                "solved but step count unavailable")
    if n_steps <= config["retro_max_steps_for_short_route"]:
        return (RouteStatus.SOLVED_SHORT.value,
                f"solved route of {n_steps} step(s) to in-stock building "
                f"blocks (<= {config['retro_max_steps_for_short_route']})")
    return (RouteStatus.SOLVED_LONG.value,
            f"solved route of {n_steps} step(s), exceeds short-route "
            f"threshold of {config['retro_max_steps_for_short_route']}")


def is_short_route(route_status: str) -> bool:
    """Tier-2 predicate. ALREADY_IN_STOCK also qualifies (it is buyable as a
    building block even if PubChem showed no vendor annotation)."""
    return route_status in (RouteStatus.SOLVED_SHORT.value,
                            RouteStatus.ALREADY_IN_STOCK.value)


# ===========================================================================
# STEP 6 -- Tier assignment
# ===========================================================================

def assign_tier(row: dict, config: dict = CONFIG) -> dict:
    """Assign one compound to a verification tier, with an audit trail.

    PRECEDENCE (first matching rule wins):
      1. purchasable                 -- vendor annotation satisfies the
                                        configured definition
      2. short_route_synthesizable   -- retro engine solved a route of
                                        <= retro_max_steps_for_short_route
                                        steps, or the target is itself an
                                        in-stock building block
      3. microed_plausible           -- NOT flagged likely-hard-to-crystallize
      4. residual                    -- everything else

    Returns dict with: tier, tier_rule_fired, tier_evidence.
    """
    trail: list[str] = []

    # ---- Rule 1: purchasable ---------------------------------------------
    purch = is_purchasable(row, config)
    defn = ("reference-standard vendor required"
            if config["require_reference_standard_vendor_for_tier1"]
            else f">= {config['min_distinct_vendors']} distinct PubChem vendor(s)")
    trail.append(
        f"R1[purchasable={purch}] status={row.get('purchase_status')} "
        f"n_distinct_vendors={row.get('n_distinct_vendors')} "
        f"ref_std_vendor={row.get('has_reference_standard_vendor')} "
        f"defn=({defn})")
    if purch:
        return {"tier": Tier.PURCHASABLE.value,
                "tier_rule_fired": "R1_vendor_annotation",
                "tier_evidence": " | ".join(trail)}

    # ---- Rule 2: short-route synthesizable -------------------------------
    rstatus = row.get("route_status", RouteStatus.NOT_RUN.value)
    short = is_short_route(rstatus)
    trail.append(
        f"R2[short_route={short}] route_status={rstatus} "
        f"n_steps_top={row.get('n_steps_top')} "
        f"max_steps={config['retro_max_steps_for_short_route']} "
        f"evidence={row.get('route_evidence', '')}")
    if short:
        return {"tier": Tier.SHORT_ROUTE.value,
                "tier_rule_fired": "R2_retrosynthesis_short_route",
                "tier_evidence": " | ".join(trail)}

    # ---- Rule 3: MicroED-plausible ---------------------------------------
    hard = bool(row.get("likely_hard_to_crystallize"))
    trail.append(
        f"R3[microed_plausible={not hard}] "
        f"cryst_n_violations={row.get('cryst_n_violations')}"
        f"(threshold {config['cryst_n_violations_for_hard_flag']}) "
        f"violations=[{row.get('cryst_violations', '')}]")
    if not hard:
        return {"tier": Tier.MICROED.value,
                "tier_rule_fired": "R3_crystallizability_proxy_pass",
                "tier_evidence": " | ".join(trail)}

    # ---- Rule 4: residual -------------------------------------------------
    trail.append("R4[residual] failed R1, R2 and R3")
    return {"tier": Tier.RESIDUAL.value,
            "tier_rule_fired": "R4_residual",
            "tier_evidence": " | ".join(trail)}


def tier_fraction_ci(n_tier: int, n_total: int, z: float = 1.96) -> tuple[float, float, float]:
    """Wilson score interval for a tier fraction. Returns (p, lo, hi).

    Used for the sampling recommendation: if the full set is not run, the
    reported tier fractions need an interval, and Wilson behaves correctly
    near p=0 or p=1 where the normal approximation does not.
    """
    if n_total == 0:
        return 0.0, 0.0, 0.0
    p = n_tier / n_total
    denom = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denom
    half = (z * ((p * (1 - p) / n_total + z**2 / (4 * n_total**2)) ** 0.5)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


def required_sample_size(margin: float = 0.03, p: float = 0.5, z: float = 1.96) -> int:
    """n for a target +/- margin on a proportion (worst case p=0.5)."""
    import math
    return math.ceil(z**2 * p * (1 - p) / margin**2)
