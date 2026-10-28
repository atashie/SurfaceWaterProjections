#!/usr/bin/env python
"""Stage the HydroShare resource folders for the version-3 products (public name "version 2": the first version-2
folders were never deposited) after the qualifier-mask + climate-frame rerun (plan
docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md, Phase H; 2026-10-28).

    signatures --run-dir RUN --prefix PREFIX --window 1993|1980 --rerun-dir R --staged-v2 <v2 staged R1|R2> --out OUT
    inputs     --rerun-dir R --staged-v2 <v2 staged R3> --base-parquet <February compile> --out OUT
    check      --out OUT --kind signatures|inputs [--window W]

Local staging only (nothing is uploaded). `signatures` copies and renames the run's artefacts to the deposit names of
2026 (summary CSV, annual parquet, provenance JSON, explorer + `_annual/` sidecar with the embedded folder name
rewritten, the validation record), adds the frozen selection manifest, the would-qualify list, the change census,
the qualification audit and the port check, and rewrites the version-2 README's counts and conventions. `inputs`
stages resource 3: the immutable February compile beside the corrected table, the correction ledger, the raw NWIS
responses, the three-dataset control report, the acceptance oracle, the selection manifests, the executable
reconstruction recipe, the updated metadata, README, abstract and dictionary. `check` verifies that every required
file is present (PAR1 footers, the explorer's embedded sidecar name) and writes a sha256 manifest. Every number in the
README and run notes is read from the files (census classes / gates, oracle, ledger, timing JSON), never typed here.
"""
import argparse, hashlib, json, os, re, shutil, subprocess, sys, tarfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def sha(p, blk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(blk), b""):
            h.update(b)
    return h.hexdigest()


def copy(src, dst):
    if os.path.isdir(src):
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)
    return os.path.basename(src), os.path.basename(dst)


def rename_sidecar_dir(html_path, old_dir, new_dir):
    with open(html_path, "rb") as fh:
        t = fh.read()
    old = f'"dir":"{old_dir}"'.encode(); new = f'"dir":"{new_dir}"'.encode()
    n = t.count(old)
    if n != 1:
        raise SystemExit(f"{html_path}: expected exactly one embedded sidecar dir {old!r}, found {n}")
    with open(html_path, "wb") as fh:
        fh.write(t.replace(old, new))


def must_replace(t, old, new, label):
    if t.count(old) != 1:
        raise SystemExit(f"README anchor not found exactly once ({t.count(old)}): {label}")
    return t.replace(old, new)


# ----------------------------------------------------------------------------------------------- signatures (R1 / R2)
def stage_signatures(a):
    W = a.window; tag = f"wy{W}_2025"; WU = f"{W}_2025"; Wd = f"{W}-2025"
    os.makedirs(a.out, exist_ok=True)
    run, p, R = a.run_dir, a.prefix, a.rerun_dir
    manifest = []
    csv_src = os.path.join(run, f"{p}_signatures.csv"); ann_src = os.path.join(run, f"{p}_signatures_annual.parquet")
    tim_src = os.path.join(run, f"{p}_timing.json"); mrg_src = os.path.join(run, f"{p}_merge.json")
    for src, dst in [(csv_src, f"hisss_signatures_{tag}.csv"), (ann_src, f"hisss_signatures_annual_{tag}.parquet"),
                     (tim_src, f"hisss_run_provenance_{tag}.json"), (mrg_src, f"hisss_run_merge_{tag}.json")]:
        if not os.path.exists(src):
            raise SystemExit(f"missing {src}")
        manifest.append(copy(src, os.path.join(a.out, dst)))
    expl = os.path.join(run, f"signature_explorer_{WU}_60pct_climate_v3.html")
    if not os.path.exists(expl):
        raise SystemExit(f"missing explorer {expl}")
    expl_dst = os.path.join(a.out, f"hisss_signature_explorer_{tag}.html"); manifest.append(copy(expl, expl_dst))
    side = expl[:-5] + "_annual"; side_dst = f"hisss_signature_explorer_{tag}_annual"
    if not os.path.isdir(side):
        raise SystemExit(f"missing sidecar folder {side}")
    manifest.append(copy(side, os.path.join(a.out, side_dst))); rename_sidecar_dir(expl_dst, os.path.basename(side), side_dst)
    # validation record (run folder + rerun folder)
    cen = os.path.join(R, f"census_{WU}"); e7 = os.path.join(R, f"e7_{WU}")
    for src, dst in [(os.path.join(run, "validation_gates.txt"), f"hisss_validation_gates_{tag}.txt"),
                     (os.path.join(run, "validate_annual_values_floor20.log"), f"hisss_validation_annual_values_{tag}.txt"),
                     (os.path.join(run, "audit_qualification.log"), f"hisss_validation_qualification_audit_{tag}.txt"),
                     (os.path.join(cen, f"change_census_v3_vs_v2_{Wd}.md"), f"hisss_validation_change_census_{tag}.md"),
                     (os.path.join(cen, f"gates_v3_vs_v2_{Wd}.json"), f"hisss_validation_change_census_gates_{tag}.json"),
                     (os.path.join(cen, f"classes_v3_vs_v2_{Wd}.csv"), f"hisss_validation_change_classes_{tag}.csv"),
                     (os.path.join(e7, "e7_annual.txt"), f"hisss_validation_port_check_annual_{tag}.txt"),
                     (os.path.join(e7, "e7_summary.txt"), f"hisss_validation_port_check_summary_{tag}.txt"),
                     (os.path.join(R, "manifests", f"selection_{Wd}.txt"), f"hisss_gage_selection_{tag}.txt")]:
        if not os.path.exists(src):
            raise SystemExit(f"missing {src}")
        manifest.append(copy(src, os.path.join(a.out, dst)))
    # the eight batch console logs, concatenated (the batch runner keeps one per batch folder)
    runs = json.load(open(os.path.join(run, "batch_runs.json")))
    with open(os.path.join(a.out, f"hisss_run_log_{tag}.log"), "w", encoding="utf-8") as fh:
        for k in sorted(runs, key=int):
            fh.write(f"===== batch {k} of {len(runs)} — {run}/b{k}/run.log =====\n")
            fh.write(open(os.path.join(run, f"b{k}", "run.log"), encoding="utf-8", errors="replace").read()); fh.write("\n")
    manifest.append((f"b1..b{len(runs)}/run.log", f"hisss_run_log_{tag}.log"))
    # dictionary (the version-2 staged file, already extended with the climate family; the stale flagged_for_high_na
    # definition of the 2026 deposit replaced) + categories unchanged
    dd = pd.read_csv(os.path.join(a.staged_v2, "hisss_data_dictionary.csv"), dtype=str, keep_default_na=False)
    m_hn = dd["column_name"] == "flagged_for_high_na"; assert m_hn.sum() == 1
    dd.loc[m_hn, "definition"] = ("Automated QA flag: TRUE when more than 30% of this gage's signature-derived fields (every column except the "
                                  "metadata and the other flags) are NA. Regenerated in this version with the shared Julia / Python / R definition; "
                                  "the 2026 deposit had computed it over the metadata columns only (that known issue is resolved here).")
    m_ice = dd["column_name"] == "ice_affected_days_total"; assert m_ice.sum() == 1
    n_ice_nonzero = int((pd.to_numeric(pd.read_csv(os.path.join(a.out, f"hisss_signatures_{tag}.csv"), usecols=["ice_affected_days_total"])["ice_affected_days_total"], errors="coerce").fillna(0) != 0).sum())
    dd.loc[m_ice, "definition"] = (f"Diagnostic intended as the total of ice-affected missing days across the record (from agency quality codes). KNOWN ISSUE: "
                                   f"in this version the canonical implementation emits 0 for every gage ({n_ice_nonzero} non-zero values in this table), so a zero does NOT "
                                   "establish the absence of ice-affected days — use the daily `flag` column of the input resource (code `Ice`, HYDAT `B`) instead.")
    dd.to_csv(os.path.join(a.out, "hisss_data_dictionary.csv"), index=False); manifest.append(("hisss_data_dictionary.csv", "hisss_data_dictionary.csv"))
    manifest.append(copy(os.path.join(a.staged_v2, "hisss_signature_categories.csv"), os.path.join(a.out, "hisss_signature_categories.csv")))
    # the would-qualify gages of this window (listed, not added — user decision 2026-10-23)
    wq = pd.read_csv(os.path.join(R, "input", "oracle", "would_qualify_outside_products.csv"), dtype={"gage_id": str})
    wq = wq[(wq["window"] == Wd) & (~wq["in_product"].astype(str).str.lower().isin(["true", "1"])) & (wq["final_qualifies"].astype(str).str.lower().isin(["true", "1"]))]
    wq_cols = ["gage_id", "window", "baseline_n_valid", "final_n_valid", "n_gained", "gained_years", "final_total_possible", "final_frac", "n_restored_days_in_window"]
    wq[wq_cols].to_csv(os.path.join(a.out, f"hisss_would_qualify_not_added_{tag}.csv"), index=False)
    n_wq = len(wq)

    # ---- numbers from the files
    summ = pd.read_csv(os.path.join(a.out, f"hisss_signatures_{tag}.csv"), dtype={"gage_id": str}, low_memory=False)
    n_gages, n_cols = summ.shape
    ann = pd.read_parquet(os.path.join(a.out, f"hisss_signatures_annual_{tag}.parquet"), columns=["signature"])
    n_rows = len(ann); n_sigs = ann.signature.nunique()
    tj = json.load(open(tim_src)); mj = json.load(open(mrg_src)); prov = tj.get("provenance", {})
    clim_na = int(summ["prcp_annual_mean"].isna().sum()); unnorm = int(summ["area_normalized"].astype(str).str.upper().isin(["FALSE", "0"]).sum())
    cl = pd.read_csv(os.path.join(cen, f"classes_v3_vs_v2_{Wd}.csv"), dtype={"gage_id": str})
    n_U = int((cl["class"] == "U").sum()); n_Q = int(cl["class"].isin(["A-Q", "A-QJ"]).sum()); n_J = int(cl["class"].isin(["A-J", "A-QJ"]).sum())
    n_A = int((cl["class"] != "U").sum())
    gained = cl["gained_years"].fillna("").astype(str).str.split().map(len); n_gained_years = int(gained.sum()); n_gained_gages = int((gained > 0).sum())
    lost = cl["lost_years"].fillna("").astype(str).str.split().map(len); n_lost_years = int(lost.sum())
    assert n_U + n_A == n_gages == len(cl), (n_U, n_A, n_gages, len(cl))
    sel = [l.strip() for l in open(os.path.join(R, "manifests", f"selection_{Wd}.txt")) if l.strip()]
    assert sel == summ["gage_id"].tolist(), "the CSV's gage order is not the selection manifest"
    led = pd.read_csv(os.path.join(R, "input", "correction_ledger.csv"), dtype={"gage_id": str, "date": str}, usecols=["gage_id", "date", "status"])
    led = led[led["status"] == "restored"]; n_led_all = int(len(led)); n_led_gages_all = int(led["gage_id"].nunique())
    d = pd.to_datetime(led["date"]); wy = d.dt.year + (d.dt.month >= 10)
    n_high_na = int(summ["flagged_for_high_na"].astype(str).str.upper().isin(["TRUE", "1"]).sum())
    inwin = led[(wy >= W) & (wy <= 2025) & led["gage_id"].isin(set(sel))]
    n_restored_days = int(len(inwin)); n_restored_gages = int(inwin["gage_id"].nunique())
    gates = json.load(open(os.path.join(cen, f"gates_v3_vs_v2_{Wd}.json")))
    glist = gates.get("gates", gates) if isinstance(gates, dict) else gates
    verdict = gates.get("verdict") if isinstance(gates, dict) else None
    if isinstance(glist, dict):
        glist = [{"gate": k, **(v if isinstance(v, dict) else {"status": v})} for k, v in glist.items()]
    n_gate_pass = sum(1 for g in glist if str(g.get("status", g.get("result", ""))).upper() == "PASS") if isinstance(glist, list) else None
    n_gates = len(glist) if isinstance(glist, list) else None
    e7a = json.load(open(os.path.join(e7, "e7_annual.json"))); e7s = json.load(open(os.path.join(e7, "e7_summary.json")))
    m_ot = re.search(r"over-tolerance=([\d,]+)", open(os.path.join(e7, "e7_annual.txt"), encoding="utf-8").read())
    if not m_ot or e7a.get("verdict") != "PASS" or e7s.get("verdict") != "PASS":
        raise SystemExit("E7 verdicts must be PASS and the annual gate text must state its over-tolerance count")
    n_e7_over = int(m_ot.group(1).replace(",", "")); n_e7_pairs = re.search(r"finite pairs=([\d,]+)", open(os.path.join(e7, "e7_annual.txt"), encoding="utf-8").read()).group(1)
    e7_sigs = e7a.get("allow_diff_signature") or []
    for k in ("gated_mismatches", "gate_state_differences", "gate_state_tie_verified", "waived", "tie_sensitive_reported", "annual_residual_reported", "n_gages_compared"):
        if k not in e7s:
            raise SystemExit(f"e7_summary.json lacks {k}")

    # ---- README: the version-2 staged README with the counts and the corrected-input conventions
    t = open(os.path.join(a.staged_v2, f"hisss_readme_{tag}.md"), encoding="utf-8").read()
    old_rows = re.search(r"`value` \(([\d,]+) rows\)", t).group(1)
    t = must_replace(t, f"`value` ({old_rows} rows)", f"`value` ({n_rows:,} rows)", "annual rows")
    t = must_replace(t, f"| `hisss_run_provenance_{tag}.json` | Machine-readable run provenance: input files with sizes, configuration, code revision, software versions |\n",
        f"| `hisss_run_provenance_{tag}.json` | Machine-readable run provenance: input files with sizes and sha256, configuration, code revision (the tagged baseline), software versions, the selection manifest and the climate-frame mode |\n"
        f"| `hisss_run_merge_{tag}.json` | How the eight gage batches were merged (counts per batch, reference order, pinned column schema) |\n"
        f"| `hisss_run_notes_{tag}.md` | Run notes: the exact commands, input hashes, what changed against the previous build and why |\n"
        f"| `hisss_run_log_{tag}.log` | Console logs of the eight batch runs, concatenated |\n"
        f"| `hisss_gage_selection_{tag}.txt` | The frozen gage set of this product (one `gage_id` per line, in the table's row order; sha256 in the provenance JSON) |\n"
        f"| `hisss_would_qualify_not_added_{tag}.csv` | Gages that qualify for this window under the corrected record but are NOT in the product ({n_wq}): listed so the cohort stays the one selected under the February 2026 compile. Columns: `gage_id`, `window`, `baseline_n_valid` / `final_n_valid` (valid water years before / after the correction), `n_gained`, `gained_years`, `final_total_possible`, `final_frac`, `n_restored_days_in_window` |\n"
        f"| `hisss_validation_gates_{tag}.txt`, `hisss_validation_annual_values_{tag}.txt`, `hisss_validation_qualification_audit_{tag}.txt` | Production gates (window, schema, gage set, annual-parquet volume, year counts), the annual-values ↔ summary consistency check, and the qualification audit (emitted = selected ⊆ qualified) |\n"
        f"| `hisss_validation_change_census_{tag}.md` (+ `_gates.json`, `hisss_validation_change_classes_{tag}.csv`) | The change census against the previous build: every changed annual value attributed to a restored day, an inserted climate day or a validity change. The classes file has one row per gage: `class` (U = untouched by both corrections; A-Q = restored days in valid years; A-J = inserted climate days in valid years; A-QJ = both), `restored_years`, `inserted_valid_years`, `n_valid_old` / `n_valid_new`, `gained_years`, `lost_years`, `oracle_gained`. In these file names `v3` is this build and `v2` the previous, undeposited build of October 2026 |\n"
        f"| `hisss_validation_port_check_annual_{tag}.txt`, `_summary_{tag}.txt` | Independent recomputation of the {e7s['n_gages_compared']:,} changed gages in a second implementation (Python port): annual key sets identical, {n_e7_over} of {n_e7_pairs} values beyond 1e-9 (all in the named threshold-tie and float-noise signatures); summary table compared cell for cell at 1e-6 with gate state — 0 gated mismatches; {e7s['waived']} finite-value cells in the named waived columns (`ice_affected_days_total`, FDC90th / FDCmid rank and Pettitt statistics of zero-heavy series) and, separately, {e7s['gate_state_tie_verified']} cells of those columns NaN in one language and finite in the other, accepted only after recomputation reproduced both sides; {e7s['tie_sensitive_reported']:,} tie-bearing rank statistics and {e7s['annual_residual_reported']:,} propagated residuals verified by recomputation |\n",
        "files table")
    t = must_replace(t,
        "- **Resolved known issue.** The 2026 release's `flagged_for_high_na` was computed over the\n  metadata columns only; this version regenerates it with the corrected definition.\n",
        "- **Corrected streamflow input (this version).** The February 2026 compile had set every USGS day\n"
        "  whose qualification code carries an approved remark (`A R` revised, `A <` / `A >` censored,\n"
        f"  `A [4]`) to NA — {n_led_all:,} approved days on {n_led_gages_all} compiled gages. Those days were re-retrieved\n"
        "  from NWIS (26 October 2026) and restored before this product was computed (resource 3 holds\n"
        f"  the compile, the corrected table and the correction ledger). In this window {n_restored_days:,} such\n"
        f"  days on {n_restored_gages:,} gages are finite again, and {n_gained_years:,} water years on {n_gained_gages:,} gages that the\n"
        f"  per-year screening had rejected for the gap they left are valid again ({n_lost_years} years lost). Record\n"
        "  ends 2013–2025 are affected most; trend statistics of those gages can move.\n"
        "- **Climate on Daymet's own calendar (this version).** Each gage's Daymet series now enters the\n"
        "  per-year screening on Daymet's dates inside the gage's record span, so a day absent from the\n"
        "  streamflow record keeps its precipitation instead of being interpolated across\n"
        f"  ({n_J:,} gages have such days in valid years). Precipitation and degree-day totals equal the\n"
        "  daily sums; the Q-to-PPT, snow and temperature families inherit the change.\n"
        f"- **Frozen gage set.** The gage set is the one selected under the February 2026 compile\n"
        f"  (`hisss_gage_selection_{tag}.txt`); {n_U:,} of the {n_gages:,} gages are untouched by both\n"
        f"  corrections and bit-identical to the previous build. {n_wq} gages that would qualify only under\n"
        f"  the corrected record are listed in `hisss_would_qualify_not_added_{tag}.csv` and not added.\n"
        "- **Resolved known issues.** (i) The 2026 release's `flagged_for_high_na` was computed over the\n"
        "  metadata columns only; this version regenerates it with the corrected definition. (ii) The\n"
        "  qualifier mask and (iii) the climate frame above.\n",
        "conventions")
    open(os.path.join(a.out, f"hisss_readme_{tag}.md"), "w", encoding="utf-8").write(t)

    # ---- run notes: commands, hashes, verification
    runs = json.load(open(os.path.join(run, "batch_runs.json"))); b1 = runs.get("1") or next(iter(runs.values()))
    cmd1 = b1.get("cmd"); env1 = {k: v for k, v in (b1.get("env") or {}).items() if str(k).startswith("STREAMFLOW_")}
    fp = b1.get("fingerprint") or prov.get("fingerprint") or {}
    split = json.load(open(os.path.join(R, f"inputs_S2_{WU}", "batches.json")))
    selj = json.load(open(os.path.join(R, "manifests", "selection_manifests.json")))[Wd]
    prev = {"1993": "processedOuts_1993_2025_climate_12oct2026", "1980": "processedOuts_1980_2025_climate_12oct2026"}[str(W)]
    def h(x): return (x or {}).get("sha256", "?")
    notes = f"""# STANDARD PRODUCT WY {W}–2025 @ 60 % — public version 2 (internal build v3, corrected input)

**Run date**: {tj.get('end_time', '?')[:10]} · **Runtime**: {tj.get('total_seconds', 0) / 60:.1f} min of Julia across {mj.get('batches') and len(mj['batches'])} gage batches
(`docs/benchmarks/run_batched_julia.py`, wrapper `docs/benchmarks/run_julia_benchmark_prod_{W}_2025_60pct_climate.jl`)
from a clean checkout of the tagged baseline `qualifier-rerun-baseline-2026-10-26` (commit `{fp.get('git_revision', '?')}`,
working tree clean: {not fp.get('git_working_tree_dirty', True)}; Julia {fp.get('julia_version', '?')}).

**Supersedes the version-2 build of 2026-10-12/13 of the same window** (`{prev}`, never deposited): same window,
qualification rule, configuration, climate input and gage set; two input corrections (README → Conventions):
the USGS qualifier mask of the February 2026 compile is undone and the climate frame is Daymet's own calendar
inside the record span (`STREAMFLOW_CLIMATE_FRAME=outer`).

| | |
|---|---|
| Gages | {n_gages:,} (= the selection manifest, sha256 `{selj.get('sha256', '?')}`) |
| Columns | {n_cols:,} ({n_cols - 32:,} signature + 20 metadata + 12 QA/QC flags; schema pinned to the previous build's header) |
| Annual parquet | {n_rows:,} rows / {n_sigs} signatures |
| Gages untouched by both corrections (bit-identical to the previous build) | {n_U:,} |
| Gages with a restored USGS day inside the window (census class A-Q / A-QJ) / with inserted Daymet days in valid years (A-J / A-QJ) | {n_Q:,} / {n_J:,} |
| Restored days inside the window (product gages) | {n_restored_days:,} on {n_restored_gages:,} gages |
| Water years valid again / lost | {n_gained_years:,} on {n_gained_gages:,} gages / {n_lost_years} |
| Would-qualify gages kept out (listed) | {n_wq} |
| Gages with NA climate summary statistics (fewer than 20 climate-valid years) | {clim_na:,} |
| Gages with raw m³/s discharge (`area_normalized = FALSE`) | {unnorm} |
| `flagged_for_high_na` TRUE (regenerated; > 30 % of the signature-derived fields NA) | {n_high_na:,} |

## Inputs (sha256)

| Input | File | sha256 |
|---|---|---|
| Streamflow (corrected) | `{os.path.basename(prov.get('streamflow', {}).get('path', '?'))}` ({prov.get('streamflow', {}).get('bytes', '?'):,} B) | `{h(prov.get('streamflow'))}` |
| Metadata (year counts recomputed) | `{os.path.basename(prov.get('metadata', {}).get('path', '?'))}` | `{h(prov.get('metadata'))}` |
| Climate | `{os.path.basename(prov.get('climate', {}).get('path', '?'))}` ({prov.get('climate', {}).get('bytes', '?'):,} B) | `{h(prov.get('climate'))}` |
| Selection manifest | `hisss_gage_selection_{tag}.txt` ({selj.get('n_gages', '?')} gages, the previous build's row order) | `{selj.get('sha256', '?')}` |
| Pinned schema | the previous build's summary CSV (2026-10-12/13, not deposited; its header is the column contract) | `{h(prov.get('qa_schema'))}` |
| Configuration | `config/signatures_config.json` | `{h(fp.get('config'))}` |

## Commands (as recorded by the batch runner)

```
# split the corrected streamflow and the Daymet input into {split.get('batches')} batches over the selection manifest ({split.get('assignment')})
run_batched_julia.py split --streamflow {os.path.basename(split.get('streamflow', '?'))} --batches {split.get('batches')} --gage-list selection_{Wd}.txt --climate daymet2025={os.path.basename(prov.get('climate', {}).get('path', '?'))}
# run (batch 1 shown; the environment is identical for every batch except the input files)
{' '.join(cmd1) if isinstance(cmd1, list) else cmd1}
{json.dumps(env1, indent=0)}
# merge in the previous build's gage order with the pinned schema
run_batched_julia.py merge --out-dir <run folder> --prefix {p} --reference <previous build's summary CSV> --qa-schema <same>
```

## Verification (files beside this one)

- Production gates 6 / 6 PASS (`hisss_validation_gates_{tag}.txt`); annual values consistent with the summary
  (`hisss_validation_annual_values_{tag}.txt`); qualification audit: emitted = selected ⊆ qualified, would-qualify
  kept out {n_wq} (`hisss_validation_qualification_audit_{tag}.txt`).
- Change census against the previous build — verdict **{verdict}** ({n_gate_pass} / {n_gates} gates PASS): every one of the
  changed annual values is attributed to a restored day, an inserted climate day or a validity change; the {n_U:,}
  untouched gages are identical at both layers; precipitation and degree-day totals equal the daily sums
  (`hisss_validation_change_census_{tag}.md`).
- Independent recomputation of the {e7s['n_gages_compared']:,} changed gages in the Python port
  (`hisss_validation_port_check_annual_{tag}.txt`, `_summary_{tag}.txt`): annual key sets identical ({e7a['n_keys_ref']:,} keys),
  {n_e7_over} of {n_e7_pairs} finite values beyond 1e-9, every one in the named signatures {', '.join(e7_sigs)}
  (±1-day threshold ties and float noise ≤ 1e-10 relative); summary table compared cell for cell at 1e-6 with gate
  state: 0 gated mismatches; {e7s['waived']} finite-value cells in the named waived columns (`ice_affected_days_total`; the
  FDC90th / FDCmid rank and Pettitt statistics of zero-heavy series, where the two languages' near-zero tails differ) and,
  separately, {e7s['gate_state_tie_verified']} cells of those columns that are NaN in one language and finite in the other,
  accepted only after the statistic recomputed from each side's own annual series reproduced that side; {e7s['tie_sensitive_reported']:,} tie-bearing rank
  statistics and {e7s['annual_residual_reported']:,} propagated annual residuals verified by recomputation.
- The code, the plan and the full verification record are in the HISSS code repository (CHANGELOG October 2026;
  `docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md` §11).
"""
    open(os.path.join(a.out, f"hisss_run_notes_{tag}.md"), "w", encoding="utf-8").write(notes)
    with open(os.path.join(a.out, f"hisss_rename_manifest_{tag}.csv"), "w") as f:
        f.write("original_name,deposit_name\n")
        for o, d_ in manifest:
            f.write(f"{o},{d_}\n")
    print(f"staged {a.out}: {n_gages:,} gages x {n_cols:,} cols; {n_rows:,} annual rows / {n_sigs} signatures; U {n_U:,}, A-Q {n_Q:,}, A-J {n_J:,}; "
          f"restored days {n_restored_days:,} on {n_restored_gages:,} gages; years +{n_gained_years:,}/−{n_lost_years}; would-qualify {n_wq}; census {verdict} {n_gate_pass}/{n_gates}")


# ----------------------------------------------------------------------------------------------- inputs (R3)
def stage_inputs(a):
    R = a.rerun_dir; os.makedirs(a.out, exist_ok=True); manifest = []
    inp = os.path.join(R, "input")
    # unchanged version-2 files
    for name in ("hisss_daymet_basin_daily.parquet", "hisss_daymet_basin_flags.csv", "hisss_daymet_basin_daily_provenance.json", "hisss_hydat_interference.csv"):
        manifest.append(copy(os.path.join(a.staged_v2, name), os.path.join(a.out, name)))
    # the immutable February compile (verified against the control provenance) and the corrected table
    cp = json.load(open(os.path.join(inp, "control_provenance.json")))
    base_sha = sha(a.base_parquet)
    if base_sha != cp["parquet"]["sha256"]:
        raise SystemExit(f"base parquet sha256 {base_sha} differs from the control provenance {cp['parquet']['sha256']}")
    manifest.append(copy(a.base_parquet, os.path.join(a.out, "hisss_streamflow_daily_as_compiled_feb2026.parquet")))
    manifest.append(copy(os.path.join(inp, "combined_streamflow_data_09feb2026_qualifier_patched.parquet"), os.path.join(a.out, "hisss_streamflow_daily.parquet")))
    manifest.append(copy(os.path.join(a.staged_v2, "hisss_gage_metadata.csv"), os.path.join(a.out, "hisss_gage_metadata_as_compiled_feb2026.csv")))
    manifest.append(copy(os.path.join(inp, "combined_watershed_metadata_09feb2026_qualifier_patched.csv"), os.path.join(a.out, "hisss_gage_metadata.csv")))
    manifest.append(copy(os.path.join(inp, "correction_ledger.csv"), os.path.join(a.out, "hisss_streamflow_qualifier_corrections.csv")))
    manifest.append(copy(os.path.join(inp, "metadata_changes.csv"), os.path.join(a.out, "hisss_gage_metadata_changes.csv")))
    for src, dst in [("control_rule_vs_vintage.csv", "hisss_qualifier_control_report.csv"), ("control_rule_vs_vintage_years.csv", "hisss_qualifier_control_report_years.csv"),
                     ("control_provenance.json", "hisss_qualifier_control_provenance.json"), ("input_audit.md", "hisss_streamflow_qualifier_audit.md"),
                     ("input_audit.json", "hisss_streamflow_qualifier_audit.json")]:
        manifest.append(copy(os.path.join(inp, src), os.path.join(a.out, dst)))
    manifest.append(copy(os.path.join(inp, "oracle", "acceptance_oracle.csv"), os.path.join(a.out, "hisss_qualifier_acceptance_oracle.csv")))
    manifest.append(copy(os.path.join(inp, "oracle", "would_qualify_outside_products.csv"), os.path.join(a.out, "hisss_would_qualify_not_added.csv")))
    for Wd in ("1993-2025", "1980-2025"):
        manifest.append(copy(os.path.join(R, "manifests", f"selection_{Wd}.txt"), os.path.join(a.out, f"hisss_gage_selection_wy{Wd.replace('-', '_')}.txt")))
    manifest.append(copy(os.path.join(R, "manifests", "selection_manifests.json"), os.path.join(a.out, "hisss_gage_selection_manifests.json")))
    manifest.append(copy(os.path.join(HERE, "reconstruct_qualifier_patched_streamflow.py"), os.path.join(a.out, "hisss_streamflow_reconstruction_recipe.py")))
    # raw NWIS responses (gzip JSON per site + fetch_meta.json) as one tarball
    tgz = os.path.join(a.out, "hisss_nwis_qualifier_responses.tar.gz"); cache = os.path.join(inp, "nwis_cache")
    with tarfile.open(tgz, "w:gz") as tf:
        for n in sorted(os.listdir(cache)):
            tf.add(os.path.join(cache, n), arcname=os.path.join("nwis_cache", n))
    manifest.append(("nwis_cache/", os.path.basename(tgz)))
    # numbers
    led = pd.read_csv(os.path.join(inp, "correction_ledger.csv"), dtype={"gage_id": str}, usecols=["gage_id", "status"])
    n_led = int((led["status"] == "restored").sum()); n_led_gages = int(led.loc[led["status"] == "restored", "gage_id"].nunique())
    n_sites = len([n for n in os.listdir(cache) if n.endswith(".json.gz")])
    mc = pd.read_csv(os.path.join(inp, "metadata_changes.csv"), dtype={"gage_id": str}); n_meta = int(mc["gage_id"].nunique()) if "gage_id" in mc.columns else len(mc)
    ctl = pd.read_csv(os.path.join(inp, "control_rule_vs_vintage.csv"), dtype={"gage_id": str})
    n_cmp = int(ctl["n_compared"].sum()); n_eq = int(ctl["n_equal_1e9"].sum())
    n_meta_changed = int(mc["values_changed"].astype(str).str.lower().isin(["true", "1"]).sum())
    import pyarrow.parquet as pq
    n_na_new = int(pq.read_table(os.path.join(inp, "combined_streamflow_data_09feb2026_qualifier_patched.parquet"), columns=["Q"]).column("Q").null_count)
    n_na_old = int(pq.read_table(a.base_parquet, columns=["Q"]).column("Q").null_count)
    assert n_na_old - n_na_new == n_led, (n_na_old, n_na_new, n_led)
    led_full = pd.read_csv(os.path.join(inp, "correction_ledger.csv"), dtype=str, usecols=["gage_id", "date", "status_today", "status"])
    led_full = led_full[led_full["status"] == "restored"]; prov_today = led_full[led_full["status_today"] == "P"]
    n_prov_today = int(len(prov_today)); prov_gages = sorted(prov_today["gage_id"].unique()); prov_span = (prov_today["date"].min(), prov_today["date"].max()) if n_prov_today else None
    meta_all = pd.read_csv(os.path.join(inp, "combined_watershed_metadata_09feb2026_qualifier_patched.csv"), dtype=str, usecols=["gage_id", "processing_status"])
    compiled = set(meta_all.loc[meta_all["processing_status"] == "success", "gage_id"].str.lstrip("0"))
    daymet_sites = set(pd.read_parquet(os.path.join(a.staged_v2, "hisss_daymet_basin_daily.parquet"), columns=["site_id"])["site_id"].astype(str).unique())
    no_forcing = sorted(g for g in compiled if g not in {x.lstrip("0") for x in daymet_sites}); n_compiled = len(compiled); n_with = n_compiled - len(no_forcing)
    in_products = {g: [w for w in ("1993-2025", "1980-2025") if g in {x.lstrip("0") for x in open(os.path.join(R, "manifests", f"selection_{w}.txt")).read().split()}] for g in no_forcing}
    m_prov = re.search(r"provisional-remark days ([\d,]+)", open(os.path.join(inp, "input_audit.md"), encoding="utf-8").read())
    prov_days = f"{int(m_prov.group(1).replace(',', '')):,} days" if m_prov else "the audit's count"
    # README
    t = open(os.path.join(a.staged_v2, "hisss_readme_inputs.md"), encoding="utf-8").read()
    t = must_replace(t,
        "| `hisss_streamflow_daily.parquet` | Daily discharge: 111,624,189 rows, 8,014 gages, 1 Oct 1979 – 1 Oct 2025. Columns: `gage_id`, `Date`, `Q`, derived calendar fields, and the agency quality `flag` |\n"
        "| `hisss_gage_metadata.csv` | One row per candidate gage (16,994): coordinates, published drainage area, retrieval extent, `processing_status`, and the `area_normalized` flag |\n",
        "| `hisss_streamflow_daily.parquet` | Daily discharge, **corrected** (the table the signature products were computed from): 111,624,189 rows, 8,014 gages, 1 Oct 1979 – 1 Oct 2025. Columns: `gage_id`, `Date`, `Q`, derived calendar fields, and the agency quality `flag`. Differs from the February 2026 compile below in exactly the "
        f"{n_led:,} `Q` cells of the correction ledger |\n"
        "| `hisss_streamflow_daily_as_compiled_feb2026.parquet` | The February 2026 compile, unchanged (same rows, order and schema): the USGS days carrying an approved remark code are NA in it (see Provenance) |\n"
        f"| `hisss_streamflow_qualifier_corrections.csv` | The correction ledger: one row per restored day ({n_led:,} rows, {n_led_gages} gages) — the February flag and value, the value and qualifiers served by NWIS on 26 Oct 2026, the conversion and the restored `Q` |\n"
        f"| `hisss_nwis_qualifier_responses.tar.gz` | The raw NWIS daily-value responses the corrections were taken from ({n_sites} sites, gzip JSON, request window 1979-10-01 – 2025-10-01) with the fetch record |\n"
        f"| `hisss_qualifier_control_report.csv`, `_years.csv`, `_provenance.json` | Three-dataset control: for every corrected gage, the February values compared with the current NWIS record ({n_eq:,} of {n_cmp:,} equal at 1e-9; the rest are agency revisions since February, reported and NOT applied), per gage and per year |\n"
        "| `hisss_qualifier_acceptance_oracle.csv` | Per gage and window: the valid water years before and after the correction under the products' screening rules (the restored years, gage for gage) |\n"
        "| `hisss_streamflow_qualifier_audit.md` (+ `.json`) | The audit of the correction (fetch completeness, cell-level diff of the two tables, conversion identity, accounting, metadata) |\n"
        "| `hisss_gage_selection_wy1993_2025.txt`, `hisss_gage_selection_wy1980_2025.txt`, `hisss_gage_selection_manifests.json` | The frozen gage sets of the two signature products (one id per line, in the tables' row order; hashes in the JSON) |\n"
        "| `hisss_would_qualify_not_added.csv` | Gages that qualify under the corrected record but are not in the products (the cohorts are frozen at the February selection) |\n"
        "| `hisss_streamflow_reconstruction_recipe.py` | Executable recipe: applies the ledger to the February compile and verifies the result against the corrected table row group by row group |\n"
        f"| `hisss_gage_metadata.csv` | One row per candidate gage (16,994): coordinates, published drainage area, retrieval extent, `processing_status`, and the `area_normalized` flag; the derived year counts (`num_years`, `start_year`, `end_year`, `num_days`) recomputed from the corrected table for the {n_meta} corrected gages — {n_meta_changed} of them change (`hisss_gage_metadata_changes.csv` lists the old and new values; `hisss_gage_metadata_as_compiled_feb2026.csv` is the February table) |\n",
        "R3 files table")
    t = must_replace(t,
        "- **Streamflow** was retrieved with the USGS dataRetrieval R package (US) and the\n  Water Survey of Canada tidyhydat R package (Canada), retaining agency quality codes\n  (including ice-affected and estimated values).",
        "- **Streamflow** was retrieved with the USGS dataRetrieval R package (US) and the\n  Water Survey of Canada tidyhydat R package (Canada), retaining agency quality codes\n  (including ice-affected and estimated values).\n"
        "- **Qualifier correction (this version).** The February 2026 compile kept a US value only when its\n"
        "  qualification code was exactly `A`, `A e`, `P` or `P e`, so every approved day carrying a remark\n"
        f"  code — `A R` (records revised), `A <` / `A >` (censored), `A [4]` — had been set to NA: {n_led:,} days on\n"
        f"  {n_led_gages} gages. The corrected table restores them from a re-retrieval of those sites on 26 October 2026\n"
        "  (status A or P, finite value, same conversion as the compile; censored values enter as the\n"
        f"  reported bound; {len(led_full) - n_prov_today:,} of the restored values are approved today and {n_prov_today} provisional — "
        f"{', '.join(prov_gages)}, {prov_span[0]} to {prov_span[1]} — kept because the compile's rule accepts a provisional value); the `flag`\n"
        "  column keeps the February code. Nothing else changed: the control\n"
        "  report compares every February value of the corrected gages with the current record, the\n"
        "  audit diffs the two tables cell by cell, and the recipe reproduces the corrected table from the\n"
        f"  compile and the ledger. Provisional-remark days (`P Ice`, `P Dis`, …; {prov_days}) stay NA: their\n"
        "  February state cannot be reconstructed from today's record.",
        "R3 provenance")
    m_na = re.search(r"with NA discharge \(~?[\d,]+ rows\)", t)
    if not m_na:
        raise SystemExit("R3 README: the NA-row sentence was not found")
    t = t.replace(m_na.group(0), f"with NA discharge ({n_na_new:,} rows in the corrected table, {n_na_old:,} in the February compile)")
    t = must_replace(t,
        "- **This resource contains raw-but-harmonized inputs**: per-water-year quality\n  screening (gap and missing-value rules, interpolation of gaps ≤ 3 days) was applied\n  downstream, inside the signature computation — not to these files.",
        "- **This resource contains raw-but-harmonized inputs**: per-water-year quality\n  screening (gap and missing-value rules, interpolation of gaps ≤ 3 days) was applied\n  downstream, inside the signature computation — not to these files. The signature products\n  (public version 2) were computed from the corrected streamflow table with each gage's Daymet\n  series joined on Daymet's own calendar inside the gage's record span.",
        "R3 screening bullet")
    nf_txt = ", ".join(f"{g} (in {'neither product' if not in_products[g] else 'product WY ' + ' and WY '.join(in_products[g])})" for g in no_forcing) if no_forcing else "none"
    t = must_replace(t, "one gage without forcing carries streamflow-only signatures in the signature products.",
                     f"compiled gage{'s' if len(no_forcing) != 1 else ''} without forcing: {nf_txt}.", "R3 no-forcing sentence")
    open(os.path.join(a.out, "hisss_readme_inputs.md"), "w", encoding="utf-8").write(t)
    # abstract
    ab = open(os.path.join(a.staged_v2, "hisss_abstract_resource3.txt"), encoding="utf-8").read()
    ab = must_replace(ab, "with agency quality codes retained;",
        f"with agency quality codes retained (deposited twice: the February 2026 compile, in which {n_led:,} approved USGS days carrying a remark code had been set to NA, and the corrected table the products were computed from, with those days restored from a re-retrieval and a correction ledger, raw responses, control report and reconstruction recipe beside it);",
        "abstract")
    ab = must_replace(ab, "extends it to every compiled basin and to 2025", f"extends it to {n_with:,} of the {n_compiled:,} compiled basins and to 2025", "abstract coverage")
    open(os.path.join(a.out, "hisss_abstract_resource3.txt"), "w", encoding="utf-8").write(ab)
    # dictionary: amend the flag / Q / year-count rows, add the new tables
    dd = pd.read_csv(os.path.join(a.staged_v2, "hisss_input_data_dictionary.csv"), dtype=str, keep_default_na=False)
    def amend(table, col, extra):
        m = (dd["table"] == table) & (dd["column"] == col); assert m.sum() == 1, (table, col)
        dd.loc[m, "definition"] = dd.loc[m, "definition"] + " " + extra
    amend("hisss_streamflow_daily.parquet", "flag", f"The flag is the February 2026 retrieval's code; in hisss_streamflow_daily.parquet the {n_led:,} US days whose approved remark code (R, <, >, [4]) had masked the value carry the restored value under that same flag (ledger: hisss_streamflow_qualifier_corrections.csv). hisss_streamflow_daily_as_compiled_feb2026.parquet has the same columns with those cells NA.")
    amend("hisss_streamflow_daily.parquet", "Q", "Corrected table: finite on the restored days of the correction ledger; otherwise identical to the February compile.")
    m_id = (dd["table"] == "hisss_gage_metadata.csv") & (dd["column"] == "gage_id"); assert m_id.sum() == 1
    dd.loc[m_id, "definition"] = ("Agency station identifier AS RETRIEVED: in this table, US site numbers are stored ZERO-STRIPPED (e.g. 1011000 for USGS site 01011000), "
                                  "while the streamflow and Daymet tables store the zero-padded form, so a direct join fails for most US gages. To join: strip leading zeros "
                                  "on BOTH sides (never re-pad — US ids are not all 8 digits, and padding to 8 leaves compiled gages unmatched); Canadian identifiers are "
                                  "unaffected. Read as a character string")
    for c in ("num_years", "start_year", "end_year", "num_days"):
        amend("hisss_gage_metadata.csv", c, f"Recomputed from the corrected streamflow table for the {n_meta} corrected gages (hisss_gage_metadata_changes.csv); hisss_gage_metadata_as_compiled_feb2026.csv keeps the February values.")
    new_rows = [
        ("hisss_streamflow_qualifier_corrections.csv", "gage_id", "-", "USGS site number, zero-padded (joins hisss_streamflow_daily.parquet)"),
        ("hisss_streamflow_qualifier_corrections.csv", "date", "date", "Calendar date of the restored day"),
        ("hisss_streamflow_qualifier_corrections.csv", "old_flag", "-", "The February 2026 qualification code (A R, A <, A >, A [4])"),
        ("hisss_streamflow_qualifier_corrections.csv", "old_Q", "mm/day", "The February value (always NA: the masked cell)"),
        ("hisss_streamflow_qualifier_corrections.csv", "value_today_cfs", "ft3/s", "Daily mean discharge served by NWIS at the re-retrieval (value_today_raw is the served string)"),
        ("hisss_streamflow_qualifier_corrections.csv", "qualifiers_today", "-", "Qualification code served at the re-retrieval; status_today is its status letter (A approved / P provisional)"),
        ("hisss_streamflow_qualifier_corrections.csv", "Q_new", "mm/day", "The restored value = value_today_cfs × conversion (the compile's cfs → mm/day factor from basin_area_km2); for the 73 unnormalized gages no US day is involved"),
        ("hisss_streamflow_qualifier_corrections.csv", "conversion", "(mm/day)/(ft3/s)", "The compile's cfs → mm/day factor, re-derived from the metadata: 0.0283168 m3/s per ft3/s × 86400 s/day × 1000 mm/m / (basin_area_km2 × 1e6 m2/km2) = 2.44657 / basin_area_km2"),
        ("hisss_streamflow_qualifier_corrections.csv", "method_id", "-", "NWIS method (time-series) identifier the value was read from (the unqualified block where a site publishes several)"),
        ("hisss_streamflow_qualifier_corrections.csv", "retrieved", "UTC timestamp", "Time of the NWIS request"),
        ("hisss_streamflow_qualifier_corrections.csv", "status", "-", "restored (every row: the day was eligible and restored)"),
        ("hisss_qualifier_control_report.csv", "gage_id", "-", "USGS site number; one row per corrected gage"),
        ("hisss_qualifier_control_report.csv", "n_compared / n_equal_1e9 / n_differ", "count", "February values with a value served today: how many equal value_today × conversion at 1e-9, how many differ (agency revisions since February, reported and not applied); max_rel_diff is the largest relative difference"),
        ("hisss_qualifier_control_report.csv", "n_restored_by_rule / n_restored_by_rule_splice_codes", "count", "Days the acceptance rule restores on this gage / of which carry one of the four corrected codes in February (the difference are plain A / P days today that were masked under another code)"),
        ("hisss_qualifier_control_report_years.csv", "gage_id, water_year, …", "-", "The same comparison per gage and water year"),
        ("hisss_qualifier_acceptance_oracle.csv", "gage_id, window", "-", "Gage and product window (1993-2025 or 1980-2025); in_product says whether the gage is in that product"),
        ("hisss_qualifier_acceptance_oracle.csv", "baseline_n_valid / final_n_valid / gained_years / lost_years", "count, years", "Valid water years under the products' screening before / after the correction, and the years gained / lost; baseline_qualifies / final_qualifies apply the 20-year and 60 % rules"),
        ("hisss_qualifier_acceptance_oracle.csv", "n_restored_days_in_window / n_source_unavailable_in_window", "count", "Restored days inside the window / masked days no current record could restore (0 everywhere)"),
        ("hisss_gage_metadata_changes.csv", "gage_id", "-", "One row per corrected gage (987): the metadata row recomputed from the corrected table"),
        ("hisss_gage_metadata_changes.csv", "old_num_years, old_start_year, old_end_year, old_num_days, old_error_message", "count, year, year, count, -", "The February values of the derived year-count columns"),
        ("hisss_gage_metadata_changes.csv", "new_num_years, new_start_year, new_end_year, new_num_days", "count, year, year, count", "The values recomputed from the corrected table with the compile's rule"),
        ("hisss_gage_metadata_changes.csv", "values_changed", "boolean", f"TRUE when at least one of the four values changed ({n_meta_changed} of {n_meta} rows; the others restored days inside years that were already counted)"),
        ("hisss_gage_selection_wy1993_2025.txt / hisss_gage_selection_wy1980_2025.txt", "gage_id", "-", "One gage id per line: the frozen gage set of the product, in the product table's row order"),
        ("hisss_would_qualify_not_added.csv", "(oracle columns)", "-", "The acceptance-oracle rows of gages that qualify under the corrected record but are not in the products"),
    ]
    dd = pd.concat([dd, pd.DataFrame(new_rows, columns=["table", "column", "units", "definition"])], ignore_index=True)
    dd.to_csv(os.path.join(a.out, "hisss_input_data_dictionary.csv"), index=False)
    with open(os.path.join(a.out, "hisss_rename_manifest_inputs.csv"), "w") as f:
        f.write("original_name,deposit_name\n")
        for o, d_ in manifest:
            f.write(f"{o},{d_}\n")
    print(f"staged {a.out}: ledger {n_led:,} rows / {n_led_gages} gages; cache {n_sites} sites; metadata rows recomputed for {n_meta} gages; control {n_eq:,} of {n_cmp:,} equal")


# ----------------------------------------------------------------------------------------------- check
REQUIRED = {
    "signatures": ["hisss_signatures_{tag}.csv", "hisss_signatures_annual_{tag}.parquet", "hisss_run_provenance_{tag}.json", "hisss_run_merge_{tag}.json",
                   "hisss_run_notes_{tag}.md", "hisss_readme_{tag}.md", "hisss_data_dictionary.csv", "hisss_signature_categories.csv",
                   "hisss_signature_explorer_{tag}.html", "hisss_signature_explorer_{tag}_annual", "hisss_gage_selection_{tag}.txt",
                   "hisss_would_qualify_not_added_{tag}.csv", "hisss_validation_gates_{tag}.txt", "hisss_validation_annual_values_{tag}.txt",
                   "hisss_validation_qualification_audit_{tag}.txt", "hisss_validation_change_census_{tag}.md", "hisss_validation_change_census_gates_{tag}.json",
                   "hisss_validation_change_classes_{tag}.csv", "hisss_validation_port_check_annual_{tag}.txt", "hisss_validation_port_check_summary_{tag}.txt",
                   "hisss_rename_manifest_{tag}.csv"],
    "inputs": ["hisss_streamflow_daily.parquet", "hisss_streamflow_daily_as_compiled_feb2026.parquet", "hisss_streamflow_qualifier_corrections.csv",
               "hisss_nwis_qualifier_responses.tar.gz", "hisss_qualifier_control_report.csv", "hisss_qualifier_control_report_years.csv",
               "hisss_qualifier_control_provenance.json", "hisss_qualifier_acceptance_oracle.csv", "hisss_streamflow_qualifier_audit.md",
               "hisss_streamflow_qualifier_audit.json", "hisss_gage_selection_wy1993_2025.txt", "hisss_gage_selection_wy1980_2025.txt",
               "hisss_gage_selection_manifests.json", "hisss_would_qualify_not_added.csv", "hisss_streamflow_reconstruction_recipe.py",
               "hisss_gage_metadata.csv", "hisss_gage_metadata_as_compiled_feb2026.csv", "hisss_gage_metadata_changes.csv",
               "hisss_daymet_basin_daily.parquet", "hisss_daymet_basin_flags.csv", "hisss_daymet_basin_daily_provenance.json",
               "hisss_hydat_interference.csv", "hisss_input_data_dictionary.csv", "hisss_readme_inputs.md", "hisss_abstract_resource3.txt",
               "hisss_rename_manifest_inputs.csv"],
}


def check(a):
    tag = f"wy{a.window}_2025" if a.window else ""
    req = [r.format(tag=tag) for r in REQUIRED[a.kind]]
    missing = [r for r in req if not os.path.exists(os.path.join(a.out, r))]
    rows = []; bad = []
    for n in sorted(os.listdir(a.out)):
        p = os.path.join(a.out, n)
        if os.path.isdir(p):
            rows.append((n + "/", sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(p) for f in fs), f"{len(os.listdir(p))} files")); continue
        if n.startswith("hisss_staging_manifest_"):
            continue
        rows.append((n, os.path.getsize(p), sha(p)))
        if n.endswith(".parquet"):
            with open(p, "rb") as fh:
                fh.seek(-4, 2); ok = fh.read(4) == b"PAR1"
            if not ok:
                bad.append(f"{n}: no PAR1 footer")
    if a.kind == "signatures":
        html = os.path.join(a.out, f"hisss_signature_explorer_{tag}.html")
        if os.path.exists(html):
            t = open(html, "rb").read(); want = f'"dir":"hisss_signature_explorer_{tag}_annual"'.encode()
            if t.count(want) != 1:
                bad.append("explorer: embedded sidecar folder name not rewritten")
            side = os.path.join(a.out, f"hisss_signature_explorer_{tag}_annual")
            if os.path.isdir(side) and not any(f.endswith(".js") for f in os.listdir(side)):
                bad.append("explorer: sidecar folder holds no .js files")
    present = set(os.listdir(a.out))
    for doc in [n for n in present if n.startswith(("hisss_readme", "hisss_run_notes"))]:
        txt = open(os.path.join(a.out, doc), encoding="utf-8").read()
        for ref in sorted(set(re.findall(r"`(hisss_[A-Za-z0-9_./-]+)`", txt))):
            base_ = ref.rstrip("/")
            if base_ not in present and not base_.startswith("hisss_staging_manifest"):
                bad.append(f"{doc} names {ref}, which is not in the folder")
    out = os.path.join(a.out, f"hisss_staging_manifest_{a.kind}{'_' + tag if tag else ''}.csv")
    pd.DataFrame(rows, columns=["file", "bytes", "sha256"]).to_csv(out, index=False)
    print(f"{a.out}: {len(rows)} entries; required present {len(req) - len(missing)} / {len(req)}; missing {missing}; problems {bad}; manifest {os.path.basename(out)}")
    sys.exit(0 if not missing and not bad else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("signatures"); s.add_argument("--run-dir", required=True); s.add_argument("--prefix", required=True)
    s.add_argument("--window", type=int, required=True, choices=[1993, 1980]); s.add_argument("--rerun-dir", required=True)
    s.add_argument("--staged-v2", required=True); s.add_argument("--out", required=True); s.set_defaults(fn=stage_signatures)
    i = sp.add_parser("inputs"); i.add_argument("--rerun-dir", required=True); i.add_argument("--staged-v2", required=True)
    i.add_argument("--base-parquet", required=True); i.add_argument("--out", required=True); i.set_defaults(fn=stage_inputs)
    c = sp.add_parser("check"); c.add_argument("--out", required=True); c.add_argument("--kind", required=True, choices=["signatures", "inputs"])
    c.add_argument("--window", type=int, default=None, choices=[1993, 1980]); c.set_defaults(fn=check)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
