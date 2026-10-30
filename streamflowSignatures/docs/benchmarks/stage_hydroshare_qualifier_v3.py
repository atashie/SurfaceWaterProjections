#!/usr/bin/env python
"""Stage the HydroShare resource folders for the version-3 products after the qualifier-mask + climate-frame rerun (plan
docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md, Phase H; 2026-10-28). PUBLIC FRAMING (user decision 2026-10-29):
the deposit is the first and only public version — no version label, no "previous build", "corrected", "regenerated"
or "this version" wording in any user-facing text; validation files that only make sense against an internal build
(change census, classes, port check) are not deposited. The README comes from `hydroshare_readme_signatures.tmpl.md`.

    signatures --run-dir RUN --prefix PREFIX --window 1993|1980 --rerun-dir R --staged-v2 <v2 staged R1|R2> --out OUT
    inputs     --rerun-dir R --staged-v2 <v2 staged R3> --out OUT [--deposit minimal|full] [--base-parquet <February compile>]
    check      --out OUT --kind signatures|inputs [--window W]

Local staging only (nothing is uploaded). `signatures` copies and renames the run's artefacts to the deposit names of
2026 (summary CSV, annual parquet, provenance JSON, explorer + `_annual/` sidecar with the embedded folder name
rewritten, the validation record), adds the frozen selection manifest, the would-qualify list, the change census,
the qualification audit and the port check, and rewrites the version-2 README's counts and conventions. `inputs`
stages resource 3: by default (`--deposit minimal`, user decision 2026-10-29) the six files of the August 2026 upload,
replaced in place — the streamflow table the products were computed from, the gage metadata, the Daymet table, the
HYDAT interference table, the dictionary and the README (abstract and manifests stay local) — the R3 deposit of 2026-10-29 used
`--companions` as well (user decision: + the correction ledger, the Daymet flags, the gage-selection lists and the would-qualify list); `--deposit full` adds
the immutable February compile beside the table, the correction ledger, the raw NWIS responses, the three-dataset
control report, the acceptance oracle, the selection manifests and the executable reconstruction recipe. `check` verifies that every required
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
        if os.path.exists(dst):
            os.remove(dst)
        r = subprocess.run(["cp", "-p", "-c", src, dst], capture_output=True)   # APFS clone (seconds for a 5 GB file)
        if r.returncode != 0:
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


def public_explorer_strings(html_path, run_csv_name, deposit_csv_name, title_drop):
    """Public strings inside the built explorer: the subtitle's source file name becomes the deposit name and the
    internal version fragment is removed from the three title occurrences (<title>, <h1>, meta.title)."""
    with open(html_path, "rb") as fh:
        t = fh.read()
    old = f'"source":"{run_csv_name}"'.encode(); new = f'"source":"{deposit_csv_name}"'.encode()
    if t.count(old) != 1:
        raise SystemExit(f"{html_path}: expected exactly one embedded source name {old!r}, found {t.count(old)}")
    t = t.replace(old, new)
    if title_drop:
        frag = title_drop.encode(); n = t.count(frag)
        if n != 3:
            raise SystemExit(f"{html_path}: expected the title fragment {title_drop!r} 3 times (<title>, <h1>, meta.title), found {n}")
        t = t.replace(frag, b"")
    with open(html_path, "wb") as fh:
        fh.write(t)


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
    public_explorer_strings(expl_dst, f"{p}_signatures.csv", f"hisss_signatures_{tag}.csv", a.explorer_title_drop)
    if a.deposit == "minimal":
        # the deposit replaces the seven files of the August 2026 upload, so the sidecar travels as a zip that unpacks
        # to the folder the HTML expects (entries sorted, mtimes = the run artefacts', so the zip is reproducible)
        import zipfile
        zpath = os.path.join(a.out, f"{side_dst}.zip")
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
            for n in sorted(os.listdir(os.path.join(a.out, side_dst))):
                if not n.startswith("."):
                    zf.write(os.path.join(a.out, side_dst, n), arcname=f"{side_dst}/{n}")
        manifest.append((side_dst + "/", os.path.basename(zpath)))
    # validation record (self-contained production checks only; the change census, the class file and the port check
    # compare against an internal build and are NOT deposited — user decision 2026-10-29)
    cen = os.path.join(R, f"census_{WU}"); e7 = os.path.join(R, f"e7_{WU}")
    for src, dst in [(os.path.join(run, "validation_gates.txt"), f"hisss_validation_gates_{tag}.txt"),
                     (os.path.join(run, "validate_annual_values_floor20.log"), f"hisss_validation_annual_values_{tag}.txt"),
                     (os.path.join(run, "audit_qualification.log"), f"hisss_validation_qualification_audit_{tag}.txt"),
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
    # dictionary: the climate-extended dictionary with three definitions reworded (no version history in the text) and
    # the rows sorted into the summary table's column order; categories unchanged
    dd = pd.read_csv(os.path.join(a.staged_v2, "hisss_data_dictionary.csv"), dtype=str, keep_default_na=False)
    m_hn = dd["column_name"] == "flagged_for_high_na"; assert m_hn.sum() == 1
    dd.loc[m_hn, "definition"] = ("Automated QA flag: TRUE when more than 30% of this gage's signature-derived fields (every column except the "
                                  "metadata and the other flags) are NA; one shared definition in the Julia, Python and R implementations.")
    m_ice = dd["column_name"] == "ice_affected_days_total"; assert m_ice.sum() == 1
    n_ice_nonzero = int((pd.to_numeric(pd.read_csv(os.path.join(a.out, f"hisss_signatures_{tag}.csv"), usecols=["ice_affected_days_total"])["ice_affected_days_total"], errors="coerce").fillna(0) != 0).sum())
    dd.loc[m_ice, "definition"] = (f"Diagnostic intended as the total of ice-affected missing days across the record (from agency quality codes). KNOWN ISSUE: "
                                   f"the implementation reports 0 for every gage ({n_ice_nonzero} non-zero values in this table), so a zero does NOT "
                                   "establish the absence of ice-affected days — use the daily `flag` column of the input-data resource (code `Ice`, HYDAT `B`) instead.")
    m_id = dd["column_name"] == "gage_id"; assert m_id.sum() == 1
    gid = dd.loc[m_id, "definition"].iloc[0]
    if gid.count("; see README)") != 1:
        raise SystemExit("dictionary gage_id definition: the README pointer was not found exactly once")
    if gid.count("some HISSS tables store long USGS ids without their leading zero") != 1:
        raise SystemExit("dictionary gage_id definition: the zero-stripped clause was not found exactly once")
    gid = gid.replace("some HISSS tables store long USGS ids without their leading zero",
                      "the gage-metadata table stores every USGS id without its leading zeros and the geometry / land-cover tables some long ones")
    dd.loc[m_id, "definition"] = gid.replace("; see README)", "; see README → Reading the data)")
    # categories: the repo's canonical sheet (docs/signature_categories.csv; user decision 2026-10-29: canonical categories
    # throughout the deposit) replaces the exploratory grouping in the dictionary's category column and in the categories file
    sheet = pd.read_csv(os.path.join(HERE, "..", "signature_categories.csv"), dtype=str, keep_default_na=False)
    cat_of = dict(zip(sheet["signature"], sheet["category"]))
    m_sig = dd["signature"] != ""
    if set(dd.loc[m_sig, "signature"]) != set(sheet["signature"]):
        raise SystemExit("dictionary signatures and the canonical category sheet do not match one to one")
    dd.loc[m_sig, "category"] = dd.loc[m_sig, "signature"].map(cat_of)
    n_scalars = int((dd["statistic"] == "per-gage scalar").sum())
    if n_scalars != int((sheet["kind"] == "scalar").sum()):
        raise SystemExit("scalar count differs between the dictionary and the canonical category sheet")
    header = pd.read_csv(os.path.join(a.out, f"hisss_signatures_{tag}.csv"), nrows=0).columns.tolist()
    if set(header) != set(dd["column_name"]) or dd["column_name"].duplicated().any():
        raise SystemExit("dictionary rows do not match the summary table's columns one to one")
    dd = dd.set_index("column_name").loc[header].reset_index()[dd.columns]
    dd.to_csv(os.path.join(a.out, "hisss_data_dictionary.csv"), index=False); manifest.append(("hisss_data_dictionary.csv", "hisss_data_dictionary.csv"))
    order = {sg: i for i, sg in reversed(list(enumerate(dd.loc[dd["signature"] != "", "signature"])))}   # first appearance in column order
    cats = sheet[["signature", "kind", "category", "function_family"]].copy()
    cats = cats.assign(_o=cats["signature"].map(order)).sort_values("_o", kind="stable").drop(columns="_o")
    cats.to_csv(os.path.join(a.out, "hisss_signature_categories.csv"), index=False)
    manifest.append(("docs/signature_categories.csv", "hisss_signature_categories.csv"))
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

    # ---- README from the template (first and only public version: every number from the files)
    from string import Template
    tmpl = Template(open(os.path.join(HERE, "hydroshare_readme_signatures.tmpl.md"), encoding="utf-8").read())
    n_qann_unnorm = int((summ["area_normalized"].astype(str).str.upper().isin(["FALSE", "0"])
                         & summ["flagged_for_qann_range"].astype(str).str.upper().isin(["TRUE", "1"])).sum())
    n_batches = len(mj["batches"]) if mj.get("batches") else len(runs)
    words = {4: "four", 6: "six", 8: "eight", 10: "ten", 12: "twelve", 16: "sixteen"}
    clim_note = ((f" {clim_na:,} gages have fewer than 20 climate-valid water years"
                  + (" (WY 1980 has no climate value for any gage: the input begins 1980-01-01)" if W == 1980 else "")
                  + ", so their climate summary statistics are NA while their climate annual values are in the parquet.")
                 if clim_na else " Every gage carries climate summary statistics.")
    nb_word = words.get(len(mj["batches"]) if mj.get("batches") else len(runs), "")
    core_rows = [
        f"| `hisss_signatures_{tag}.csv` | Signature summary table: one row per gage, {n_cols:,} columns (gage metadata, {n_sigs} signatures × 16 statistics, 21 per-gage scalars, 12 QA flags) |",
        f"| `hisss_signatures_annual_{tag}.parquet` | The annual values behind every statistic, long format: `gage_id`, `signature`, `water_year`, `value` ({n_rows:,} rows) |",
        "| `hisss_data_dictionary.csv` | Definition, units, category, and statistic for every column of the summary table, in the table's column order |",
        f"| `hisss_signature_categories.csv` | Category and computing function family of each of the {n_sigs} annual signatures and the {n_scalars} per-gage scalars (`kind` = annual / scalar); the nine categories are the eight of the HISSS data paper plus Climate; `ice_affected_days_total` carries no category (a preprocessing diagnostic) |",
    ]
    if a.deposit == "minimal":
        files_rows = core_rows + [
            f"| `hisss_signature_explorer_{tag}.html` + `hisss_signature_explorer_{tag}_annual.zip` | Interactive map explorer: unzip `hisss_signature_explorer_{tag}_annual.zip` beside the HTML (it unpacks to the folder `hisss_signature_explorer_{tag}_annual/`, which the per-gage time-series panels load from) and open the HTML in any browser — it needs an internet connection for the map library and the basemap tiles |",
        ]
        gage_set_files = (f"The gage list (`hisss_gage_selection_{tag}.txt`) and the list of those {n_wq} gages "
                          f"(`hisss_would_qualify_not_added.csv`) are in the HISSS input-data resource.")
    else:
        files_rows = core_rows + [
            f"| `hisss_signature_explorer_{tag}.html` + `hisss_signature_explorer_{tag}_annual/` | Interactive map explorer (open the HTML in any browser — it needs an internet connection for the map library and the basemap tiles; keep the `_annual/` folder beside it — the per-gage time-series panels load from it) |",
            f"| `hisss_run_provenance_{tag}.json` | Machine-readable run provenance: input files with sizes and sha256, configuration, code revision, software versions, the selection manifest and the climate-frame mode |",
            f"| `hisss_run_merge_{tag}.json` | How the {nb_word} gage batches were merged (counts per batch, reference order, pinned column schema) |",
            f"| `hisss_run_notes_{tag}.md` | Run notes: the exact commands, input hashes and the verification record |",
            f"| `hisss_run_log_{tag}.log` | Console logs of the {nb_word} batch runs, concatenated |",
            f"| `hisss_gage_selection_{tag}.txt` | The gage set of this product (one `gage_id` per line, in the table's row order; sha256 in the provenance JSON) |",
            f"| `hisss_would_qualify_not_added_{tag}.csv` | Gages that meet the qualification rule under the corrected streamflow table but are NOT in the product ({n_wq}): the gage set was fixed under the February 2026 compile (see Conventions → Gage set). Columns: `gage_id`, `window`, `baseline_n_valid` / `final_n_valid` (valid water years under the February compile / under the corrected table), `n_gained`, `gained_years`, `final_total_possible`, `final_frac`, `n_restored_days_in_window` |",
            f"| `hisss_validation_gates_{tag}.txt`, `hisss_validation_annual_values_{tag}.txt`, `hisss_validation_qualification_audit_{tag}.txt` | Production gates (window, schema, gage set, annual-parquet volume, year counts), the annual-values ↔ summary consistency check, and the qualification audit (emitted = selected ⊆ qualified) |",
            f"| `hisss_rename_manifest_{tag}.csv` | Maps each file's name in the processing run to its name in this resource |",
            f"| `hisss_staging_manifest_signatures_{tag}.csv` | Byte size and sha256 of every file in this resource (the sidecar folder as one aggregate row) |",
        ]
        gage_set_files = (f"The gage list is `hisss_gage_selection_{tag}.txt` and the {n_wq} gages are listed in "
                          f"`hisss_would_qualify_not_added_{tag}.csv`, both in this resource.")
    ctx = dict(files_table="\n".join(files_rows), gage_set_files=gage_set_files, W=W, tag=tag, companion=1980 if W == 1993 else 1993, n_gages=f"{n_gages:,}", n_cols=f"{n_cols:,}", n_rows=f"{n_rows:,}",
               n_sigs=n_sigs, n_sig_cols=f"{n_cols - 32:,}", n_batches_word=words.get(n_batches, str(n_batches)), n_wq=n_wq,
               n_led_all=f"{n_led_all:,}", n_led_gages_all=n_led_gages_all, n_restored_days=f"{n_restored_days:,}",
               n_restored_gages=f"{n_restored_gages:,}", n_gained_years=f"{n_gained_years:,}", n_gained_gages=f"{n_gained_gages:,}",
               n_J=f"{n_J:,}", unnorm=unnorm, n_qann_unnorm=n_qann_unnorm, n_high_na=f"{n_high_na:,}", clim_note=clim_note,
               wy1980_note=" WY 1980 has none (the input begins 1980-01-01);" if W == 1980 else "",
               pettitt_note="; WY 2025 is not included." if W == 1980 else f" — for this product that is effectively WY {W}–2024; WY 2025 is not included.")
    readme_text = tmpl.substitute(ctx)   # render first: a template error must not truncate the staged file
    open(os.path.join(a.out, f"hisss_readme_{tag}.md"), "w", encoding="utf-8").write(readme_text)

    # ---- run notes: commands, hashes, verification (no previous-build text)
    runs = json.load(open(os.path.join(run, "batch_runs.json"))); b1 = runs.get("1") or next(iter(runs.values()))
    cmd1 = b1.get("cmd"); env1 = {k: v for k, v in (b1.get("env") or {}).items() if str(k).startswith("STREAMFLOW_")}
    fp = b1.get("fingerprint") or prov.get("fingerprint") or {}
    split = json.load(open(os.path.join(R, f"inputs_S2_{WU}", "batches.json")))
    selj = json.load(open(os.path.join(R, "manifests", "selection_manifests.json")))[Wd]
    def h(x): return (x or {}).get("sha256", "?")
    jv = str(fp.get("julia_version", "?")).replace("julia version ", "")
    notes = f"""# STANDARD PRODUCT WY {W}–2025 @ 60 % — run record

**Run date**: {tj.get('end_time', '?')[:10]} · **Runtime**: {tj.get('total_seconds', 0) / 60:.1f} min of Julia across {n_batches} gage batches
(`docs/benchmarks/run_batched_julia.py`, wrapper `docs/benchmarks/run_julia_benchmark_prod_{W}_2025_60pct_climate.jl`)
from a clean checkout of the HISSS code at commit `{fp.get('git_revision', '?')}` (working tree clean:
{not fp.get('git_working_tree_dirty', True)}; Julia {jv}).

**Inputs**: the corrected daily streamflow table and the Daymet 1980–2025 table of the HISSS input-data resource
(README → Conventions → Streamflow input, Climate forcing); the climate frame is Daymet's own calendar inside each
gage's record span (`STREAMFLOW_CLIMATE_FRAME=outer`); a gage qualifies with 20 valid water years that are 60 % of
the window's years.

| | |
|---|---|
| Gages | {n_gages:,} (= the selection manifest, sha256 `{selj.get('sha256', '?')}`) |
| Columns | {n_cols:,} ({n_cols - 32:,} signature + 20 metadata + 12 QA/QC flags; schema pinned to the HISSS column contract) |
| Annual parquet | {n_rows:,} rows / {n_sigs} signatures |
| Gages with a restored USGS day inside the window / with Daymet days absent from the streamflow record in valid years | {n_Q:,} / {n_J:,} |
| Restored days inside the window (product gages) | {n_restored_days:,} on {n_restored_gages:,} gages |
| Water years valid only because of the restored days | {n_gained_years:,} on {n_gained_gages:,} gages |
| Would-qualify gages not included (listed) | {n_wq} |
| Gages with NA climate summary statistics (fewer than 20 climate-valid years) | {clim_na:,} |
| Gages with raw m³/s discharge (`area_normalized = FALSE`) | {unnorm} |
| `flagged_for_high_na` TRUE (> 30 % of the signature-derived fields NA) | {n_high_na:,} |

## Inputs (sha256)

| Input | File | sha256 |
|---|---|---|
| Streamflow (corrected table) | `{os.path.basename(prov.get('streamflow', {}).get('path', '?'))}` ({prov.get('streamflow', {}).get('bytes', '?'):,} B) | `{h(prov.get('streamflow'))}` |
| Metadata (year counts from the corrected table) | `{os.path.basename(prov.get('metadata', {}).get('path', '?'))}` | `{h(prov.get('metadata'))}` |
| Climate | `{os.path.basename(prov.get('climate', {}).get('path', '?'))}` ({prov.get('climate', {}).get('bytes', '?'):,} B) | `{h(prov.get('climate'))}` |
| Selection manifest | `hisss_gage_selection_{tag}.txt` ({selj.get('n_gages', '?')} gages, the table's row order) | `{selj.get('sha256', '?')}` |
| Column contract | the reference summary header the merge pinned the schema to | `{h(prov.get('qa_schema'))}` |
| Configuration | `config/signatures_config.json` | `{h(fp.get('config'))}` |

## Commands (as recorded by the batch runner)

```
# split the streamflow and the Daymet input into {split.get('batches')} batches over the selection manifest ({split.get('assignment')})
run_batched_julia.py split --streamflow {os.path.basename(split.get('streamflow', '?'))} --batches {split.get('batches')} --out-dir <inputs folder> --gage-list selection_{Wd}.txt --climate daymet2025={os.path.basename(prov.get('climate', {}).get('path', '?'))}
# run (batch 1 shown; the environment is identical for every batch except the input files)
{' '.join(cmd1) if isinstance(cmd1, list) else cmd1}
{json.dumps(env1, indent=0)}
# merge in the selection manifest's gage order with the pinned schema
run_batched_julia.py merge --out-dir <run folder> --prefix {p} --reference <reference summary CSV: gage order> --qa-schema <same: column contract>
```

## Verification (files beside this one)

- Production gates 6 / 6 PASS (`hisss_validation_gates_{tag}.txt`); annual values consistent with the summary
  (`hisss_validation_annual_values_{tag}.txt`); qualification audit: emitted = selected ⊆ qualified, would-qualify
  not included {n_wq} (`hisss_validation_qualification_audit_{tag}.txt`).
- The code, the configuration, the input audit and the cross-implementation validation (the canonical Julia
  implementation against its Python and R ports) are in the HISSS code repository, https://github.com/CZ-Sync/HISSS .
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
    if not a.base_parquet:
        raise SystemExit("--base-parquet (the February 2026 compile) is required for --deposit full")
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
        "- **Qualifier correction.** The February 2026 compile kept a US value only when its\n"
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
        "- **This resource contains raw-but-harmonized inputs**: per-water-year quality\n  screening (gap and missing-value rules, interpolation of gaps ≤ 3 days) was applied\n  downstream, inside the signature computation — not to these files. The signature products\n  were computed from the corrected streamflow table with each gage's Daymet\n  series joined on Daymet's own calendar inside the gage's record span.",
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
    ab = must_replace(ab, "(version 2, aggregated with", "(aggregated with", "abstract version label")
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


# ----------------------------------------------------------------------------------------------- inputs (R3), six files
HYDAT_RELEASE = "14 October 2025"   # the quarterly HYDAT release current at the 7 February 2026 retrieval (docs/DATA_SOURCES.md, 2026-09-10)
DAYMET_REPRODUCED_BASINS = "5,965"  # basins of the co-authors' original aggregation that the Daymet build reproduces (EO_data_processing/README_DAYMET.md)
DAYMET_LARGE_BASINS = 53            # basins above 100,000 km² added by that build (same README)
COMPILE_MIN_DAYS, COMPILE_MIN_Q = 30, 0.0001   # the compile's year rule (splice_nwis_qualifiers.compile_year_stats; R/helperFunctions.R)
REQUEST_WINDOW = ("1 Oct 1979", "1 Oct 2025")  # config.R DEFAULT_START_DATE / DEFAULT_END_DATE of the ingestion (retrieved 7 February 2026)
# reproduction of the co-authors' 1980–2023 series on the 5,965 shared basins: exhaustive daily maxima measured 2026-10-29 (blind review B),
# rounded up; README_DAYMET.md quotes 0.015 / 0.003 / 0.26 / 0.26 / 0.09 and ±0.012 %
DAYMET_TOL = "daily differences within 0.016 mm for precipitation, 0.003 °C for temperature, 0.27 mm for snow water equivalent, 0.26 Pa for vapour pressure and 0.093 W/m² for radiation; annual precipitation totals within ±0.012 %"


def _rg_nulls(pf, col):
    i = pf.schema_arrow.names.index(col); n = 0
    for g in range(pf.metadata.num_row_groups):
        s = pf.metadata.row_group(g).column(i).statistics
        if s is None or not s.has_null_count:
            raise SystemExit(f"{col}: row-group null statistics missing (row group {g})")
        n += s.null_count
    return n


def _rg_minmax(pf, col):
    i = pf.schema_arrow.names.index(col); mn = mx = None
    for g in range(pf.metadata.num_row_groups):
        s = pf.metadata.row_group(g).column(i).statistics
        if s is None or not s.has_min_max:
            raise SystemExit(f"{col}: row-group min/max statistics missing (row group {g})")
        mn = s.min if mn is None or s.min < mn else mn; mx = s.max if mx is None or s.max > mx else mx
    return mn, mx


def stage_inputs_minimal(a):
    """Resource 3 as the six files of the August 2026 upload, replaced in place (user decision 2026-10-29: a deposit
    replaces only the files already in the resource; first and only public version — no version wording): the
    streamflow table the products were computed from, the gage metadata, the Daymet table, the HYDAT interference
    table, the dictionary and the README. The abstract, the rename manifest and the staging manifest stay local."""
    import glob
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    R = a.rerun_dir; os.makedirs(a.out, exist_ok=True); manifest = []
    inp = os.path.join(R, "input")
    q_src = os.path.join(inp, "combined_streamflow_data_09feb2026_qualifier_patched.parquet")
    m_src = os.path.join(inp, "combined_watershed_metadata_09feb2026_qualifier_patched.csv")
    c_src = os.path.join(a.staged_v2, "hisss_daymet_basin_daily.parquet")
    h_src = os.path.join(a.staged_v2, "hisss_hydat_interference.csv")
    for s in (q_src, m_src, c_src, h_src):
        if not os.path.exists(s):
            raise SystemExit(f"missing {s}")
    # the three data inputs must be the files every product recorded (sha256 in each product's merged timing JSON)
    want = {"streamflow": sha(q_src), "metadata": sha(m_src), "climate": sha(c_src)}
    tjs = sorted(glob.glob(os.path.join(R, "processedOuts_*", "*_timing.json")))
    if not tjs:
        raise SystemExit(f"no product timing JSON under {R}/processedOuts_*")
    for tj in tjs:
        p = json.load(open(tj))["provenance"]
        for k, v in want.items():
            if p[k]["sha256"] != v:
                raise SystemExit(f"{k}: the staged file ({v[:16]}…) is not the input recorded by {os.path.basename(tj)} ({p[k]['sha256'][:16]}…)")
    manifest.append(copy(q_src, os.path.join(a.out, "hisss_streamflow_daily.parquet")))
    manifest.append(copy(m_src, os.path.join(a.out, "hisss_gage_metadata.csv")))
    manifest.append(copy(c_src, os.path.join(a.out, "hisss_daymet_basin_daily.parquet")))
    manifest.append(copy(h_src, os.path.join(a.out, "hisss_hydat_interference.csv")))
    key = lambda s: s.lstrip("0") if s.isdigit() else s.upper()
    # --- streamflow table: rows, gages, span, NA rows, last Canadian day (row-group statistics + one pass over gage_id / Date)
    pf = pq.ParquetFile(q_src); n_rows_q = pf.metadata.num_rows
    n_na_q = _rg_nulls(pf, "Q"); q_d0, q_d1 = _rg_minmax(pf, "Date")
    gq = set(); can_last = None; n_neg = 0; n_us_null_flag = 0; n_us_null_flag_finite = 0; neg_gages = set()
    for g in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(g, columns=["gage_id", "Date", "Q", "flag"])
        gq.update(t.column("gage_id").unique().to_pylist())
        is_ca = pc.match_substring_regex(t.column("gage_id"), "[A-Za-z]")
        ca = t.filter(is_ca)
        if ca.num_rows:
            d = pc.max(ca.column("Date")).as_py(); can_last = d if can_last is None or d > can_last else can_last
        neg = pc.fill_null(pc.less(t.column("Q"), 0), False)
        n_neg += pc.sum(neg).as_py() or 0
        neg_gages.update(t.filter(neg).column("gage_id").unique().to_pylist())
        us_null = pc.and_(pc.invert(is_ca), pc.is_null(t.column("flag")))
        n_us_null_flag += pc.sum(us_null).as_py() or 0
        n_us_null_flag_finite += pc.sum(pc.and_(us_null, pc.is_valid(t.column("Q")))).as_py() or 0
    n_gages_q = len(gq)
    if a.base_parquet:   # optional accounting against the February compile: its NA rows = this table's + the ledger
        cp_ = json.load(open(os.path.join(inp, "control_provenance.json")))
        if sha(a.base_parquet) != cp_["parquet"]["sha256"]:
            raise SystemExit("--base-parquet is not the February compile of the control provenance")
        n_na_base = _rg_nulls(pq.ParquetFile(a.base_parquet), "Q")
    # --- Daymet table
    pfc = pq.ParquetFile(c_src); n_rows_c = pfc.metadata.num_rows; c_d0, c_d1 = _rg_minmax(pfc, "Date")
    for v in ("prcp", "tmin", "tmax", "swe", "vp", "srad"):
        if _rg_nulls(pfc, v):
            raise SystemExit(f"Daymet table: {v} has missing values")
    sites = set()
    for g in range(pfc.metadata.num_row_groups):
        sites.update(pfc.read_row_group(g, columns=["site_id"]).column("site_id").unique().to_pylist())
    n_sites = len(sites)
    # --- metadata
    meta = pd.read_csv(m_src, dtype=str, keep_default_na=False); n_meta_rows = len(meta)
    st = meta["processing_status"].value_counts().to_dict()
    compiled = meta[meta["processing_status"] == "success"]; n_compiled = len(compiled)
    an_true = int((meta["area_normalized"].str.upper() == "TRUE").sum()); an_false = int((meta["area_normalized"].str.upper() == "FALSE").sum())
    n_unnorm = int((compiled["area_normalized"].str.upper() == "FALSE").sum())
    comp_keys = set(compiled["gage_id"].map(key)); q_keys = {key(x) for x in gq}
    if comp_keys != q_keys:
        raise SystemExit(f"the compiled gages of the metadata ({len(comp_keys)}) and the gages of the streamflow table ({len(q_keys)}) differ")
    site_keys = {key(x) for x in sites}; pad = {key(x): x for x in gq}
    no_forcing = sorted(comp_keys - site_keys); n_with = n_compiled - len(no_forcing); n_extra = len(site_keys - comp_keys)
    sel = {Wd: {key(x) for x in open(os.path.join(R, "manifests", f"selection_{Wd}.txt")).read().split()} for Wd in ("1993-2025", "1980-2025")}
    meta_ids = set(compiled["gage_id"]); n_direct = len(gq & meta_ids)   # compiled gages whose metadata id equals the parquet id as stored
    n_pad = sum(1 for x in meta_ids if (x.zfill(8) if x.isdigit() and len(x) <= 8 else x) not in gq)   # compiled gages a pad-to-8 join misses
    unnorm_keys = set(compiled.loc[compiled["area_normalized"].str.upper() == "FALSE", "gage_id"].map(key))
    n_unnorm_p = {w: len(unnorm_keys & sel[w]) for w in sel}
    in_products = {g: [w for w in ("1993-2025", "1980-2025") if g in sel[w]] for g in no_forcing}
    nf_txt = ", ".join(f"`{pad[g]}` (in {'neither product' if not in_products[g] else 'the WY ' + ' and the WY '.join(in_products[g]) + ' product'})"
                       for g in no_forcing) if no_forcing else "none"
    # --- the qualifier facts (ledger + masked-day census of the rerun) and the HYDAT table
    led = pd.read_csv(os.path.join(inp, "correction_ledger.csv"), dtype=str, usecols=["gage_id", "status", "status_today"])
    led = led[led["status"] == "restored"]; n_led = len(led); n_led_gages = int(led["gage_id"].nunique()); n_prov_today = int((led["status_today"] == "P").sum())
    mk = pd.read_csv(os.path.join(inp, "masked_days_by_gage_code.csv"), dtype=str)
    pr = mk[mk["category"] == "provisional_remark"]; n_prem = int(pr["n_days"].astype(int).sum()); n_prem_gages = int(pr["gage_id"].nunique())
    if a.base_parquet and n_na_base - n_na_q != n_led:
        raise SystemExit(f"NA accounting: base {n_na_base:,} − table {n_na_q:,} ≠ ledger {n_led:,}")
    hy = pd.read_csv(h_src, dtype=str, keep_default_na=False); n_hydat = len(hy); n_reg_empty = int((hy["REGULATED"] == "").sum())
    fl = pd.read_csv(os.path.join(a.staged_v2, "hisss_daymet_basin_flags.csv"), dtype=str, keep_default_na=False)
    tb = lambda c: fl[c].str.upper() == "TRUE"
    n_flag = len(fl); n_lowconf = int(tb("low_confidence").sum()); n_hb = int(tb("hydrobasins_fallback").sum()); n_am = int(tb("area_mismatch").sum())
    n_nar = int(tb("no_area_reference").sum()); n_lps = int(tb("low_pixel_support").sum()); n_both = int((tb("low_confidence") & tb("low_pixel_support")).sum())
    if n_lowconf + n_lps - n_both != n_flag:
        raise SystemExit("Daymet flags: the two groups do not account for every flagged basin")
    both_txt = f"; {n_both} basin{'s are' if n_both != 1 else ' is'} in both groups" if n_both else ""
    overlap_txt = "; the reasons overlap" if n_hb + n_am + n_nar > n_lowconf else ""
    n_p1 = int(tb("in_p1").sum()); n_p2 = int(tb("in_p2").sum()); n_either = int((tb("in_p1") | tb("in_p2")).sum())
    fmt_d = lambda d: f"{d.day} {d.strftime('%b')} {d.year}"
    # --- optional companions (--companions): the small files the signature readmes point readers to
    comp = bool(getattr(a, "companions", False))
    if comp:
        manifest.append(copy(os.path.join(inp, "correction_ledger.csv"), os.path.join(a.out, "hisss_streamflow_qualifier_corrections.csv")))
        manifest.append(copy(os.path.join(a.staged_v2, "hisss_daymet_basin_flags.csv"), os.path.join(a.out, "hisss_daymet_basin_flags.csv")))
        manifest.append(copy(os.path.join(inp, "oracle", "would_qualify_outside_products.csv"), os.path.join(a.out, "hisss_would_qualify_not_added.csv")))
        for Wd in ("1993-2025", "1980-2025"):
            manifest.append(copy(os.path.join(R, "manifests", f"selection_{Wd}.txt"), os.path.join(a.out, f"hisss_gage_selection_wy{Wd.replace('-', '_')}.txt")))
        wq = pd.read_csv(os.path.join(inp, "oracle", "would_qualify_outside_products.csv"), dtype=str)
        n_wq = {w: int((wq["window"] == w).sum()) for w in ("1993-2025", "1980-2025")}
        n_sel = {w: len(sel[w]) for w in sel}
    comp_rows = ("" if not comp else
        f"| `hisss_streamflow_qualifier_corrections.csv` | One row per USGS day retrieved on 26 October 2026 ({n_led:,} rows, {n_led_gages} gages): the February code and NA value, the value and qualifiers served by NWIS at the retrieval, the conversion and the stored `Q` |\n"
        f"| `hisss_daymet_basin_flags.csv` | One row per flagged basin ({n_flag}): low-confidence outlines (HydroBASINS fallback, polygon-vs-published area mismatch, no area reference) and basins covered by few Daymet cells; the reason, both areas, cell counts and product membership. Flagged, not removed — filter on it if your analysis is sensitive to basin delineation |\n"
        f"| `hisss_gage_selection_wy1993_2025.txt`, `hisss_gage_selection_wy1980_2025.txt` | The gage sets of the two signature products ({n_sel['1993-2025']:,} and {n_sel['1980-2025']:,} ids, one per line, in the product tables' row order) |\n"
        f"| `hisss_would_qualify_not_added.csv` | Gages that meet the products' qualification rule under this streamflow table but are not in the products ({n_wq['1993-2025']} for WY 1993–2025, {n_wq['1980-2025']} for WY 1980–2025): the gage sets were fixed before the remark-code days were retrieved. Columns: `gage_id`, `window`, the valid water years before / after those days (`baseline_n_valid` / `final_n_valid`, `baseline_valid_years` / `final_valid_years`), `gained_years`, `baseline_qualifies` / `final_qualifies`, `n_restored_days_in_window` |\n")
    dict_scope = "tables above" if comp else "four tables above"
    kept_txt = f"none was removed ({n_either} of them are in a signature product: {n_p1} in WY 1993–2025, {n_p2} in WY 1980–2025)"
    flags_txt = (f"{kept_txt}; `hisss_daymet_basin_flags.csv` lists them with the reason, both areas, the cell counts and product membership." if comp else
                 f"{kept_txt}; the list with the reason for each basin is available from the authors.")
    ledger_txt = (f" `hisss_streamflow_qualifier_corrections.csv` lists those days with the values and qualifiers served at the retrieval." if comp else "")
    gagesets_bullet = ("" if not comp else
        f"- **Gage sets.** A gage enters a signature product with ≥ 20 valid water years that are ≥ 60 % of the\n"
        f"  window's years. The two gage lists are the products' row order ({n_sel['1993-2025']:,} and {n_sel['1980-2025']:,} gages); the\n"
        f"  {n_wq['1993-2025']} + {n_wq['1980-2025']} gages of `hisss_would_qualify_not_added.csv` meet the rule under this table only because of\n"
        f"  the days retrieved on 26 October 2026 and are not in the products, whose gage sets were fixed before that retrieval.\n")
    # --- README
    readme = f"""# HISSS Harmonized Input Data

Part of the **Hydrologic Information Signatures and Summary Statistics (HISSS)**
dataset. This resource contains the harmonized daily inputs from which the two HISSS
streamflow-signature products (WY 1993–2025 and WY 1980–2025; DOIs to be assigned at
publication) were computed: daily discharge with the agencies' quality codes for
{n_gages_q:,} gaged watersheds across the United States and Canada, the gage-level
processing metadata for all {n_meta_rows:,} candidate gages, basin-averaged daily Daymet
meteorological and snow forcing for {n_sites:,} basins ({n_with:,} of the {n_compiled:,} compiled gages and
{n_extra} basins whose gages were not compiled; calendar years 1980–2025), and Canadian station
interference metadata.
Depositing the inputs alongside the signature products makes the full HISSS workflow
reproducible end-to-end with the archived code.

## Files

| File | Description |
|---|---|
| `hisss_streamflow_daily.parquet` | Daily discharge: {n_rows_q:,} rows, {n_gages_q:,} gages, {fmt_d(q_d0)} – {fmt_d(q_d1)}. Columns: `gage_id`, `Date`, `Q`, derived calendar fields, and the agency quality `flag` |
| `hisss_gage_metadata.csv` | One row per candidate gage ({n_meta_rows:,}): coordinates, published drainage area, retrieval extent, `processing_status`, and the `area_normalized` flag |
| `hisss_daymet_basin_daily.parquet` | Basin-averaged daily Daymet V4 R1 forcing: {n_rows_c:,} rows, {n_sites:,} basins, calendar years {c_d0.year}–{c_d1.year}, no missing days. Variables: `prcp`, `tmin`, `tmax`, `swe`, `vp`, `srad` |
| `hisss_hydat_interference.csv` | RHBN (reference network) and regulation status for {n_hydat:,} Canadian stations, exported from the HYDAT database (STATIONS and STN_REGULATION tables) via tidyhydat |
{comp_rows}| `hisss_input_data_dictionary.csv` | Definition and units for every column of the {dict_scope} |

## Reading the data

```python
import pandas as pd
q    = pd.read_parquet("hisss_streamflow_daily.parquet")
clim = pd.read_parquet("hisss_daymet_basin_daily.parquet")
met  = pd.read_csv("hisss_gage_metadata.csv", dtype={{"gage_id": str}})
# metadata stores US ids zero-stripped; join on a stripped key on BOTH sides (do not pad):
key = lambda s: s.lstrip("0") if s.isdigit() else s.upper()
met["key"] = met.gage_id.map(key); q["key"] = q.gage_id.map(key)
joined = q.merge(met, on="key", how="left", suffixes=("", "_meta"))
```

```r
q    <- arrow::read_parquet("hisss_streamflow_daily.parquet")
clim <- arrow::read_parquet("hisss_daymet_basin_daily.parquet")
met  <- readr::read_csv("hisss_gage_metadata.csv",
                        col_types = readr::cols(gage_id = readr::col_character()))
key <- function(x) ifelse(grepl("^[0-9]+$", x), sub("^0+", "", x), toupper(x))
met$key <- key(met$gage_id); q$key <- key(q$gage_id)
joined <- dplyr::left_join(q, met, by = "key")
```

Both parquet tables are large ({n_rows_q / 1e6:,.1f} M and {n_rows_c / 1e6:,.1f} M rows; 0.6 GB and 5.1 GB on disk): read
the columns, gages or years you need (pyarrow row groups, DuckDB) rather than whole tables.

**⚠ Join-key convention.** The two parquet tables store `gage_id` / `site_id`
**zero-padded** (USGS site 01011000 appears as `01011000`) and join to each other
directly on (`gage_id` = `site_id`, `Date`). The metadata CSV, however, stores every US id
**as retrieved — zero-stripped** (`1011000`): a naive join to the parquets silently
matches only {n_direct:,} of the {n_compiled:,} compiled gages ({100 * n_direct / n_compiled:.0f} %). **Strip leading zeros on both sides** as shown above;
Canadian identifiers are unaffected. Do NOT pad the metadata ids to 8 digits instead:
9- and 10-digit USGS site numbers that begin with a zero are stored here as 8–9 digits,
so padding to 8 leaves them wrong and silently mismatches {n_pad} of the {n_compiled:,} compiled gages.
In all cases read identifiers as character strings, never numeric.

## Provenance and conventions

- **Streamflow** was retrieved on 7 February 2026 with the USGS dataRetrieval R package
  (US) and the Water Survey of Canada tidyhydat R package (Canada), retaining agency
  quality codes (including ice-affected and estimated values). A US daily value is kept
  when its qualification code begins with the status letter `A` (approved) or `P`
  (provisional) and the value is finite; the remark codes that may follow the letter
  (`e` estimated, `R` revised, `<` / `>` censored at a reporting limit, `[4]`, `Ice`,
  `Bkw`, …) never mask a value and are kept verbatim in `flag`. The {n_led:,} approved days
  carrying a remark code (`A R`, `A <`, `A >`, `A [4]`; {n_led_gages} gages) were retrieved
  from NWIS on 26 October 2026, the February retrieval having left them without a value;
  {n_prov_today} of them are provisional in today's record.{ledger_txt} The {n_prem:,} provisional days
  carrying a remark code (`P Ice`, `P Dis`, `P Rat`, …; {n_prem_gages} gages) have NA
  discharge in this table. Missing days appear either as rows with NA discharge
  ({n_na_q:,} rows) or as absent rows — never zero-filled; zero flow is a valid
  observation, and negative values are kept as served ({n_neg:,} rows on {len(neg_gages)} gages).
  Requests covered {REQUEST_WINDOW[0]} – {REQUEST_WINDOW[1]}, so a
  single day of water year 2026 ({fmt_d(q_d1)}) is present.
- **HYDAT.** The Canadian records come from the HYDAT release current at the retrieval —
  the release of {HYDAT_RELEASE}, inferred from the retrieval date and the tidyhydat version
  used rather than read from the database — and end on {fmt_d(can_last)}, so no Canadian
  gage has a water year 2025 in either signature product. The interference table (RHBN
  membership and regulation status) was exported from the same release; `REGULATED` is
  empty for the {n_reg_empty} stations without a recorded regulation status.
- **Units**: `Q` is mm/day, converted using agency-published drainage areas. For the
  **{n_unnorm} compiled gages** with `area_normalized = FALSE` (Canadian stations with no
  published drainage area — mostly canals, diversions, dam outflows, and channel
  splits), `Q` remains in native m³/s; {n_unnorm_p['1993-2025']} of them are in the WY 1993–2025 product and
  {n_unnorm_p['1980-2025']} in the WY 1980–2025 product. Filter on `area_normalized` before comparing
  flow magnitudes across gages. (The metadata table flags {an_false:,} candidate stations
  FALSE in total, but {an_false - n_unnorm:,} of those were not compiled — no data returned, or fewer
  than 20 years — and have no rows in the streamflow table.)
- **Gage metadata.** `num_years`, `start_year`, `end_year` and `num_days` describe the
  compiled record: the calendar years with at least {COMPILE_MIN_DAYS} days of non-missing discharge
  above {COMPILE_MIN_Q} (mm/day, or m³/s for the gages that are not area-normalized), and the
  number of daily rows in those years; they agree with the streamflow table as deposited.
- **This resource contains raw-but-harmonized inputs**: per-water-year quality
  screening (gap and missing-value rules, interpolation of gaps ≤ 3 days) was applied
  downstream, inside the signature computation — not to these files. Consult the data
  paper's methods for those rules. The signature products were computed from these
  tables, with each gage's Daymet series joined on Daymet's own calendar inside the
  gage's record span (a day absent from the discharge record keeps its forcing).
- **Daymet** (built September 2026) is the project's area-weighted aggregation of the
  1 km Daymet V4 R1 grids to {n_sites:,} gaged basins: fractional cell-coverage (true-area)
  weights over the full-resolution basin polygons, six variables, calendar years
  1980–2025, no missing days. On the {DAYMET_REPRODUCED_BASINS} basins covered by the co-authors'
  original aggregation (gdptools, 1980–2023) it reproduces that series ({DAYMET_TOL}); it adds the
  basins that lacked one, the {DAYMET_LARGE_BASINS} basins above 100,000 km², and 2024–2025. {n_flag}
  basins have outlines of lower confidence or few Daymet cells ({n_lowconf} low-confidence
  outlines — {n_hb} HydroBASINS fallbacks, {n_am} polygon-vs-published area mismatches,
  {n_nar} without an area reference{overlap_txt} — and {n_lps} basins covered by few cells{both_txt}); {flags_txt}
  Two properties to know: (i) **Daymet publishes a 365-day
  calendar** — leap years keep February 29 and omit December 31, so every year has
  exactly 365 rows per basin; (ii) `swe` is a model product with documented biases,
  strongest in mountain terrain. Both signature products were computed from this file.
- **Coverage relationships**: {n_with:,} of the {n_compiled:,} compiled gages have Daymet
  forcing (the Daymet table additionally covers {n_extra} basins whose gages were not
  compiled); compiled gage{'s' if len(no_forcing) != 1 else ''} without forcing: {nf_txt}. The
  {n_meta_rows:,}-row metadata table records why non-compiled gages were excluded
  (`processing_status`).
{gagesets_bullet}
## License and citation

License: CC-BY 4.0. Please cite the HISSS data paper and this resource (DOI: to be
assigned at publication). Source data credits: USGS National Water Information System
(discharge); Water Survey of Canada HYDAT (discharge, station metadata); Daymet V4 R1
(Thornton et al. 2022, ORNL DAAC; basin aggregation with the archived HISSS Daymet
toolchain, fractional cell-coverage weighting over the basin polygons);
GAGES-II (Falcone 2011; US drainage areas).

This resource is part of the HISSS collection on HydroShare:
https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/
(collection DOI upon publication: https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768).
"""
    open(os.path.join(a.out, "hisss_readme_inputs.md"), "w", encoding="utf-8").write(readme)
    # --- abstract (local; pasted into HydroShare by hand)
    ab = open(os.path.join(a.staged_v2, "hisss_abstract_resource3.txt"), encoding="utf-8").read()
    ab = must_replace(ab, "(version 2, aggregated with", "(aggregated with", "abstract version label")
    ab = must_replace(ab, ", with a companion table flagging 101 basins whose outlines are low-confidence or covered by few Daymet cells", "", "abstract flags clause")
    ab = must_replace(ab, "extends it to every compiled basin and to 2025", f"extends it to {n_with:,} of the {n_compiled:,} compiled basins and to 2025", "abstract coverage")
    ab = must_replace(ab, "A data dictionary defines every column of the tables.",
                      ("The gage sets of the two signature products, the gages that qualify under this table but are not included, the ledger of the USGS days retrieved in October 2026 and the list of flagged Daymet basins accompany the tables; a data dictionary defines every column." if comp else
                       "A data dictionary defines every column of the four tables."), "abstract dictionary")
    ab = must_replace(ab, f"for {n_sites:,} gaged basins over calendar years 1980-2025", f"for {n_sites:,} gaged basins over calendar years {c_d0.year}-{c_d1.year}", "abstract Daymet span")
    open(os.path.join(a.out, "hisss_abstract_resource3.txt"), "w", encoding="utf-8").write(ab)
    # --- dictionary: the rows of the four tables, four definitions brought up to date
    dd = pd.read_csv(os.path.join(a.staged_v2, "hisss_input_data_dictionary.csv"), dtype=str, keep_default_na=False)
    dd = dd[dd["table"].isin(["hisss_streamflow_daily.parquet", "hisss_daymet_basin_daily.parquet", "hisss_gage_metadata.csv", "hisss_hydat_interference.csv"])].reset_index(drop=True)
    def setdef(table, col, text=None, extra=None):
        m = (dd["table"] == table) & (dd["column"] == col); assert m.sum() == 1, (table, col)
        dd.loc[m, "definition"] = text if text is not None else dd.loc[m, "definition"] + extra   # extra starts with its own punctuation
    setdef("hisss_streamflow_daily.parquet", "flag", extra=(
        f". A US value is kept when the status letter is A or P and the value is finite; the {n_led:,} approved days carrying a remark code "
        f"(A R, A <, A >, A [4]; {n_led_gages} gages) were retrieved on 26 October 2026 and carry their value; the {n_prem:,} provisional days "
        f"carrying a remark code (P Ice, P Dis, P Rat, …; {n_prem_gages} gages) have NA Q"))
    setdef("hisss_streamflow_daily.parquet", "Date", extra=(
        f". Requests covered {q_d0.isoformat()} – {q_d1.isoformat()} (retrieved 7 February 2026); Canadian records end {can_last.isoformat()} "
        f"(HYDAT release of {HYDAT_RELEASE}, inferred from the retrieval date and the tidyhydat version)"))
    setdef("hisss_streamflow_daily.parquet", "Q", extra=f"; negative values are kept as served ({n_neg:,} rows on {len(neg_gages)} gages)")
    m_fl = (dd["table"] == "hisss_streamflow_daily.parquet") & (dd["column"] == "flag"); d_fl = dd.loc[m_fl, "definition"].item()
    for old_, new_ in (("US rows always carry a status letter", "US rows carry a status letter"),
                       ("(NULL also on ~12,000 US rows)", f"(NULL also on {n_us_null_flag:,} US rows, {'all' if n_us_null_flag_finite == 0 else n_us_null_flag_finite} without a value)")):
        if d_fl.count(old_) != 1:
            raise SystemExit(f"dictionary flag row: anchor not found once: {old_!r}")
        d_fl = d_fl.replace(old_, new_)
    dd.loc[m_fl, "definition"] = d_fl
    setdef("hisss_gage_metadata.csv", "gage_id", text=(
        "Agency station identifier AS RETRIEVED: in this table, every US site number is stored ZERO-STRIPPED (e.g. 1011000 for USGS site 01011000), "
        "while the streamflow and Daymet tables store the zero-padded form, so a direct join fails for most US gages. To join: strip leading zeros "
        "on BOTH sides (never re-pad — US ids are not all 8 digits, and padding to 8 leaves compiled gages unmatched); Canadian identifiers are "
        "unaffected. Read as a character string"))
    rule = f"calendar years of the compiled record with at least {COMPILE_MIN_DAYS} days of non-missing discharge above {COMPILE_MIN_Q} (the compile's year rule; consistent with the streamflow table)"
    setdef("hisss_gage_metadata.csv", "num_years", text=f"Number of {rule}")
    setdef("hisss_gage_metadata.csv", "start_year", text=f"First of the {rule}")
    setdef("hisss_gage_metadata.csv", "end_year", text=f"Last of the {rule}")
    setdef("hisss_gage_metadata.csv", "num_days", text=f"Number of daily rows in the {rule}")
    setdef("hisss_hydat_interference.csv", "STATION_NUMBER", extra=f"; HYDAT release of {HYDAT_RELEASE} (inferred from the retrieval date and the tidyhydat version)")
    for table, col, needles in [("hisss_gage_metadata.csv", "processing_status", [f"success ({st.get('success', 0):,} gages", f"no_data ({st.get('no_data', 0):,})", f"insufficient_years ({st.get('insufficient_years', 0)}", f"processing ({st.get('processing', 0)}"]),
                                ("hisss_gage_metadata.csv", "area_normalized", [f"TRUE ({an_true:,} candidate stations)", f"FALSE ({an_false:,} Canadian candidate stations", f"Only {n_unnorm} of the FALSE stations"]),
                                ("hisss_streamflow_daily.parquet", "Q", [f"for the {n_unnorm} compiled gages with area_normalized = FALSE"]),
                                ("hisss_hydat_interference.csv", "REGULATED", [f"empty for {n_reg_empty} stations"])]:
        m = (dd["table"] == table) & (dd["column"] == col); d_ = dd.loc[m, "definition"].item()
        for nd in needles:
            if nd not in d_:
                raise SystemExit(f"dictionary {table}/{col}: expected '{nd}' in the definition — recount")
    if comp:
        setdef("hisss_streamflow_daily.parquet", "flag", extra=" (ledger: hisss_streamflow_qualifier_corrections.csv)")
        dd0 = pd.read_csv(os.path.join(a.staged_v2, "hisss_input_data_dictionary.csv"), dtype=str, keep_default_na=False)
        flag_rows = dd0[dd0["table"] == "hisss_daymet_basin_flags.csv"]
        if len(flag_rows) != 19:
            raise SystemExit(f"flags dictionary rows: expected 19, found {len(flag_rows)}")
        comp_rows_dd = [
            ("hisss_streamflow_qualifier_corrections.csv", "gage_id", "-", "USGS site number, zero-padded (joins hisss_streamflow_daily.parquet)"),
            ("hisss_streamflow_qualifier_corrections.csv", "date", "date", "Calendar date of the retrieved day"),
            ("hisss_streamflow_qualifier_corrections.csv", "old_flag", "-", "The qualification code of the February 2026 retrieval (A R, A <, A >, A [4]); the flag column of the streamflow table keeps it"),
            ("hisss_streamflow_qualifier_corrections.csv", "old_Q", "mm/day", "The February value (always NA)"),
            ("hisss_streamflow_qualifier_corrections.csv", "value_today_cfs", "ft3/s", "Daily mean discharge served by NWIS on 26 October 2026, as a number"),
            ("hisss_streamflow_qualifier_corrections.csv", "value_today_raw", "-", "The served value as a string, before conversion to a number (equal to value_today_cfs on every row)"),
            ("hisss_streamflow_qualifier_corrections.csv", "qualifiers_today", "-", "Qualification code served on 26 October 2026"),
            ("hisss_streamflow_qualifier_corrections.csv", "status_today", "-", "Status letter of qualifiers_today: A (approved) or P (provisional)"),
            ("hisss_streamflow_qualifier_corrections.csv", "Q_new", "mm/day", "The stored value = value_today_cfs × conversion (every ledger gage is area-normalized)"),
            ("hisss_streamflow_qualifier_corrections.csv", "conversion", "(mm/day)/(ft3/s)", "The ingestion's ft3/s → mm/day factor from the metadata drainage area: 86400 s/day × 1e6 / (3280.84^3 × basin_area_km2) = 2.4465753 / basin_area_km2 (3280.84 ft per km)"),
            ("hisss_streamflow_qualifier_corrections.csv", "basin_area_km2", "km2", "The drainage area used in the conversion (= basin_area of the gage-metadata table)"),
            ("hisss_streamflow_qualifier_corrections.csv", "unit_today", "-", "Unit of the served value (ft3/s on every row)"),
            ("hisss_streamflow_qualifier_corrections.csv", "method_id", "-", "NWIS method (time-series) identifier the value was read from (the unqualified block where a site publishes several)"),
            ("hisss_streamflow_qualifier_corrections.csv", "retrieved", "UTC timestamp", "Time of the NWIS request"),
            ("hisss_streamflow_qualifier_corrections.csv", "status", "-", "restored (every row: the day was eligible and its value stored)"),
            ("hisss_gage_selection_wy1993_2025.txt / hisss_gage_selection_wy1980_2025.txt", "gage_id", "-", "One gage id per line, zero-padded: the gage set of the product, in the product table's row order"),
            ("hisss_would_qualify_not_added.csv", "gage_id, window", "-", "Gage (zero-padded) and product window (1993-2025 or 1980-2025)"),
            ("hisss_would_qualify_not_added.csv", "in_product", "boolean", "Whether the gage is in that product (False on every row: these gages are not in the products)"),
            ("hisss_would_qualify_not_added.csv", "baseline_n_valid / final_n_valid / gained_years / n_gained / lost_years / n_lost", "count, years", "Valid water years under the products' screening without / with the days retrieved on 26 October 2026, and the years gained / lost"),
            ("hisss_would_qualify_not_added.csv", "baseline_total_possible / final_total_possible / baseline_frac / final_frac / baseline_qualifies / final_qualifies", "count, fraction, boolean", "The gage's possible water years inside the window (gage-specific, bounded by its record), the valid fraction of them and the qualification verdict (20 years and 60 %) without / with those days"),
            ("hisss_would_qualify_not_added.csv", "n_restored_days_in_window / n_source_unavailable_in_window / baseline_valid_years / final_valid_years", "count, years", "Retrieved days inside the window; masked days no current record could supply (0 everywhere); the valid years without / with those days"),
        ]
        dd = pd.concat([dd, pd.DataFrame(comp_rows_dd, columns=["table", "column", "units", "definition"]), flag_rows], ignore_index=True)
        setdef("hisss_daymet_basin_flags.csv", "canon_id", extra=" (zero-stripped for USGS ids — join on site_id, not on this column)")
    dd.to_csv(os.path.join(a.out, "hisss_input_data_dictionary.csv"), index=False)
    with open(os.path.join(a.out, "hisss_rename_manifest_inputs.csv"), "w") as f:
        f.write("original_name,deposit_name\n")
        for o, d_ in manifest:
            f.write(f"{o},{d_}\n")
    print(f"staged {a.out} ({'six files + companions' if comp else 'six files'}): streamflow {n_rows_q:,} rows / {n_gages_q:,} gages / {n_na_q:,} NA rows ({q_d0} – {q_d1}; Canadian records end {can_last}); "
          f"Daymet {n_rows_c:,} rows / {n_sites:,} sites ({c_d0} – {c_d1}); metadata {n_meta_rows:,} rows, compiled {n_compiled:,}, without forcing {no_forcing}; "
          f"ledger {n_led:,} days / {n_led_gages} gages ({n_prov_today} provisional today); provisional-remark days {n_prem:,} / {n_prem_gages} gages; "
          f"HYDAT {n_hydat:,} stations ({n_reg_empty} without regulation status); Daymet flags {n_flag} ({n_lowconf} low-confidence, {n_lps} few cells, {n_both} both; in products {n_p1} / {n_p2} / either {n_either}); "
          f"negative Q rows {n_neg:,} on {len(neg_gages)} gages; US rows with NULL flag {n_us_null_flag:,} ({n_us_null_flag_finite:,} with a value); direct-id matches {n_direct:,}; pad-to-8 misses {n_pad}; unnormalized in products {n_unnorm_p}")


# ----------------------------------------------------------------------------------------------- check
REQUIRED = {
    "signatures": ["hisss_signatures_{tag}.csv", "hisss_signatures_annual_{tag}.parquet", "hisss_run_provenance_{tag}.json", "hisss_run_merge_{tag}.json",
                   "hisss_run_notes_{tag}.md", "hisss_readme_{tag}.md", "hisss_data_dictionary.csv", "hisss_signature_categories.csv",
                   "hisss_signature_explorer_{tag}.html", "hisss_signature_explorer_{tag}_annual", "hisss_gage_selection_{tag}.txt",
                   "hisss_would_qualify_not_added_{tag}.csv", "hisss_validation_gates_{tag}.txt", "hisss_validation_annual_values_{tag}.txt",
                   "hisss_validation_qualification_audit_{tag}.txt", "hisss_run_log_{tag}.log", "hisss_rename_manifest_{tag}.csv"],
    "inputs_minimal": ["hisss_streamflow_daily.parquet", "hisss_gage_metadata.csv", "hisss_daymet_basin_daily.parquet", "hisss_hydat_interference.csv",
                       "hisss_input_data_dictionary.csv", "hisss_readme_inputs.md", "hisss_abstract_resource3.txt", "hisss_rename_manifest_inputs.csv"],
    "inputs_companions": ["hisss_streamflow_qualifier_corrections.csv", "hisss_daymet_basin_flags.csv", "hisss_would_qualify_not_added.csv",
                          "hisss_gage_selection_wy1993_2025.txt", "hisss_gage_selection_wy1980_2025.txt"],
    "inputs_full": ["hisss_streamflow_daily.parquet", "hisss_streamflow_daily_as_compiled_feb2026.parquet", "hisss_streamflow_qualifier_corrections.csv",
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
    req = [r.format(tag=tag) for r in REQUIRED[a.kind if a.kind == "signatures" else f"inputs_{a.deposit}"]]
    if a.kind == "inputs" and a.deposit == "minimal" and getattr(a, "companions", False):
        req += REQUIRED["inputs_companions"]
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
            for frag in (b"version 2", b"version 3", b"_v3_", b"26oct2026"):
                if frag in t:
                    bad.append(f"explorer: public text contains {frag!r}")
            side = os.path.join(a.out, f"hisss_signature_explorer_{tag}_annual")
            if os.path.isdir(side) and not any(f.endswith(".js") for f in os.listdir(side)):
                bad.append("explorer: sidecar folder holds no .js files")
            if os.path.exists(side + ".zip"):
                import zipfile
                names = set(zipfile.ZipFile(side + ".zip").namelist())
                want_names = {f"hisss_signature_explorer_{tag}_annual/{f}" for f in os.listdir(side) if not f.startswith(".")}
                if names != want_names:
                    bad.append(f"explorer: sidecar zip entries differ from the folder ({len(names)} vs {len(want_names)})")
    present = set(os.listdir(a.out))
    for doc in [n for n in present if n.startswith(("hisss_readme", "hisss_run_notes", "hisss_abstract"))]:
        txt = open(os.path.join(a.out, doc), encoding="utf-8").read()
        for phrase in ("version 2", "version 3", "previous build", "supersede", "this version", "earlier release", "2026 release", "regenerated", "undeposited", "not deposited"):
            if phrase in txt.lower():
                bad.append(f"{doc}: public text contains {phrase!r}")
        for ref in sorted(set(re.findall(r"`(hisss_[A-Za-z0-9_./-]+)`", txt))):
            base_ = ref.rstrip("/")
            sibling = os.path.join(os.path.dirname(os.path.abspath(a.out)), "resource3_input_data")
            if base_ not in present and not base_.startswith("hisss_staging_manifest") and not os.path.exists(os.path.join(sibling, base_)):
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
    s.add_argument("--staged-v2", required=True); s.add_argument("--out", required=True)
    s.add_argument("--deposit", choices=["minimal", "full"], default="minimal",
                   help="minimal (default, user decision 2026-10-29): the README describes the seven deposited files (summary CSV, annual parquet, dictionary, categories, explorer HTML + sidecar zip) and the sidecar zip is built; full: every staged file")
    s.add_argument("--explorer-title-drop", default=" (version 2, October 2026: corrected USGS qualifier mask and climate frame)",
                   help="exact fragment removed from the explorer's <title>, <h1> and meta.title (empty string = none)")
    s.set_defaults(fn=stage_signatures)
    i = sp.add_parser("inputs"); i.add_argument("--rerun-dir", required=True); i.add_argument("--staged-v2", required=True)
    i.add_argument("--base-parquet", default=None, help="the February 2026 compile: required for --deposit full; with minimal, an optional NA-row accounting check")
    i.add_argument("--out", required=True)
    i.add_argument("--deposit", choices=["minimal", "full"], default="minimal",
                   help="minimal (default, user decision 2026-10-29): the six files of the August 2026 upload, replaced in place; full: the February compile, ledger, raw responses, control report, oracle, selection lists and recipe beside them")
    i.add_argument("--companions", action="store_true",
                   help="with --deposit minimal: also stage the correction ledger, the Daymet flags table, the two gage-selection lists and the would-qualify list (the files the signature readmes point readers to). DECISION (user, 2026-10-29): the R3 deposit carries them — eleven files")
    i.set_defaults(fn=lambda a: stage_inputs_minimal(a) if a.deposit == "minimal" else stage_inputs(a))
    c = sp.add_parser("check"); c.add_argument("--out", required=True); c.add_argument("--kind", required=True, choices=["signatures", "inputs"])
    c.add_argument("--window", type=int, default=None, choices=[1993, 1980]); c.add_argument("--deposit", choices=["minimal", "full"], default="minimal")
    c.add_argument("--companions", action="store_true"); c.set_defaults(fn=check)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
