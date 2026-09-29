---
name: sync-docs
description: Session-start sync of the two collaborative Google Docs that govern this project — the signature guidelines (methodology ground truth) and the HISSS manuscript draft — against their repo snapshots, followed by the reconciliation review. Run at the start of every session, and whenever asked to sync, re-check or reconcile the guidelines or the manuscript.
argument-hint: [guidelines|manuscript|all]
---

# Sync the guidelines and manuscript Google Docs

The guidelines doc is the methodology ground truth (declared 2026-08-31): domain
experts write plain-English definitions and QA rules there, and the code implements
them. The manuscript (Scientific Data, "HISSS", submission target 2026-11-09) must
stay consistent with both the code and the repo docs. The Google Doc URLs live in
`sync_google_docs.py` next to this file (single source).

## 1. Fetch and diff — do NOT read the snapshots first

```bash
python .claude/skills/sync-docs/sync_google_docs.py             # both docs
python .claude/skills/sync-docs/sync_google_docs.py --doc guidelines
```

The script fetches the published pages (raw HTML, the `<div id="contents">` block,
scripts stripped), diffs them paragraph-by-paragraph against
`docs/SIGNATURE_GUIDELINES.md` and `docs/MANUSCRIPT_DRAFT.md`, and prints only the
changed paragraphs. Exit 0 = unchanged, 1 = changed, 2 = fetch error. Never use
WebFetch for this (it paraphrases) and never hand-edit a snapshot body.

**Both unchanged** → report "guidelines unchanged, manuscript unchanged" in one line
and stop.

## 2. Guidelines changed

1. `python .claude/skills/sync-docs/sync_google_docs.py --doc guidelines --write`
   rewrites the body and the `Last synced` date. Then edit the header blockquote of
   `docs/SIGNATURE_GUIDELINES.md` by hand to a SHORT note of what changed (keep the
   header under ~15 lines; the detail goes in the reconciliation log below).
2. Classify each changed paragraph: new/changed signature definition, QA flag or
   threshold, function requirement or parameter, or a colleague's comment.
3. Verify every changed claim against canonical Julia (`julia/src/`) and
   `config/signatures_config.json`, then decide the direction of fix: doc-side (code is
   right → relay note for the user; the doc cannot be edited from here) or code-side
   (implementation TODO).
4. Add a dated entry, newest first, with checkboxes, to
   `docs/reconciliation/guidelines_todos.md`.
5. Report: "Guidelines document has X new/changed items. Would you like to review and
   implement?" — grouped doc-side vs code-side.
6. When implementation is approved: Julia first → `julia --project=julia
   julia/test/runtests.jl` → `/run-benchmark` → port to Python and rpkg
   (`/cross-language-alignment` if outputs diverge) → tick the checkbox → CHANGELOG entry.

## 3. Manuscript changed

1. `python .claude/skills/sync-docs/sync_google_docs.py --doc manuscript --write`;
   shorten the header note of `docs/MANUSCRIPT_DRAFT.md` by hand as above.
2. Reconciliation review of every changed METHODS claim, in three directions:
   - manuscript vs code (Julia canonical + config): correctly and precisely implemented?
   - manuscript vs repo docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`,
     `docs/SIGNATURE_GUIDELINES.md`): correctly and precisely documented?
   - direction of fix: code/docs right and manuscript wrong → relay to the user;
     manuscript states the agreed methodology and code/docs lag → implementation TODO.
3. Counts in the text are checkable against the product CSVs and
   `docs/benchmarks/qualification_census.jl`; say explicitly which numbers were NOT
   verified.
4. Add a dated entry, newest first, to `docs/reconciliation/manuscript_log.md`.
5. Present the findings grouped by direction of fix.

## 4. Afterwards

- Update `docs/STATUS.md` if a decision, an open item, or a sync date changed.
- CHANGELOG gets an entry only when code or repo docs change as a result; the sync
  itself is recorded in the reconciliation files, not the changelog.
