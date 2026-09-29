---
paths:
  - "CHANGELOG.md"
  - "changelog-old.md"
  - "docs/CHANGELOG_ARCHIVE.md"
  - "docs/STATUS.md"
  - "docs/reconciliation/**"
---

# Changelog, status file and reconciliation logs

Four files, four jobs (restructured 2026-09-29; only STATUS.md is auto-loaded):

| File | Holds | Loaded |
|---|---|---|
| `docs/STATUS.md` | ONE-LINERS with pointers: delivered products, pending user decisions, live known issues, deferred fixes, doc-sync dates | every session (the only import in CLAUDE.md) — keep under 60 lines |
| `CHANGELOG.md` | `[Unreleased]` Planned + Known Issues in full; the current month in full; condensed headline summaries of the two months before it | on demand |
| `changelog-old.md` | full text of closed months, newest first; resolved or superseded `[Unreleased]` entries. Dec 2025 – Apr 2026 detail: `docs/CHANGELOG_ARCHIVE.md` | on demand — never `@`-reference it |
| `docs/reconciliation/guidelines_todos.md`, `manuscript_log.md` | dated sync and reconciliation entries, newest first, checkboxes for open items | by `/sync-docs` |

Conventions:
- Document every code change: date-based sections `[Month Year]`, severity labels
  HIGH/MEDIUM/LOW on fixes, user DECISIONS stated as such with the date. File-level
  change lists belong in `git log`; analysis and benchmark tables belong in
  `docs/SIGNATURES.md` / `docs/DEVELOPMENT.md` and are linked, not re-hosted.
- When a month closes: condense it in CHANGELOG.md to headline bullets + pointers and
  move the full text VERBATIM to `changelog-old.md`; move completed Planned items,
  resolved Known Issues (leave a one-line product caveat if one still applies) and
  superseded reconciliation entries there too.
- When an open item's state changes, update its STATUS.md one-liner in the same edit;
  when it closes, delete the line.
- Guidelines and manuscript syncs are logged in `docs/reconciliation/`, not in
  CHANGELOG — CHANGELOG gets an entry only when code or repo docs change as a result.
- Keep `claude-skill/streamflow-signatures.md` (the user-facing interpretation skill)
  current when methodology, output formats, validation or cross-language status change.
