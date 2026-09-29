#!/usr/bin/env python3
"""Fetch the two published Google Docs that govern this project and diff them
against their repo snapshots at paragraph level.

    python .claude/skills/sync-docs/sync_google_docs.py            # check both, print changes
    python .claude/skills/sync-docs/sync_google_docs.py --doc guidelines
    python .claude/skills/sync-docs/sync_google_docs.py --write    # also overwrite the snapshot bodies

Why a script: the WebFetch tool paraphrases through a small model, and reading
both snapshots into context costs ~20k tokens even when nothing changed. This
script does the fetch, the text extraction (the ``<div id="contents">`` block,
scripts/styles stripped) and the diff deterministically, so the session only
sees the paragraphs that actually changed.

Exit codes: 0 = every requested doc unchanged, 1 = at least one doc changed,
2 = fetch/parse error. Stdlib only, so it runs on the Windows laptop and the
Mac without a virtualenv.

Snapshot layout (both files): a hand-written header block, a line that is
exactly ``---``, then the extracted body. ``--write`` replaces the body and the
date in the header's ``**Last synced**: YYYY-MM-DD`` field only; the rest of the
header note is for Claude to edit by hand.
"""
from __future__ import annotations

import argparse
import datetime as dt
import difflib
import html
import re
import sys
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

DOCS = {
    "guidelines": {
        "title": "Signature guidelines (methodology ground truth)",
        "url": "https://docs.google.com/document/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub",
        "snapshot": ROOT / "docs" / "SIGNATURE_GUIDELINES.md",
    },
    "manuscript": {
        "title": "HISSS manuscript draft",
        "url": "https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub",
        "snapshot": ROOT / "docs" / "MANUSCRIPT_DRAFT.md",
    },
}

BLOCK_TAGS = {"p", "h1", "h2", "h3", "h4", "h5", "h6", "li", "tr"}
HEADING_LEVEL = {"h1": 1, "h2": 2, "h3": 3, "h4": 4, "h5": 5, "h6": 6}


class ContentsExtractor(HTMLParser):
    """Collect (heading_level, text) blocks from inside <div id="contents">."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_contents = False
        self.div_depth = 0          # depth relative to the contents div
        self.skip_depth = 0         # >0 while inside <script>/<style>
        self.block_tag: str | None = None
        self.buf: list[str] = []
        self.blocks: list[tuple[int, str]] = []

    def handle_starttag(self, tag, attrs):
        if tag == "div":
            if not self.in_contents and dict(attrs).get("id") == "contents":
                self.in_contents = True
                self.div_depth = 1
                return
            if self.in_contents:
                self.div_depth += 1
        if not self.in_contents:
            return
        if tag in ("script", "style"):
            self.skip_depth += 1
            return
        if self.skip_depth:
            return
        if tag in BLOCK_TAGS:
            self._flush()
            self.block_tag = tag
        elif tag in ("td", "th") and self.block_tag == "tr" and self.buf:
            self.buf.append(" | ")
        elif tag == "br":
            # a soft line break inside a block is a paragraph boundary in the
            # snapshots (keeps multi-line formula blocks as separate paragraphs)
            self._split()

    def handle_endtag(self, tag):
        if not self.in_contents:
            return
        if tag in ("script", "style"):
            self.skip_depth = max(0, self.skip_depth - 1)
            return
        if tag == "div":
            self.div_depth -= 1
            if self.div_depth == 0:
                self._flush()
                self.in_contents = False
            return
        if self.skip_depth:
            return
        if tag in BLOCK_TAGS and tag == self.block_tag:
            self._flush()

    def handle_data(self, data):
        if self.in_contents and not self.skip_depth and self.block_tag:
            self.buf.append(data)

    def _split(self):
        if self.block_tag is None:
            return
        text = normalize("".join(self.buf))
        if text:
            self.blocks.append((HEADING_LEVEL.get(self.block_tag, 0), text))
        self.buf = []

    def _flush(self):
        if self.block_tag is None:
            return
        text = normalize("".join(self.buf))
        if text:
            self.blocks.append((HEADING_LEVEL.get(self.block_tag, 0), text))
        self.buf = []
        self.block_tag = None


def normalize(s: str) -> str:
    s = html.unescape(s).replace("\xa0", " ")
    return re.sub(r"\s+", " ", s).strip()


def fetch(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (streamflowsignatures sync)"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read()
    return raw.decode("utf-8", errors="replace")


def extract_blocks(page_html: str) -> list[tuple[int, str]]:
    p = ContentsExtractor()
    p.feed(page_html)
    if not p.blocks:
        raise RuntimeError('no text found inside <div id="contents"> - page layout changed?')
    return p.blocks


def split_snapshot(path: Path) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8")
    m = re.search(r"^---\s*$", text, flags=re.M)
    if not m:
        raise RuntimeError(f"{path}: no '---' header/body separator")
    return text[: m.end()], text[m.end():]


def snapshot_paragraphs(body: str) -> list[str]:
    out = []
    for para in re.split(r"\n\s*\n", body):
        t = normalize(" ".join(line.strip() for line in para.splitlines()))
        t = re.sub(r"^#{1,6}\s+", "", t)
        if t:
            out.append(t)
    return out


def render_body(blocks: list[tuple[int, str]]) -> str:
    lines = []
    for level, text in blocks:
        lines.append(("#" * level + " " + text) if level else text)
    return "\n\n" + "\n\n".join(lines) + "\n"


def report_diff(name: str, old: list[str], new: list[str], context: int, width: int) -> int:
    sm = difflib.SequenceMatcher(a=old, b=new, autojunk=False)
    changes = [op for op in sm.get_opcodes() if op[0] != "equal"]
    print(f"[{name}] snapshot={len(old)} paragraphs, live={len(new)} paragraphs, changed regions={len(changes)}")
    for k, (tag, i1, i2, j1, j2) in enumerate(changes, 1):
        print(f"  --- region {k}: {tag} (snapshot paras {i1 + 1}-{i2}, live paras {j1 + 1}-{j2})")
        for c in range(max(0, i1 - context), i1):
            print(f"      {clip(old[c], width)}")
        for i in range(i1, i2):
            print(f"    - {clip(old[i], width)}")
        for j in range(j1, j2):
            print(f"    + {clip(new[j], width)}")
    return len(changes)


def clip(s: str, width: int) -> str:
    return s if len(s) <= width else s[: width - 3] + "..."


def update_last_synced(header: str, today: str) -> str:
    return re.sub(r"(\*\*Last synced\*\*:\s*)\d{4}-\d{2}-\d{2}", rf"\g<1>{today}", header, count=1)


def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):  # Windows consoles default to cp1252
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--doc", choices=[*DOCS, "all"], default="all")
    ap.add_argument("--write", action="store_true", help="overwrite the snapshot body + Last synced date when changed")
    ap.add_argument("--force-write", action="store_true", help="rewrite the snapshot body even if unchanged (reformat)")
    ap.add_argument("--context", type=int, default=1, help="unchanged paragraphs to print before each region")
    ap.add_argument("--width", type=int, default=400, help="max characters printed per paragraph")
    ap.add_argument("--save-html", type=Path, help="directory to dump the raw fetched HTML into (debugging)")
    args = ap.parse_args(argv)

    names = list(DOCS) if args.doc == "all" else [args.doc]
    today = dt.date.today().isoformat()
    any_changed = False
    for name in names:
        spec = DOCS[name]
        try:
            page = fetch(spec["url"])
            blocks = extract_blocks(page)
        except Exception as exc:  # noqa: BLE001
            print(f"[{name}] ERROR fetching/parsing: {exc}", file=sys.stderr)
            return 2
        if args.save_html:
            args.save_html.mkdir(parents=True, exist_ok=True)
            (args.save_html / f"{name}.html").write_text(page, encoding="utf-8")
        header, body = split_snapshot(spec["snapshot"])
        old = snapshot_paragraphs(body)
        new = [t for _, t in blocks]
        n = report_diff(name, old, new, args.context, args.width)
        changed = n > 0
        any_changed |= changed
        print(f"[{name}] {'CHANGED' if changed else 'UNCHANGED'} vs snapshot {spec['snapshot'].relative_to(ROOT)}")
        if (changed and args.write) or args.force_write:
            spec["snapshot"].write_text(update_last_synced(header, today) + render_body(blocks), encoding="utf-8", newline="\n")
            print(f"[{name}] snapshot body rewritten, Last synced -> {today} (edit the header note by hand)")
    return 1 if any_changed else 0


if __name__ == "__main__":
    sys.exit(main())
