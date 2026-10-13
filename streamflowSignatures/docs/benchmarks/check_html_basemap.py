#!/usr/bin/env python
"""Headless check that an explorer / dashboard HTML shows a basemap (Leaflet tiles load).

Added 2026-10-13 after CARTO's tiles began returning an "API KEY REQUIRED" watermark for
keyless requests: the builders now use `hisssBasemap()` (Esri gray canvas, OpenStreetMap
fallback, layers control), and this script is the post-build check — it loads each page in
headless Chromium (Playwright, in the repo's `.venv`), waits for the tiles, and reports the
tile requests per host and HTTP status, the number of tiles Leaflet marked loaded, the
attribution text and which base layer is selected; a screenshot is written next to each page.

Usage:
    python docs/benchmarks/check_html_basemap.py [--wait=MS] page.html [page2.html ...]

Exit status 1 when any page loads no tile or shows no attribution. A watermark tile cannot
be detected from the DOM (it is a valid image) — look at the screenshot.
"""
import collections
import json
import os
import sys

from playwright.sync_api import sync_playwright


def run(paths, wait_ms=9000, timeout=240000):
    ok = True
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for path in paths:
            path = os.path.abspath(path)
            pg = browser.new_page(viewport={"width": 1280, "height": 820})
            resp, failed, console = collections.Counter(), collections.Counter(), []
            pg.on("response", lambda r: resp.update({(r.url.split("/")[2], r.status): 1}))
            pg.on("requestfailed", lambda r: failed.update({r.url.split("/")[2]: 1}))
            pg.on("console", lambda m: console.append(m.type + ": " + m.text) if m.type in ("error", "warning") else None)
            pg.on("pageerror", lambda e: console.append("pageerror: " + str(e)))
            pg.goto("file://" + path, wait_until="load", timeout=timeout)
            pg.wait_for_timeout(wait_ms)
            info = pg.evaluate("""() => ({
                tiles: document.querySelectorAll('img.leaflet-tile').length,
                loaded: document.querySelectorAll('img.leaflet-tile-loaded').length,
                attribution: (document.querySelector('.leaflet-control-attribution')||{innerText:''}).innerText,
                bases: Array.from(document.querySelectorAll('.leaflet-control-layers-base label')).map(l => l.innerText.trim() + (l.querySelector('input').checked ? ' [on]' : '')),
                controls: Array.from(document.querySelectorAll('.leaflet-control-layers')).map(c => c.parentElement.className)
            })""")
            pg.screenshot(path=path.rsplit(".", 1)[0] + "_basemap_check.png")
            print("PAGE", path)
            print("  dom:", json.dumps(info))
            print("  responses:", dict(resp))
            print("  failed:", dict(failed))
            print("  console:", console[:5])
            if info["loaded"] == 0 or not info["attribution"].strip():
                ok = False
                print("  FAIL: no tile loaded or no attribution")
            pg.close()
        browser.close()
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    wait = 9000
    if args and args[0].startswith("--wait="):
        wait = int(args[0][len("--wait="):])
        args = args[1:]
    if not args:
        print(__doc__)
        sys.exit(2)
    sys.exit(0 if run(args, wait_ms=wait) else 1)
