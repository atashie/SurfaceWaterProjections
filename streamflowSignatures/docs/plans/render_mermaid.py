#!/usr/bin/env python
"""Render a mermaid .mmd file to a high-resolution PNG via Playwright Chromium."""
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

SCRATCH = Path(__file__).parent
src_path = Path(sys.argv[1])
out_path = Path(sys.argv[2])
scale = int(sys.argv[3]) if len(sys.argv) > 3 else 3

diagram = src_path.read_text(encoding="utf-8")
mermaid_js = (SCRATCH / "mermaid.min.js").read_text(encoding="utf-8")

# Styling injected here rather than as an in-source %%{init}%% directive:
# GitHub's pinned mermaid crashes at render time on init directives combined
# with labeled edges (mermaid-js/mermaid#6452/#6022), so the committed source
# stays directive-free and the PNG render carries the theme config instead.
# Concrete font face because the headless shell mishandles ui-sans-serif/
# system-ui stacks (falls back to serif).
MERMAID_CONFIG = {
    "startOnLoad": False,
    "securityLevel": "loose",
    "theme": "base",
    "themeVariables": {
        "fontFamily": "Helvetica Neue, Helvetica, Arial, sans-serif",
        "fontSize": "13px",
        "lineColor": "#7f8ea4",
    },
    "flowchart": {"htmlLabels": True, "nodeSpacing": 40, "rankSpacing": 52, "curve": "basis"},
}

html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
  body {{ margin: 0; background: #ffffff; }}
  /* Headless-shell quirk: mermaid's injected font rules never reach the
     foreignObject label divs, which then inherit the browser default (Times).
     Force one concrete face everywhere, before mermaid measures text. */
  body, body * {{ font-family: "Helvetica Neue", Helvetica, Arial, sans-serif !important; }}
  #wrap {{ display: inline-block; padding: 28px; background: #ffffff; }}
</style>
<script>{mermaid_js}</script>
</head><body>
<div id="wrap"><div id="target"></div></div>
<script>
  const src = {json.dumps(diagram)};
  mermaid.initialize({json.dumps(MERMAID_CONFIG)});
  mermaid.render("theSvg", src).then(({{ svg }}) => {{
    document.getElementById("target").innerHTML = svg;
    const el = document.querySelector("#target svg");
    const vb = el.viewBox.baseVal;
    el.style.maxWidth = "none";
    el.setAttribute("width", Math.ceil(vb.width));
    el.setAttribute("height", Math.ceil(vb.height));
    document.title = "RENDERED";
  }}).catch(e => {{ document.title = "ERROR: " + e.message; }});
</script>
</body></html>"""

html_path = SCRATCH / "mermaid_render.html"
html_path.write_text(html, encoding="utf-8")

with sync_playwright() as p:
    browser = p.chromium.launch()
    ctx = browser.new_context(device_scale_factor=scale, viewport={"width": 1600, "height": 1200})
    page = ctx.new_page()
    page.goto(html_path.as_uri())
    page.wait_for_function("document.title === 'RENDERED' || document.title.startsWith('ERROR')", timeout=30000)
    title = page.title()
    if title.startswith("ERROR"):
        print(title, file=sys.stderr)
        sys.exit(1)
    box = page.locator("#wrap").bounding_box()
    page.set_viewport_size({"width": int(box["width"]) + 10, "height": int(box["height"]) + 10})
    page.locator("#wrap").screenshot(path=str(out_path))
    print(f"svg css size: {box['width']:.0f} x {box['height']:.0f} (scale {scale}x)")
    browser.close()

print(f"wrote {out_path}")
