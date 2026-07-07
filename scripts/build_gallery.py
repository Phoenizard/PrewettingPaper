#!/usr/bin/env python
"""Build a static, offline HTML gallery of our results vs the result/ reference.

Reads $RESULTS_DB/verify (ours, synced from the server) and pairs each case with
result/ (theirs, the local reference). Writes $RESULTS_DB/index.html — open it in
a browser (no server needed): pick chi(stage) / om / chibb to see the two figures
side by side, plus the pw_line.csv path (the data itself is not rendered).

  RESULTS_DB=./database conda run -n numenv python scripts/build_gallery.py

Paths in the HTML are absolute file:// URLs, so it works even when the database
lives on an external drive while result/ stays on the Mac.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import cases  # noqa: E402

RESULTS_DB = os.path.abspath(os.environ.get("RESULTS_DB", os.path.join(ROOT, "database")))
VERIFY = os.path.join(RESULTS_DB, "verify")


def _file_url(path):
    return "file://" + path if os.path.exists(path) else None


def _iter_db_cases():
    if not os.path.isdir(VERIFY):
        return
    for dirpath, _, files in os.walk(VERIFY):
        if "pw_line.csv" in files:
            rel = os.path.relpath(dirpath, VERIFY).split(os.sep)
            if len(rel) == 3:
                yield tuple(rel)


def _load_summary():
    """rel -> {n_pw, status} from database/verify/SUMMARY.csv (if present)."""
    path = os.path.join(VERIFY, "SUMMARY.csv")
    out = {}
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                out[(r["chi_dir"], r["om_dir"], r["chibb_dir"])] = {
                    "n_pw": r.get("n_pw", ""), "status": r.get("status", "")}
    return out


def _link_reference():
    """Symlink $RESULTS_DB/result -> repo result/ so the reference figures live
    inside the gallery's own tree. Browsers (esp. Safari) block file:// images
    outside the HTML's directory, so ours (under the DB) loaded but result/ (a
    sibling of the DB) did not; routing it through this symlink fixes that."""
    link = os.path.join(RESULTS_DB, "result")
    repo_result = os.path.join(ROOT, "result")
    if os.path.isdir(repo_result) and not os.path.lexists(link):
        try:
            os.symlink(repo_result, link)
        except OSError:
            pass


def build():
    _link_reference()
    summ = _load_summary()
    catalog = []
    for rel in sorted(_iter_db_cases()):
        _, surf = cases.parse_case(*rel)
        ours = os.path.join(VERIFY, *rel, "overlay.png")
        # via the in-DB symlink so the URL stays inside the gallery tree
        theirs = os.path.abspath(cases.result_overlay(rel, root=RESULTS_DB))
        meta = summ.get(rel, {})
        catalog.append({
            "chi": rel[0], "om": rel[1], "chibb": rel[2],
            "om1": surf.w1, "om2": surf.w2,
            "cbb1": surf.cbb1, "cbb2": surf.cbb2, "cbb12": surf.cbb12,
            "n_pw": meta.get("n_pw", ""), "status": meta.get("status", ""),
            "ours": _file_url(ours),
            "theirs": _file_url(theirs),
            "data": os.path.join(VERIFY, *rel, "pw_line.csv"),
        })
    out_path = os.path.join(RESULTS_DB, "index.html")
    os.makedirs(RESULTS_DB, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(_HTML.replace("/*CATALOG*/", json.dumps(catalog)))
    n_theirs = sum(1 for c in catalog if c["theirs"])
    print(f"wrote {out_path}  ({len(catalog)} cases, {n_theirs} with result/ reference)")


_HTML = r"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><title>Prewetting 结果画廊</title>
<style>
 body{font-family:-apple-system,Helvetica,Arial,sans-serif;margin:1.2rem;color:#222}
 h1{font-size:1.1rem} .bar{display:flex;gap:.8rem;align-items:center;flex-wrap:wrap;margin:.6rem 0}
 select{font-size:1rem;padding:.2rem} label{font-size:.85rem;color:#555}
 .meta{font-size:.85rem;color:#333;margin:.4rem 0}
 .figs{display:flex;gap:1rem;flex-wrap:wrap} .fig{flex:1;min-width:320px}
 .fig h3{font-size:.9rem;margin:.2rem 0} img{max-width:100%;border:1px solid #ddd;background:#fafafa}
 .miss{color:#a00;font-size:.85rem;padding:2rem;border:1px dashed #ccc;text-align:center}
 code{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px;word-break:break-all}
 .count{color:#888;font-size:.8rem}
</style></head><body>
<h1>Prewetting 结果画廊 <span class="count" id="count"></span></h1>
<div class="bar">
 <label>stage (chi)</label><select id="chi"></select>
 <label>om</label><select id="om"></select>
 <label>chibb</label><select id="chibb"></select>
</div>
<div class="meta" id="meta"></div>
<div class="figs">
 <div class="fig"><h3>ours (我们, 服务器算 → 同步)</h3><div id="ours"></div></div>
 <div class="fig"><h3>result (参考)</h3><div id="theirs"></div></div>
</div>
<div class="meta">数据: <code id="data"></code></div>
<script>
const C = /*CATALOG*/;
const $ = id => document.getElementById(id);
$("count").textContent = "(" + C.length + " cases)";
const uniq = a => [...new Set(a)].sort();
function fill(sel, vals, fmt){ sel.innerHTML=""; vals.forEach(v=>{const o=document.createElement("option");
  o.value=v; o.textContent=fmt?fmt(v):v; sel.appendChild(o);}); }
function fmtOm(om){ const c=C.find(x=>x.om===om); return c?`om1=${c.om1}, om2=${c.om2}`:om; }
function fmtCb(cb){ const c=C.find(x=>x.chibb===cb); return c?`cbb=(${c.cbb1},${c.cbb2},${c.cbb12})`:cb; }
function img(url){ return url?`<img loading="lazy" src="${url}">`:`<div class="miss">缺图</div>`; }
function refreshOm(){ const chi=$("chi").value;
  fill($("om"), uniq(C.filter(x=>x.chi===chi).map(x=>x.om)), fmtOm); refreshCb(); }
function refreshCb(){ const chi=$("chi").value, om=$("om").value;
  fill($("chibb"), uniq(C.filter(x=>x.chi===chi&&x.om===om).map(x=>x.chibb)), fmtCb); render(); }
function render(){ const c=C.find(x=>x.chi===$("chi").value&&x.om===$("om").value&&x.chibb===$("chibb").value);
  if(!c){return;}
  $("meta").innerHTML = `${c.chi} | ${c.om} | ${c.chibb} &nbsp; n_pw=${c.n_pw} status=${c.status}`;
  $("ours").innerHTML = img(c.ours); $("theirs").innerHTML = img(c.theirs);
  $("data").textContent = c.data; }
$("chi").onchange=refreshOm; $("om").onchange=refreshCb; $("chibb").onchange=render;
fill($("chi"), uniq(C.map(x=>x.chi))); refreshOm();
</script></body></html>"""


if __name__ == "__main__":
    build()
