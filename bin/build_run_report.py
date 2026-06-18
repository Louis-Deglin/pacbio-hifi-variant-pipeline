#!/usr/bin/env python3
# Assemble the run-level aggregated report from the per-interval fragments:
#   - general_summary.html : every (sample x interval) fragment concatenated, with a
#     plain-text table of contents (sample/interval names, no links).
#
# Reads a manifest TSV (one row per interval: sample<TAB>chrom<TAB>interval). The
# fragment filenames are reconstructed from those fields using the same naming
# convention as render_interval.py / IntervalReport, so no fragile filename parsing
# is needed. Pure Python stdlib. ASCII source.

import argparse
import html
import os
import re

CSS = """
body{font-family:Arial,Helvetica,sans-serif;margin:20px;color:#222;}
h1{border-bottom:2px solid #444;padding-bottom:4px;}
h2{margin-top:28px;border-bottom:1px solid #999;padding-bottom:3px;}
h3{margin-top:18px;}
table{border-collapse:collapse;margin:10px 0;font-size:13px;}
th,td{border:1px solid #ccc;padding:3px 7px;text-align:left;vertical-align:top;}
th{background:#f0f0f0;}
.counts{background:#f8f8f8;padding:8px 12px;border:1px solid #ddd;display:inline-block;}
details{margin:3px 0;}
summary{cursor:pointer;color:#06c;}
pre{background:#f8f8f8;border:1px solid #ddd;padding:8px;overflow:auto;font-size:12px;}
.disabled{color:#888;font-style:italic;}
.note{color:#666;font-size:12px;margin-top:6px;}
.toc{background:#f8f8f8;border:1px solid #ddd;padding:10px 16px;display:inline-block;}
.toc div{margin:2px 0;}
.toc a{color:#06c;text-decoration:none;}
.toc a:hover{text-decoration:underline;}
.toc .sub{margin:1px 0 6px 20px;font-size:12px;}
.toc .sub a{margin-right:12px;color:#39c;}
section{margin-bottom:40px;}
.seq{max-width:240px;overflow-x:auto;white-space:nowrap;display:block;}
.search{margin:8px 6px 2px 0;padding:3px 6px;width:300px;font-size:13px;}
.tablewrap{margin:8px 0;}
.toggle{cursor:pointer;background:#eee;border:1px solid #bbb;border-radius:3px;padding:3px 10px;font-size:13px;}
.toggle:hover{background:#e0e0e0;}
.pager{margin:6px 0;font-size:13px;color:#444;}
.pager button{cursor:pointer;margin:0 4px;}
"""

# Same table behaviour (toggle / AND search / pagination) as render_interval.py, so
# the per-fragment controls work inside the aggregated page too. Keep in sync.
TABLE_JS = r"""
var PAGE_SIZE = 50;
var tableState = {};

function rowMatches(text, tokens){
  for(var j=0;j<tokens.length;j++){
    if(text.indexOf(tokens[j]) === -1){ return false; }
  }
  return true;
}

function applyTable(id){
  var tb = document.getElementById(id);
  if(!tb || !tb.tBodies.length){ return; }
  var input = document.getElementById(id + '-search');
  var tokens = input ? input.value.toLowerCase().split(/\s+/).filter(Boolean) : [];
  var rows = tb.tBodies[0].rows;
  var matched = [];
  for(var i=0;i<rows.length;i++){
    if(rowMatches(rows[i].textContent.toLowerCase(), tokens)){
      matched.push(rows[i]);
    } else {
      rows[i].style.display = 'none';
    }
  }
  var st = tableState[id] || (tableState[id] = {page:0});
  var pages = Math.max(1, Math.ceil(matched.length / PAGE_SIZE));
  if(st.page >= pages){ st.page = pages - 1; }
  if(st.page < 0){ st.page = 0; }
  for(var k=0;k<matched.length;k++){
    matched[k].style.display =
      (k >= st.page*PAGE_SIZE && k < (st.page+1)*PAGE_SIZE) ? '' : 'none';
  }
  renderPager(id, matched.length, pages, st.page);
}

function renderPager(id, total, pages, page){
  var pg = document.getElementById(id + '-pager');
  if(!pg){ return; }
  if(pages <= 1){ pg.textContent = total + ' row(s)'; return; }
  var start = total ? page*PAGE_SIZE + 1 : 0;
  var end = Math.min((page+1)*PAGE_SIZE, total);
  pg.innerHTML =
    `<button onclick="pageTable('${id}',-1)" ${page<=0?'disabled':''}>Prev</button> ` +
    `Rows ${start}-${end} of ${total} (page ${page+1}/${pages}) ` +
    `<button onclick="pageTable('${id}',1)" ${page>=pages-1?'disabled':''}>Next</button>`;
}

function pageTable(id, delta){
  var st = tableState[id] || (tableState[id] = {page:0});
  st.page += delta;
  applyTable(id);
}

function searchTable(id){
  var st = tableState[id] || (tableState[id] = {page:0});
  st.page = 0;
  applyTable(id);
}

function toggleTable(id){
  var body = document.getElementById(id + '-body');
  var btn = document.getElementById(id + '-btn');
  if(!body){ return; }
  if(body.style.display === 'none'){
    body.style.display = '';
    if(btn){ btn.textContent = 'Hide table'; }
    applyTable(id);
  } else {
    body.style.display = 'none';
    if(btn){ btn.textContent = 'Show table'; }
  }
}
"""


def esc(value):
    return html.escape("" if value is None else str(value))


def anchor_id(sample, chrom, interval):
    # Must match render_interval.anchor_id so the plan links resolve to the
    # <section id="..."> emitted in each fragment.
    raw = "%s-%s-%s" % (sample, chrom, interval)
    return re.sub(r"[^A-Za-z0-9_-]", "-", raw)


def read_manifest(path):
    rows = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rows.append((parts[0], parts[1], parts[2]))
    # Stable, predictable order.
    rows.sort()
    return rows


def fragment_name(sample, chrom, interval):
    return "%s_fragment_chr%s_%s.html" % (sample, chrom, interval)


def build_summary(rows, run_name):
    toc = ['<div class="toc"><b>Contents</b>']
    bodies = []
    for sample, chrom, interval in rows:
        aid = anchor_id(sample, chrom, interval)
        label = "%s &mdash; chr%s:%s" % (esc(sample), esc(chrom), esc(interval))
        # Clickable plan entry: jumps to the matching <section id> in this page, with
        # nested links to each section (the -snv-h/-sv-h/-meth-h header ids always
        # exist in the fragment, even when SV/methylation are disabled).
        toc.append(
            '<div><a href="#%s">%s</a>'
            '<div class="sub">'
            '<a href="#%s-snv-h">Variants</a>'
            '<a href="#%s-sv-h">Structural variants</a>'
            '<a href="#%s-meth-h">Methylation</a>'
            "</div></div>" % (aid, label, aid, aid, aid)
        )
        fname = fragment_name(sample, chrom, interval)
        if os.path.isfile(fname):
            with open(fname, "rt", encoding="utf-8", errors="replace") as fh:
                bodies.append(fh.read())
        else:
            bodies.append('<section><h1>%s</h1>'
                          '<p class="disabled">Fragment missing.</p></section>' % label)
    toc.append("</div>")

    title = "General summary &mdash; %s" % esc(run_name) if run_name else "General summary"
    return (
        "<!DOCTYPE html>\n<html><head><meta charset=\"utf-8\">"
        "<title>%s</title>\n<style>%s</style>\n<script>%s</script></head>\n<body>\n"
        "<h1>%s</h1>\n%s\n%s\n</body></html>\n"
        % (title, CSS, TABLE_JS, title, "\n".join(toc), "\n".join(bodies))
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--run-name", default="")
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    rows = read_manifest(args.manifest)

    with open(args.out_summary, "wt", encoding="utf-8") as fh:
        fh.write(build_summary(rows, args.run_name))


if __name__ == "__main__":
    main()
