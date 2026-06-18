#!/usr/bin/env python3
# Render one (sample x interval) HTML report from the annotated SNV VCF, the
# annotated SV VCF and the methylation summary. Produces two files:
#   - a standalone HTML page (full <html> doc, inline CSS)
#   - a reusable HTML fragment (just the <section>...</section>, no CSS) that the
#     run-level report concatenates into full_report.html.
#
# Pure Python stdlib (gzip + html). The CSQ subfield order is read from the VCF
# header (##INFO=<ID=CSQ,...Format: Allele|Consequence|IMPACT|...">), so no
# bcftools +split-vep dependency. ASCII source; output stays ASCII-safe.

import argparse
import gzip
import html
import re

# Curated CSQ subfields surfaced as table columns (only those actually present in
# the VCF's CSQ Format are shown). Everything else goes in a collapsible <details>.
CURATED_SNV = ["Consequence", "IMPACT", "SYMBOL", "Gene", "Feature",
               "BIOTYPE", "CANONICAL", "HGVSc", "HGVSp"]
CURATED_SV = ["Consequence", "IMPACT", "SYMBOL", "Gene", "Feature", "BIOTYPE"]

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
.plan{background:#f8f8f8;border:1px solid #ddd;padding:10px 16px;display:inline-block;margin:6px 0;}
.plan div{margin:2px 0;}
.plan a{color:#06c;text-decoration:none;}
.plan a:hover{text-decoration:underline;}
section{margin-bottom:40px;}
.seq{max-width:240px;overflow-x:auto;white-space:nowrap;display:block;}
.search{margin:8px 6px 2px 0;padding:3px 6px;width:300px;font-size:13px;}
.tablewrap{margin:8px 0;}
.toggle{cursor:pointer;background:#eee;border:1px solid #bbb;border-radius:3px;padding:3px 10px;font-size:13px;}
.toggle:hover{background:#e0e0e0;}
.pager{margin:6px 0;font-size:13px;color:#444;}
.pager button{cursor:pointer;margin:0 4px;}
"""

# Client-side table behaviour shared by the standalone and aggregated pages:
#   - toggleTable(id): collapse/expand a table (collapsed by default).
#   - searchTable(id): per-table filter; the query is split on whitespace and EVERY
#     token must be present in a row (AND search). Resets to page 1.
#   - pagination: at most PAGE_SIZE (50) matching rows are shown at once, with
#     Prev/Next controls. Only this table's own tbody is touched (nested <details>
#     tables live in their own tbody, so they're untouched).
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


def table_block(table_id, header_cells, body_rows):
    # Collapsible (closed by default) + searchable (AND) + paginated table wrapper.
    tid = esc(table_id)
    return (
        '<div class="tablewrap">'
        '<button class="toggle" id="%s-btn" onclick="toggleTable(\'%s\')">Show table</button>'
        '<div class="tablebody" id="%s-body" style="display:none">'
        '<input class="search" id="%s-search" onkeyup="searchTable(\'%s\')" '
        'placeholder="Search (space = AND)...">'
        '<table id="%s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>'
        '<div class="pager" id="%s-pager"></div>'
        "</div></div>"
        % (tid, tid, tid, tid, tid, tid, header_cells, body_rows, tid)
    )


def esc(value):
    return html.escape("" if value is None else str(value))


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def anchor_id(sample, chrom, interval):
    raw = "%s-%s-%s" % (sample, chrom, interval)
    return re.sub(r"[^A-Za-z0-9_-]", "-", raw)


def parse_csq_format(headers):
    for line in headers:
        if line.startswith("##INFO=<ID=CSQ"):
            idx = line.find("Format:")
            if idx != -1:
                tail = line[idx + len("Format:"):].strip()
                tail = tail.rstrip(">").rstrip('"').strip()
                return tail.split("|")
    return []


def read_vcf(path):
    headers = []
    records = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                headers.append(line)
                continue
            if not line.strip():
                continue
            records.append(line.split("\t"))
    return headers, records


def info_dict(info):
    out = {}
    if info == ".":
        return out
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def fmt_sample(fmt_col, sample_col):
    keys = fmt_col.split(":")
    vals = sample_col.split(":")
    return dict(zip(keys, vals))


def split_csq(raw_csq, csq_fields):
    """Return a list of transcripts, each a dict {field_name: value}."""
    transcripts = []
    if not raw_csq:
        return transcripts
    for tr in raw_csq.split(","):
        parts = tr.split("|")
        d = {}
        for i, name in enumerate(csq_fields):
            d[name] = parts[i] if i < len(parts) else ""
        transcripts.append(d)
    return transcripts


def csq_details_html(transcripts, csq_fields):
    """Collapsible full-CSQ table (all subfields, all transcripts)."""
    if not transcripts:
        return ""
    head = "".join("<th>%s</th>" % esc(f) for f in csq_fields)
    rows = []
    for t in transcripts:
        cells = "".join("<td>%s</td>" % esc(t.get(f, "")) for f in csq_fields)
        rows.append("<tr>%s</tr>" % cells)
    return ("<details><summary>full VEP (%d transcript(s))</summary>"
            "<table><tr>%s</tr>%s</table></details>"
            % (len(transcripts), head, "".join(rows)))


def classify_gt(gt):
    g = gt.replace("|", "/")
    if g in ("0/0",):
        return "hom_ref"
    if g in ("0/1", "1/0"):
        return "het"
    if g in ("1/1",):
        return "hom_alt"
    return "other"


def build_snv_section(snv_vcf, table_id):
    headers, records = read_vcf(snv_vcf)
    csq_fields = parse_csq_format(headers)
    curated = [f for f in CURATED_SNV if f in csq_fields]

    total = len(records)
    counts = {"hom_ref": 0, "het": 0, "hom_alt": 0, "other": 0}
    dp_sum = 0.0
    dp_n = 0
    rows = []

    for fields in records:
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        info = info_dict(fields[7])
        gt = dp = ps = "."
        if len(fields) >= 10:
            sd = fmt_sample(fields[8], fields[9])
            gt = sd.get("GT", ".")
            dp = sd.get("DP", ".")
            ps = sd.get("PS", ".")
        counts[classify_gt(gt)] += 1
        try:
            dp_sum += float(dp)
            dp_n += 1
        except (ValueError, TypeError):
            pass

        transcripts = split_csq(info.get("CSQ", ""), csq_fields)
        first = transcripts[0] if transcripts else {}
        curated_cells = "".join("<td>%s</td>" % esc(first.get(f, "")) for f in curated)
        details = csq_details_html(transcripts, csq_fields)

        rows.append(
            "<tr><td>%s</td><td>%s</td>"
            '<td><div class="seq">%s</div></td><td><div class="seq">%s</div></td>'
            "<td>%s</td><td>%s</td><td>%s</td><td>%s</td>%s<td>%s</td></tr>"
            % (esc(chrom), esc(pos), esc(ref), esc(alt), esc(gt),
               esc(ps if ps else "."), esc(dp), esc(qual),
               curated_cells, details)
        )

    avg_dp = "%.1f" % (dp_sum / dp_n) if dp_n else "N/A"

    parts = []
    parts.append('<h2 id="%s-h">Single nucleotide variants &amp; indels</h2>'
                 % esc(table_id))
    parts.append(
        '<div class="counts">Total: <b>%d</b> &nbsp;|&nbsp; '
        "Hom-ref (0/0): %d &nbsp;|&nbsp; Het (0/1): %d &nbsp;|&nbsp; "
        "Hom-alt (1/1): %d &nbsp;|&nbsp; Other: %d &nbsp;|&nbsp; "
        "Mean DP: %s</div>"
        % (total, counts["hom_ref"], counts["het"], counts["hom_alt"],
           counts["other"], esc(avg_dp))
    )

    if total == 0:
        parts.append('<p class="disabled">No variants found in this interval.</p>')
        return "\n".join(parts)

    curated_headers = "".join("<th>%s</th>" % esc(f) for f in curated)
    header_cells = ("<th>CHROM</th><th>POS</th><th>REF</th><th>ALT</th>"
                    "<th>GT</th><th>PHASE_SET</th><th>DP</th><th>QUAL</th>"
                    + curated_headers + "<th>VEP detail</th>")
    parts.append(table_block(table_id, header_cells, "\n".join(rows)))
    parts.append(
        '<p class="note">PHASE_SET (PS) is the HiPhase block ID: variants sharing the '
        "same PS are phased together, so their GT (0|1 / 1|0) are comparable (cis/trans). "
        "'.' means unphased. Curated CSQ columns shown; expand 'VEP detail' for all "
        "transcripts/subfields.</p>"
    )
    return "\n".join(parts)


def build_sv_section(sv_vcf, run_sv, table_id):
    h2 = '<h2 id="%s-h">Structural variants (sawfish)</h2>' % esc(table_id)
    if not run_sv:
        return (h2 + '<p class="disabled">Structural variant calling disabled '
                "(enable with --run_sv).</p>")

    headers, records = read_vcf(sv_vcf)
    csq_fields = parse_csq_format(headers)
    curated = [f for f in CURATED_SV if f in csq_fields]

    parts = [h2]
    if not records:
        parts.append('<p class="disabled">No structural variants overlapping this '
                     "region.</p>")
        return "\n".join(parts)

    by_type = {}
    rows = []
    for fields in records:
        pos = fields[1]
        filt = fields[6]
        info = info_dict(fields[7])
        svtype = info.get("SVTYPE", ".")
        endp = info.get("END", ".")
        svlen = info.get("SVLEN", ".")
        gt = "."
        if len(fields) >= 10:
            gt = fmt_sample(fields[8], fields[9]).get("GT", ".")
        by_type[svtype] = by_type.get(svtype, 0) + 1

        transcripts = split_csq(info.get("CSQ", ""), csq_fields)
        first = transcripts[0] if transcripts else {}
        curated_cells = "".join("<td>%s</td>" % esc(first.get(f, "")) for f in curated)
        details = csq_details_html(transcripts, csq_fields)

        rows.append(
            "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
            "<td>%s</td>%s<td>%s</td></tr>"
            % (esc(svtype), esc(pos), esc(endp), esc(svlen), esc(filt),
               esc(gt), curated_cells, details)
        )

    summary = " &nbsp;|&nbsp; ".join("%s: %d" % (esc(t), n) for t, n in sorted(by_type.items()))
    parts.append('<div class="counts">Total: <b>%d</b> &nbsp;|&nbsp; %s</div>'
                 % (len(records), summary))
    curated_headers = "".join("<th>%s</th>" % esc(f) for f in curated)
    header_cells = ("<th>TYPE</th><th>POS</th><th>END</th><th>SVLEN</th>"
                    "<th>FILTER</th><th>GT</th>" + curated_headers + "<th>VEP detail</th>")
    parts.append(table_block(table_id, header_cells, "\n".join(rows)))
    return "\n".join(parts)


def build_meth_section(meth_summary, run_meth, table_id):
    # The SliceMethylation summary is a fixed-width text block: a few header stat
    # lines (CpG count, mean % combined/hap1/hap2), a blank line, then a per-CpG
    # table (POS METH% COV HAP1% HAP2%). Reparse it into the same HTML table shape
    # as the variant sections (counts box + searchable table).
    parts = ['<h2 id="%s-h">Methylation (CpG)</h2>' % esc(table_id)]
    if not run_meth:
        parts.append('<p class="disabled">Methylation profiling disabled '
                     "(enable with --run_methylation).</p>")
        return "\n".join(parts)
    try:
        with open(meth_summary, "rt", encoding="utf-8", errors="replace") as fh:
            lines = fh.read().split("\n")
    except OSError:
        lines = []

    if not any(l.strip() for l in lines):
        parts.append('<p class="disabled">No methylation summary available.</p>')
        return "\n".join(parts)
    if lines[0].startswith("No methylation"):
        parts.append('<p class="disabled">No methylation profile found.</p>')
        return "\n".join(parts)

    # Header stat lines = everything up to the first blank line.
    stat_lines = []
    i = 0
    while i < len(lines) and lines[i].strip():
        stat_lines.append(lines[i].strip())
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    # Skip the column-header line and the dashes separator, keep the data rows.
    data_rows = []
    if i < len(lines):
        for line in lines[i + 2:]:
            if line.strip():
                data_rows.append(line.split())

    parts.append('<div class="counts">'
                 + " &nbsp;|&nbsp; ".join(esc(s) for s in stat_lines)
                 + "</div>")
    if not data_rows:
        return "\n".join(parts)

    header_cells = ("<th>POS</th><th>METH%</th><th>COV</th>"
                    "<th>HAP1%</th><th>HAP2%</th>")
    body_rows = []
    for r in data_rows:
        cells = (r + ["", "", "", "", ""])[:5]
        body_rows.append("<tr>" + "".join("<td>%s</td>" % esc(c) for c in cells) + "</tr>")
    parts.append(table_block(table_id, header_cells, "\n".join(body_rows)))
    return "\n".join(parts)


def build_fragment(sample, chrom, interval, snv_html, sv_html, meth_html):
    aid = anchor_id(sample, chrom, interval)
    title = "%s &mdash; chr%s:%s" % (esc(sample), esc(chrom), esc(interval))
    return (
        '<section id="%s">\n'
        "<h1>%s</h1>\n%s\n%s\n%s\n</section>"
        % (aid, title, snv_html, sv_html, meth_html)
    )


def standalone_page(sample, chrom, interval, fragment):
    title = "Report %s chr%s:%s" % (sample, chrom, interval)
    aid = anchor_id(sample, chrom, interval)
    # Small clickable plan jumping to each section of this report.
    nav = ('<div class="plan"><b>Sections</b>'
           '<div><a href="#%s-snv-h">Single nucleotide variants &amp; indels</a></div>'
           '<div><a href="#%s-sv-h">Structural variants</a></div>'
           '<div><a href="#%s-meth-h">Methylation</a></div></div>'
           % (aid, aid, aid))
    return (
        "<!DOCTYPE html>\n<html><head><meta charset=\"utf-8\">"
        "<title>%s</title>\n<style>%s</style>\n<script>%s</script></head>\n"
        "<body>\n%s\n%s\n</body></html>\n"
        % (esc(title), CSS, TABLE_JS, nav, fragment)
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--interval", required=True)
    ap.add_argument("--snv-vcf", required=True)
    ap.add_argument("--sv-vcf", required=True)
    ap.add_argument("--meth-summary", required=True)
    ap.add_argument("--run-sv", default="false")
    ap.add_argument("--run-methylation", default="false")
    ap.add_argument("--out-standalone", required=True)
    ap.add_argument("--out-fragment", required=True)
    args = ap.parse_args()

    run_sv = args.run_sv == "true"
    run_meth = args.run_methylation == "true"

    aid = anchor_id(args.sample, args.chrom, args.interval)
    snv_html = build_snv_section(args.snv_vcf, aid + "-snv")
    sv_html = build_sv_section(args.sv_vcf, run_sv, aid + "-sv")
    meth_html = build_meth_section(args.meth_summary, run_meth, aid + "-meth")

    fragment = build_fragment(args.sample, args.chrom, args.interval,
                              snv_html, sv_html, meth_html)

    with open(args.out_fragment, "wt", encoding="utf-8") as fh:
        fh.write(fragment + "\n")
    with open(args.out_standalone, "wt", encoding="utf-8") as fh:
        fh.write(standalone_page(args.sample, args.chrom, args.interval, fragment))


if __name__ == "__main__":
    main()
