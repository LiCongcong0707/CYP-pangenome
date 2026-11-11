```python
#!/usr/bin/env python3
"""
Summarise FIMO promoter scans (one sub-directory per sample).

Expected layout (run from ~/cyp/promoter):
    ./<sample>.bed                       promoter regions (chr start0 end gene . strand)
    ./fimo_out/<sample>/fimo.tsv         FIMO output (genomic coordinates, i.e. --parse-genomic-coord)
    ./fimo_out/<sample>/best_site.narrowPeak

Usage:
    python fimo_summary.py --promoter-dir . --fimo-dir fimo_out
    python fimo_summary.py --qvalue 0.05 --bin 200 --family-map my_family_map.tsv
    # use JASPAR's own TFClass family for the *.family.* tables:
    python fimo_summary.py --jaspar-meta jaspar_meta.tsv --cluster-by jaspar-family

Only the Python standard library is required.
"""
import argparse
import bisect
import os
import re
import sys
from collections import defaultdict


# --------------------------------------------------------------------------
# motif family clustering
# --------------------------------------------------------------------------
VAR_RE = re.compile(r"\(var\.\d+\)$")
TAIL_RE = re.compile(r"[-_.]?\d+$")


def auto_family(alt):
    """BPC5 -> BPC, ARF2(var.2) -> ARF, DOF5.8 -> DOF, MYB3R1 -> MYB3R."""
    a = VAR_RE.sub("", alt).strip()
    prev = None
    while prev != a:
        prev = a
        b = TAIL_RE.sub("", a)
        if b:
            a = b
    return a


def load_family_map(path):
    m = {}
    if path:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and not line.startswith("#"):
                    m[p[0]] = p[1]
    return m


def get_family(alt, famap):
    if alt in famap:
        return famap[alt]
    auto = auto_family(alt)
    return famap.get(auto, auto)


def load_jaspar_meta(path):
    """TSV from fetch_jaspar_meta.py: matrix_id, name, class, family, ..."""
    m = {}
    if path:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            ci, fi = header.index("class"), header.index("family")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) > max(ci, fi):
                    m[p[0]] = (p[ci], p[fi])
    return m


class Clusterer:
    """Decides the grouping label ('family' columns) and holds JASPAR class/family."""

    def __init__(self, famap, meta, mode):
        self.famap, self.meta, self.mode = famap, meta, mode
        if mode in ("jaspar-family", "jaspar-class") and not meta:
            sys.exit("--cluster-by %s needs --jaspar-meta" % mode)

    def jaspar(self, mid):
        return self.meta.get(mid, ("", ""))

    def label(self, mid, alt):
        cls, fam = self.jaspar(mid)
        if self.mode == "jaspar-family" and fam:
            return fam
        if self.mode == "jaspar-class" and cls:
            return cls
        return get_family(alt, self.famap)  # fall back to name-based clustering


# --------------------------------------------------------------------------
# promoter bed
# --------------------------------------------------------------------------
class PromoterIndex:
    def __init__(self, path):
        self.info = {}  # gene -> (chrom, start0, end, strand)
        recs = defaultdict(list)
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.rstrip("\n").split("\t")
                chrom, s0, e, gene = p[0], int(p[1]), int(p[2]), p[3]
                strand = p[5] if len(p) > 5 and p[5] in ("+", "-") else "+"
                self.info[gene] = (chrom, s0, e, strand)
                recs[chrom].append((s0 + 1, e, gene))  # 1-based inclusive
        self.recs, self.starts, self.maxlen = {}, {}, {}
        for chrom, lst in recs.items():
            lst.sort()
            self.recs[chrom] = lst
            self.starts[chrom] = [r[0] for r in lst]
            self.maxlen[chrom] = max(r[1] - r[0] + 1 for r in lst)

    def query(self, chrom, a, b):
        """genes whose promoter overlaps [a, b] (1-based inclusive)"""
        lst = self.recs.get(chrom)
        if not lst:
            return []
        starts = self.starts[chrom]
        j = bisect.bisect_right(starts, b) - 1
        out = []
        while j >= 0 and starts[j] >= a - self.maxlen[chrom]:
            s, e, g = lst[j]
            if e >= a:
                out.append(g)
            j -= 1
        return sorted(out)

    def rel_pos(self, gene, start, stop):
        """Distance of site centre to the TSS (negative = upstream).
        Assumes the promoter bed ends at the TSS (upstream-only region)."""
        chrom, s0, e, strand = self.info[gene]
        centre = (start + stop) / 2.0
        if strand == "+":
            return int(round(centre - e))
        return int(round((s0 + 1) - centre))


# --------------------------------------------------------------------------
# fimo.tsv
# --------------------------------------------------------------------------
def read_fimo(path):
    rows, header = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if header is None:
                header = p
                continue
            rows.append(dict(zip(header, p)))
    return rows


def parse_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build_hits(rows, idx, cl, pmax, qmax):
    """Return (hits, fimo_lookup, n_unassigned). One hit per (site, gene)."""
    hits, lookup, unassigned = [], {}, 0
    for d in rows:
        chrom = d["sequence_name"]
        start, stop = int(d["start"]), int(d["stop"])
        alt = d.get("motif_alt_id", "")
        if alt in ("", "."):
            alt = d["motif_id"]
        lookup[(chrom, start, d["strand"], d["motif_id"])] = (alt, start, stop)
        pv, qv = parse_float(d.get("p-value")), parse_float(d.get("q-value"))
        if pmax is not None and pv is not None and pv > pmax:
            continue
        if qmax is not None and qv is not None and qv > qmax:
            continue
        genes = idx.query(chrom, start, stop)
        if not genes:
            unassigned += 1
            continue
        fam = cl.label(d["motif_id"], alt)
        jc, jf = cl.jaspar(d["motif_id"])
        for g in genes:
            hits.append(dict(
                gene=g, chrom=chrom, start=start, stop=stop, strand=d["strand"],
                motif_id=d["motif_id"], alt=alt, family=fam, jclass=jc, jfamily=jf,
                score=d.get("score", ""),
                pvalue=d.get("p-value", ""), qvalue=d.get("q-value", ""),
                rel=idx.rel_pos(g, start, stop)))
    hits.sort(key=lambda h: (h["gene"], h["start"], h["stop"], h["alt"]))
    return hits, lookup, unassigned


# --------------------------------------------------------------------------
# narrowPeak annotation (task 1 + 2)
# --------------------------------------------------------------------------
def annotate_narrowpeak(path, out, lookup, idx, cl):
    n_total = n_nogene = 0
    with open(path) as f, open(out, "w") as o:
        o.write("\t".join(["gene", "chrom", "chromStart", "chromEnd", "motif_id",
                           "motif_alt_id", "family", "jaspar_class", "jaspar_family", "score",
                           "strand", "signalValue",
                           "pValue", "qValue", "peak"]) + "\n")
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            chrom, ns, ne, mid, strand = p[0], int(p[1]), int(p[2]), p[3], p[5]
            hit = None
            for s in (ns, ns + 1, ns - 1):  # tolerate 0/1-based offset
                hit = lookup.get((chrom, s, strand, mid))
                if hit:
                    break
            if hit:
                alt, a, b = hit
            else:
                alt, a, b = mid, ns, ne - 1
            genes = idx.query(chrom, a, b)
            n_total += 1
            if not genes:
                n_nogene += 1
            fam = cl.label(mid, alt)
            jc, jf = cl.jaspar(mid)
            o.write("\t".join([",".join(genes) if genes else "NA", chrom, p[1], p[2], mid,
                               alt, fam, jc, jf, p[4], strand] + p[6:10]) + "\n")
    return n_total, n_nogene


# --------------------------------------------------------------------------
# tables
# --------------------------------------------------------------------------
def write_table1(path, hits, key, sample=None):
    with open(path, "w") as o:
        cols = ["gene", "chr", "start", "stop", "motif_alt_id"]
        o.write("\t".join((["sample"] if sample else []) + cols) + "\n")
        for h in hits:
            row = [h["gene"], h["chrom"], str(h["start"]), str(h["stop"]), h[key]]
            o.write("\t".join(([sample] if sample else []) + row) + "\n")


def write_table2(path, hits, key, all_genes, sample=None):
    per = defaultdict(list)
    for h in hits:
        per[h["gene"]].append("%s:%d-%d" % (h[key], h["start"], h["stop"]))
    with open(path, "w") as o:
        o.write("\t".join((["sample"] if sample else []) +
                          ["gene", "motif_alt_id:start-stop", "motif_count"]) + "\n")
        for g in all_genes:
            sites = per.get(g, [])
            o.write("\t".join(([sample] if sample else []) +
                              [g, ";".join(sites), str(len(sites))]) + "\n")


def write_detail(path, hits, sample=None):
    with open(path, "w") as o:
        cols = ["gene", "chr", "start", "stop", "strand", "motif_id", "motif_alt_id",
                "family", "jaspar_class", "jaspar_family", "score", "p-value", "q-value",
                "rel_pos_to_TSS"]
        o.write("\t".join((["sample"] if sample else []) + cols) + "\n")
        for h in hits:
            row = [h["gene"], h["chrom"], str(h["start"]), str(h["stop"]), h["strand"],
                   h["motif_id"], h["alt"], h["family"], h["jclass"], h["jfamily"],
                   h["score"], h["pvalue"], h["qvalue"], str(h["rel"])]
            o.write("\t".join(([sample] if sample else []) + row) + "\n")


# --------------------------------------------------------------------------
# extra analyses
# --------------------------------------------------------------------------
def merge_count(intervals):
    """number of non-overlapping sites after merging overlapping intervals"""
    n, cur_end = 0, None
    for s, e in sorted(intervals):
        if cur_end is None or s > cur_end:
            n += 1
            cur_end = e
        else:
            cur_end = max(cur_end, e)
    return n


def gene_family_counts(hits):
    """{(gene, family): (n_hits, n_nonredundant_sites)}"""
    iv = defaultdict(list)
    for h in hits:
        iv[(h["gene"], h["family"])].append((h["start"], h["stop"]))
    return {k: (len(v), merge_count(v)) for k, v in iv.items()}


def write_family_summary(path, per_sample, with_sample):
    """per_sample: list of (sample, hits, n_genes_total)"""
    with open(path, "w") as o:
        o.write("\t".join((["sample"] if with_sample else []) +
                          ["family", "n_hits", "n_sites_nonredundant", "n_genes",
                           "pct_genes_with_site"]) + "\n")
        for sample, hits, n_genes in per_sample:
            gfc = gene_family_counts(hits)
            agg = defaultdict(lambda: [0, 0, 0])
            for (g, fam), (nh, nr) in gfc.items():
                agg[fam][0] += nh
                agg[fam][1] += nr
                agg[fam][2] += 1
            for fam, (nh, nr, ng) in sorted(agg.items(), key=lambda x: -x[1][2]):
                pct = 100.0 * ng / n_genes if n_genes else 0
                o.write("\t".join(([sample] if with_sample else []) +
                                  [fam, str(nh), str(nr), str(ng), "%.1f" % pct]) + "\n")


def write_matrix(path, per_sample, all_genes_by_sample, with_sample):
    fams = sorted({h["family"] for _, hits, _ in per_sample for h in hits})
    with open(path, "w") as o:
        o.write("\t".join((["sample"] if with_sample else []) + ["gene"] + fams) + "\n")
        for sample, hits, _ in per_sample:
            gfc = gene_family_counts(hits)
            for g in all_genes_by_sample[sample]:
                row = [str(gfc.get((g, f), (0, 0))[1]) for f in fams]
                o.write("\t".join(([sample] if with_sample else []) + [g] + row) + "\n")


def write_position_dist(path, per_sample, binsize, with_sample):
    with open(path, "w") as o:
        o.write("\t".join((["sample"] if with_sample else []) +
                          ["family", "bin_to_TSS", "n_hits"]) + "\n")
        for sample, hits, _ in per_sample:
            cnt = defaultdict(int)
            for h in hits:
                r = h["rel"]
                if r > 0:
                    lab = "downstream(>0)"
                else:
                    b = int((-r) // binsize)
                    lab = "-%d~-%d" % ((b + 1) * binsize, b * binsize)
                cnt[(h["family"], lab)] += 1
            for (fam, lab), n in sorted(cnt.items()):
                o.write("\t".join(([sample] if with_sample else []) +
                                  [fam, lab, str(n)]) + "\n")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--promoter-dir", default=".", help="directory holding <sample>.bed")
    ap.add_argument("--fimo-dir", default="fimo_out")
    ap.add_argument("--samples", nargs="*", help="default: every sub-dir with fimo.tsv")
    ap.add_argument("--bed-suffix", default=".bed")
    ap.add_argument("--family-map", help="optional 2-column tsv: name<TAB>family")
    ap.add_argument("--jaspar-meta", help="TSV from fetch_jaspar_meta.py (class/family)")
    ap.add_argument("--cluster-by", default="name",
                    choices=["name", "jaspar-family", "jaspar-class"],
                    help="grouping used for the *.family.* tables (default: strip digits)")
    ap.add_argument("--pvalue", type=float, help="extra p-value cut-off")
    ap.add_argument("--qvalue", type=float, help="extra q-value cut-off")
    ap.add_argument("--bin", type=int, default=200, help="bin size for TSS distance")
    a = ap.parse_args()

    famap = load_family_map(a.family_map)
    cl = Clusterer(famap, load_jaspar_meta(a.jaspar_meta), a.cluster_by)
    samples = a.samples or sorted(
        d for d in os.listdir(a.fimo_dir)
        if os.path.isfile(os.path.join(a.fimo_dir, d, "fimo.tsv")))
    if not samples:
        sys.exit("no sample with fimo.tsv found in %s" % a.fimo_dir)

    per_sample, genes_by_sample = [], {}
    for s in samples:
        sdir = os.path.join(a.fimo_dir, s)
        bed = os.path.join(a.promoter_dir, s + a.bed_suffix)
        if not os.path.isfile(bed):
            print("[skip] %s: %s not found" % (s, bed), file=sys.stderr)
            continue
        idx = PromoterIndex(bed)
        rows = read_fimo(os.path.join(sdir, "fimo.tsv"))
        hits, lookup, unassigned = build_hits(rows, idx, cl, a.pvalue, a.qvalue)
        all_genes = sorted(idx.info)

        np_in = os.path.join(sdir, "best_site.narrowPeak")
        msg = ""
        if os.path.isfile(np_in):
            n, nog = annotate_narrowpeak(
                np_in, os.path.join(sdir, "best_site.annotated.narrowPeak"),
                lookup, idx, cl)
            msg = " narrowPeak=%d(no gene: %d)" % (n, nog)

        write_table1(os.path.join(sdir, "table1.sites.tsv"), hits, "alt")
        write_table2(os.path.join(sdir, "table2.gene_summary.tsv"), hits, "alt", all_genes)
        write_table1(os.path.join(sdir, "table1.sites.family.tsv"), hits, "family")
        write_table2(os.path.join(sdir, "table2.gene_summary.family.tsv"), hits,
                     "family", all_genes)
        write_detail(os.path.join(sdir, "sites.detail.tsv"), hits)
        one = [(s, hits, len(all_genes))]
        write_family_summary(os.path.join(sdir, "family_summary.tsv"), one, False)
        write_matrix(os.path.join(sdir, "gene_family_matrix.tsv"), one,
                     {s: all_genes}, False)
        write_position_dist(os.path.join(sdir, "position_distribution.tsv"), one,
                            a.bin, False)

        n_with = len({h["gene"] for h in hits})
        print("[%s] promoters=%d with_sites=%d hits=%d outside_promoters=%d%s" %
              (s, len(all_genes), n_with, len(hits), unassigned, msg))
        per_sample.append((s, hits, len(all_genes)))
        genes_by_sample[s] = all_genes

    if not per_sample:
        sys.exit("nothing processed")

    # ---- combined tables in fimo_dir (task 4) ----
    out = a.fimo_dir
    all_t1 = [(s, h) for s, hits, _ in per_sample for h in hits]
    for key, tag in (("alt", ""), ("family", ".family")):
        with open(os.path.join(out, "all_samples.table1.sites%s.tsv" % tag), "w") as o:
            o.write("sample\tgene\tchr\tstart\tstop\tmotif_alt_id\n")
            for s, h in all_t1:
                o.write("\t".join([s, h["gene"], h["chrom"], str(h["start"]),
                                   str(h["stop"]), h[key]]) + "\n")
        with open(os.path.join(out, "all_samples.table2.gene_summary%s.tsv" % tag), "w") as o:
            o.write("sample\tgene\tmotif_alt_id:start-stop\tmotif_count\n")
            for s, hits, _ in per_sample:
                per = defaultdict(list)
                for h in hits:
                    per[h["gene"]].append("%s:%d-%d" % (h[key], h["start"], h["stop"]))
                for g in genes_by_sample[s]:
                    o.write("\t".join([s, g, ";".join(per.get(g, [])),
                                       str(len(per.get(g, [])))]) + "\n")
    with open(os.path.join(out, "all_samples.sites.detail.tsv"), "w") as o:
        o.write("sample\tgene\tchr\tstart\tstop\tstrand\tmotif_id\tmotif_alt_id\tfamily\t"
                "jaspar_class\tjaspar_family\tscore\tp-value\tq-value\trel_pos_to_TSS\n")
        for s, h in all_t1:
            o.write("\t".join([s, h["gene"], h["chrom"], str(h["start"]), str(h["stop"]),
                               h["strand"], h["motif_id"], h["alt"], h["family"],
                               h["jclass"], h["jfamily"], h["score"], h["pvalue"], h["qvalue"], str(h["rel"])]) + "\n")
    write_family_summary(os.path.join(out, "all_samples.family_summary.tsv"), per_sample, True)
    write_matrix(os.path.join(out, "all_samples.gene_family_matrix.tsv"), per_sample,
                 genes_by_sample, True)
    write_position_dist(os.path.join(out, "all_samples.position_distribution.tsv"),
                        per_sample, a.bin, True)
    print("done -> %s/all_samples.*" % out)


if __name__ == "__main__":
    main()
```


```PYTHON
##获取JASPAR的元数据
#!/usr/bin/env python3
"""
Build a matrix_id -> class/family table for the motifs in a JASPAR MEME file,
using the JASPAR REST API (one request per matrix, so the version always matches).

Usage:
    python fetch_jaspar_meta.py JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt jaspar_meta.tsv

Needs internet access (run it on a machine/server that can reach jaspar.elixir.no).
Results are cached in <out>.cache.json, so an interrupted run can simply be restarted.
"""
import json
import os
import sys
import time
import urllib.request

HOSTS = ["https://jaspar.elixir.no", "https://jaspar.genereg.net"]


def meme_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith("MOTIF"):
                p = line.split()
                ids.append((p[1], p[2] if len(p) > 2 else ""))
    return ids


def fetch(mid, retries=3):
    last = None
    for host in HOSTS:
        url = "%s/api/v1/matrix/%s/?format=json" % (host, mid)
        for k in range(retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "fetch_jaspar_meta"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    return json.loads(r.read().decode())
            except Exception as e:  # network / HTTP errors
                last = e
                time.sleep(1 + k)
    raise RuntimeError("%s: %s" % (mid, last))


def as_text(x):
    if x is None:
        return ""
    if isinstance(x, (list, tuple)):
        return "::".join(str(i) for i in x)  # dimers have 2 entries
    return str(x)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    meme, out = sys.argv[1], sys.argv[2]
    cache_path = out + ".cache.json"
    cache = json.load(open(cache_path)) if os.path.isfile(cache_path) else {}
    motifs = meme_ids(meme)
    failed = []
    for i, (mid, name) in enumerate(motifs, 1):
        if mid in cache:
            continue
        try:
            d = fetch(mid)
            cache[mid] = {k: as_text(d.get(k)) for k in ("name", "class", "family", "tax_group")}
        except RuntimeError as e:
            failed.append(mid)
            print("[fail]", e, file=sys.stderr)
        if i % 50 == 0:
            json.dump(cache, open(cache_path, "w"))
            print("%d/%d" % (i, len(motifs)), file=sys.stderr)
        time.sleep(0.05)
    json.dump(cache, open(cache_path, "w"))

    with open(out, "w") as o:
        o.write("matrix_id\tname\tclass\tfamily\ttax_group\n")
        for mid, name in motifs:
            c = cache.get(mid)
            if c:
                o.write("\t".join([mid, c["name"] or name, c["class"], c["family"],
                                   c["tax_group"]]) + "\n")
    print("written %d/%d motifs -> %s" % (len(motifs) - len(failed), len(motifs), out))
    if failed:
        print("missing:", " ".join(failed), file=sys.stderr)


if __name__ == "__main__":
    main()
```