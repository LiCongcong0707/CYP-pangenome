血红素结合区（heme-binding region）、I-helix、PERF 基序和 EXXR 基序存在与位置

```python
import re
from Bio import SeqIO

# 血红素结合区（PROSITE PS00086）及较宽松的备用写法
HEME_STRICT = re.compile(r"[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAP]")
HEME_LOOSE  = re.compile(r"F..G.[RHPT].C.[GAP]")
EXXR   = re.compile(r"E..R")
PERF   = re.compile(r"P.[EKR][FH]")
IHELIX = re.compile(r"[AG]..?[ED]T[TS]")   # (A/G)GX(E/D)T(T/S)

# 各基序相对血红素区起点的上游距离窗口 (min, max)，偏宽松，可按需调整
WIN_PERF   = (10, 100)
WIN_EXXR   = (50, 250)
WIN_IHELIX = (80, 300)

def in_window(pattern, seq, heme_pos, win):
    lo, hi = max(0, heme_pos - win[1]), max(0, heme_pos - win[0])
    return pattern.search(seq, lo, hi) is not None

def check_cyp_protein(fasta_file, out_tsv="cyp_check.tsv"):
    rows = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        seq = str(rec.seq).upper().replace("*", "")
        m = HEME_STRICT.search(seq)
        level = "strict"
        if not m:
            m = HEME_LOOSE.search(seq)
            level = "loose"

        if not m:
            rows.append((rec.id, len(seq), "-", "-", "-", "-", "-", "CHECK (no heme motif)"))
            continue

        pos = m.start()
        has_perf = in_window(PERF, seq, pos, WIN_PERF)
        has_exxr = in_window(EXXR, seq, pos, WIN_EXXR)
        has_ihel = in_window(IHELIX, seq, pos, WIN_IHELIX)

        if level == "strict" and has_exxr and (has_perf or has_ihel):
            verdict = "PASS"
        elif level == "strict" and (has_perf or has_exxr or has_ihel):
            verdict = "LIKELY (strict heme, partial others)"
        elif has_perf or has_exxr or has_ihel:
            verdict = "CHECK (partial)"
        else:
            verdict = "CHECK (heme only)"
        rows.append((rec.id, len(seq), f"{pos+1}:{m.group()}({level})",
                     int(has_perf), int(has_exxr), int(has_ihel), "-", verdict))

    header = ["ID", "Length", "Heme_pos:motif", "PERF", "EXXR", "I-helix", "-", "Verdict"]
    with open(out_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    n_pass = sum(r[-1] == "PASS" for r in rows)
    print(f"共 {len(rows)} 条序列，PASS: {n_pass}，需人工检查: {len(rows) - n_pass}")
    print(f"详细结果已写入 {out_tsv}")
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("用法: python3 check_cyp.py <input.fasta> [output.tsv]")
        sys.exit(1)
    out = sys.argv[2] if len(sys.argv) > 2 else "cyp_check.tsv"
    check_cyp_protein(sys.argv[1], out)
# check_cyp_protein("your_hmmer_hits.fasta")
```

用法：python3 check_cyp.py
