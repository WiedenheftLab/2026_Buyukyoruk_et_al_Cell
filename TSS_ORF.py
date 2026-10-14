#!/software/mambaforge/bin/python
import argparse
from typing import Optional, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# --------------------------------------------------------------
# Tunable parameters
# --------------------------------------------------------------
OPTIMAL_SPACER            = 17      # ideal bp between end of −35 and start of −10
SPACER_PENALTY_WEIGHT     = 0.2     # penalty per bp deviation from OPTIMAL_SPACER

OPT_START_DIST            = 20      # ideal bp from TSS to start‐codon
START_DIST_PENALTY_WEIGHT = 0.25    # penalty per bp deviation from OPT_START_DIST

FIRST_REPEAT_PENALTY_WEIGHT = 1.0   # penalty if the ORF start codon falls in the first repeat
# --------------------------------------------------------------

# Consensus motifs & codons
CONS_35      = "TTGACA"
CONS_10      = "TATAAT"
CONS_RBS     = "AGGAGG"
START_CODONS = ("ATG", "GTG", "TTG")
STOP_CODONS  = ("TAA", "TAG", "TGA")

# Start‐codon intrinsic scores
START_WEIGHTS = {
    "ATG": 5,   # strong
    "GTG": 3,   # moderate
    "TTG": 1    # weak
}

# Weights for the −35 motif: first 3 positions count double
WEIGHTS_35 = [2, 2, 2, 1, 1, 1]


def score35(box35: str) -> int:
    return sum(w for w, a, b in zip(WEIGHTS_35, box35, CONS_35) if a == b)


def score10(box10: str) -> int:
    return sum(1 for a, b in zip(box10, CONS_10) if a == b)


def score_rbs(window: str) -> tuple[int, Optional[int]]:
    best, offset = -1, None
    for i in range(len(window) - len(CONS_RBS) + 1):
        sc = sum(1 for a, b in zip(window[i:i+6], CONS_RBS) if a == b)
        if sc > best:
            best, offset = sc, i
    return best, offset


def find_promoters(seq: str) -> list[dict]:
    hits = []
    L = len(seq)
    for i in range(L - 5):
        seq35 = seq[i:i+6]
        sc35  = score35(seq35)
        for spacer in (15, 16, 17, 18, 19):
            j = i + 6 + spacer
            if j + 6 > L:
                continue
            seq10      = seq[j:j+6]
            sc10       = score10(seq10)
            raw_score  = sc35 + sc10
            δ_spacer   = abs(spacer - OPTIMAL_SPACER)
            pen_spacer = δ_spacer * SPACER_PENALTY_WEIGHT
            hits.append({
                "pos35": i+1, "end35": i+6, "seq35": seq35, "sc35": sc35,
                "spacer": spacer, "delta_spacer": δ_spacer, "penalty_spacer": pen_spacer,
                "pos10": j+1, "end10": j+6, "seq10": seq10, "sc10": sc10,
                "raw_prom_score": raw_score,
                "prom_score": raw_score - pen_spacer
            })
    return hits


def find_start_codons(seq: str, tss0: int) -> list[dict]:
    L, cands = len(seq), []
    for cod in START_CODONS:
        for i in range(tss0, L-2):
            if seq[i:i+3] != cod:
                continue
            ws, we     = max(0, i-20), max(0, i-5)
            window     = seq[ws:we]
            rbs_sc, off = score_rbs(window)
            if off is not None:
                rbs1    = ws + off + 1
                seq_rbs = seq[rbs1-1:rbs1+5]
                end_rbs = rbs1 + 5
            else:
                rbs1 = end_rbs = seq_rbs = None
            cands.append({
                "start1": i+1, "end_start": i+3,
                "codon": cod,
                "rbs1": rbs1, "end_rbs": end_rbs,
                "seq_rbs": seq_rbs, "rbs_score": rbs_sc
            })
    return sorted(cands, key=lambda x: (x["rbs_score"], -x["start1"]), reverse=True)


def extract_orf(seq: str, start0: int) -> Optional[str]:
    cods = []
    for i in range(start0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c in STOP_CODONS:
            return None
        cods.append(c)
    return "".join(cods)


def translate_dna(dna: str) -> str:
    return str(Seq(dna).translate(to_stop=False))


def trim_partial_repeat(aa: str) -> str:
    n   = len(aa)
    lps = [0]*n
    for i in range(1, n):
        j = lps[i-1]
        while j > 0 and aa[i] != aa[j]:
            j = lps[j-1]
        if aa[i] == aa[j]:
            j += 1
        lps[i] = j
    p = n - lps[-1]
    if lps[-1] > 0 and p < n:
        return aa[: (n//p)*p ]
    return aa


def list_stopless_frames(seq: str) -> List[int]:
    """
    Return list of reading frames (0,1,2) that contain no in-frame stop codons
    across the entire sequence.
    """
    frames = []
    for frame in (0,1,2):
        has_stop = False
        for i in range(frame, len(seq)-2, 3):
            if seq[i:i+3] in STOP_CODONS:
                has_stop = True
                break
        if not has_stop:
            frames.append(frame)
    return frames


def predict_orfs(seq: str, repeat_len: int, allowed_frames: List[int]) -> list[dict]:
    scenarios = []
    for prom in find_promoters(seq):
        tss   = prom["end10"] + 10
        tss0  = tss - 1
        for sc in find_start_codons(seq, tss0):
            frame = (sc["start1"] - 1) % 3
            if frame not in allowed_frames:
                continue

            orf_nt = extract_orf(seq, sc["start1"]-1)
            if orf_nt is None:
                continue

            dist        = sc["start1"] - tss
            δ_start     = abs(dist - OPT_START_DIST)
            pen_start   = δ_start * START_DIST_PENALTY_WEIGHT

            in_first    = sc["start1"] <= repeat_len
            pen_first   = FIRST_REPEAT_PENALTY_WEIGHT if in_first else 0.0

            sc_start    = START_WEIGHTS[sc["codon"]]
            cum_score   = (
                prom["prom_score"]
              + sc["rbs_score"]
              + sc_start
              - pen_start
              - pen_first
            )

            scenarios.append({
                **prom, **sc,
                "tss":            tss,
                "orf_seq":        orf_nt,
                "orf_len":        len(orf_nt),
                "score_start":    sc_start,
                "start_dist":     dist,
                "delta_start":    δ_start,
                "penalty_start":  pen_start,
                "in_first_repeat":in_first,
                "penalty_first":  pen_first,
                "frame":          frame,
                "cum_score":      cum_score
            })

    return sorted(scenarios,
                  key=lambda s: (s["cum_score"], s["orf_len"], s["tss"]),
                  reverse=True)


def has_stopless_frame(seq: str) -> bool:
    return bool(list_stopless_frames(seq))


def main():
    p = argparse.ArgumentParser(
        description="Predict TSS‑driven ORFs on an N‑concatemer of each sequence’s reverse complement"
    )
    p.add_argument("-i", "--input",   required=True, help="Input FASTA file")
    p.add_argument("-n", "--concat",  type=int, default=1,
                   help="Repeat each reverse‑complement N times before prediction")
    p.add_argument("-v", "--verbose", type=int, default=0,
                   help="1=print all scenarios (with full ORF) to stdout")
    p.add_argument("-o", "--out",     required=True,
                   help="Write top ORF(s) as amino‑acid FASTA to this file")
    args = p.parse_args()

    total_records = sum(1 for _ in SeqIO.parse(args.input, "fasta"))
    out_recs = []

    for record in tqdm(SeqIO.parse(args.input, "fasta"),
                       total=total_records,
                       desc="Processing records"):
        rc         = record.seq.reverse_complement()
        concatseq  = str(rc).upper() * args.concat
        repeat_len = len(rc)

        allowed_frames = list_stopless_frames(concatseq)
        if not allowed_frames:
            continue

        scens = predict_orfs(concatseq, repeat_len, allowed_frames)

        if args.verbose >= 1:
            print(f">{record.id}  (revcomp ×{args.concat})")
            if not scens:
                print("  (no valid scenarios)\n")
            else:
                for i, s in enumerate(scens, 1):
                    print(f"Scenario {i}: cumulative_score={s['cum_score']:.1f}")
                    # ... [rest of verbose print unchanged] ...

        if scens:
            best = scens[0]
            aa   = translate_dna(best["orf_seq"])
            aa   = trim_partial_repeat(aa)
            rec = SeqRecord(
                Seq(aa),
                id=f"{record.id}_NEO",
                description=f"score={best['cum_score']:.1f};len={len(aa)}"
            )
            out_recs.append(rec)

    if out_recs:
        from Bio.SeqIO import FastaIO
        fasta_out = FastaIO.FastaWriter(args.out, wrap=None)
        fasta_out.write_file(out_recs)


if __name__ == "__main__":
    main()

