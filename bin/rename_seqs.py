#!/usr/bin/env python3
import argparse
from ete3 import NCBITaxa
import os

# ----- Helpers -----
def load_taxid_map(path):
    print(f"Loading taxid map from: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        mapping = {parts[0]: parts[1] for parts in (line.strip().split("\t") for line in fh)}
    print(f"Loaded {len(mapping)} entries from taxid map.")
    return mapping

def parse_fasta(path):
    print(f"Parsing FASTA file: {path}")
    header, seq = None, []
    count = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    count += 1
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            count += 1
            yield header, "".join(seq)
    print(f"Parsed {count} sequences from FASTA.")

def extract_gene(desc):
    parts = desc.split("[gene=", 1)
    return parts[1].split("]", 1)[0] if len(parts) > 1 else "unknowngene"

def make_code(scname):
    if not scname:
        return ""
    toks = scname.strip().split()
    n = len(toks)
    lower = [t.lower() for t in toks]
    if n > 3 and 'x' in lower:
        i = lower.index('x')
        left = toks[max(0, i-2):i]
        right = toks[i+1:i+3]
        def side_code(side):
            if not side: return ""
            a = side[0][:2].upper()
            b = side[1][:1].upper() if len(side) > 1 else ""
            return a + b
        return side_code(left) + side_code(right)
    if n == 2:
        return (toks[0][:3] + toks[1][:3]).upper()
    if n == 3:
        return (toks[0][:3] + toks[1][:2] + toks[2][:1]).upper()
    return "".join(toks)[:6].upper()

# ----- Main -----
def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers using taxid mapping and species codes (ETE3), with deduplication.")
    parser.add_argument("--taxidmap", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    ncbi = NCBITaxa()

    out_fasta = f"{args.prefix}_filtered.fasta"
    out_codes = f"{args.prefix}_speciescodes.txt"

    taxid_map = load_taxid_map(args.taxidmap)

    taxid_to_seqs = {}
    total_input_seqs = 0
    total_retained_seqs = 0

    print("Starting FASTA parsing and deduplication...")
    for header, seq in parse_fasta(args.input):
        total_input_seqs += 1
        seqid = header.split()[0]
        tid = taxid_map.get(seqid)
        if not tid:
            print(f"Warning: No taxid found for sequence ID {seqid}")
            continue
        if tid not in taxid_to_seqs:
            taxid_to_seqs[tid] = {}
        if seq not in taxid_to_seqs[tid]:
            taxid_to_seqs[tid][seq] = (header, seq)
            total_retained_seqs += 1
        else:
            print(f"Duplicate sequence found for taxid {tid}, skipping.")

    print(f"Total input sequences: {total_input_seqs}")
    print(f"Total retained (unique) sequences: {total_retained_seqs}")
    print(f"Unique taxids found: {len(taxid_to_seqs)}")

    taxid_to_name = ncbi.get_taxid_translator([int(t) for t in taxid_to_seqs if t.isdigit()])
    print(f"Retrieved {len(taxid_to_name)} scientific names from NCBI.")

    species_codes = {}
    seen_codes = {}
    connector = "|"
    total_written = 0

    print(f"Writing deduplicated sequences to: {out_fasta}")
    with open(out_fasta, "w", encoding="utf-8") as fh_out:
        for tid, seqs in taxid_to_seqs.items():
            scname = taxid_to_name.get(int(tid), "") if tid.isdigit() else ""
            code_base = make_code(scname) if scname else (tid or "UNK")

            if scname:
                if scname in species_codes:
                    code = species_codes[scname]
                else:
                    code = code_base or "UNK"
                    if code in seen_codes and seen_codes[code] != scname:
                        suffix = 1
                        newcode = f"{code}{suffix}"
                        while newcode in seen_codes:
                            suffix += 1
                            newcode = f"{code}{suffix}"
                        code = newcode
                    species_codes[scname] = code
                    seen_codes[code] = scname
            else:
                code = code_base
                species_codes[scname] = code

            for seq, (header, sequence) in seqs.items():
                seqid = header.split()[0]
                gene = extract_gene(header)
                new_header = f"{code}{connector}{gene}{connector}{seqid}"
                fh_out.write(f">{new_header}\n{sequence}\n")
                total_written += 1

    print(f"Finished writing {total_written} sequences.")

    print(f"Writing species code map to: {out_codes}")
    with open(out_codes, "w", encoding="utf-8") as fh_codes:
        fh_codes.write("Scientific name\tCode\n")
        for name, code in species_codes.items():
            fh_codes.write(f"{name}\t{code}\n")

    print("All tasks completed successfully.")

if __name__ == "__main__":
    main()