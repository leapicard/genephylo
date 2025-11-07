#!/usr/bin/env python
import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO, Entrez

# ----- Helpers -----
def load_taxid_map(path):
    """Load seqid -> taxid mapping from a tab-delimited file."""
    with open(path, "r", encoding="utf-8") as fh:
        return {parts[0]: parts[1] for parts in (line.strip().split("\t") for line in fh)}

def chunked(iterable, n):
    it = iter(iterable)
    while True:
        chunk = []
        try:
            for _ in range(n):
                chunk.append(next(it))
        except StopIteration:
            if chunk:
                yield chunk
            break
        yield chunk

def fetch_taxid_to_name(taxids):
    """Fetch scientific names for given taxids using NCBI Entrez."""
    ids = ",".join(taxids)
    handle = Entrez.esummary(db="taxonomy", id=ids)
    records = Entrez.read(handle)
    return {str(rec.get("Id")): rec.get("ScientificName") or rec.get("TaxId") or "" for rec in records}

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
    parser = argparse.ArgumentParser(description="Rename FASTA headers using taxid mapping and species codes.")
    parser.add_argument("--taxidmap", required=True, help="Tab-delimited file: seqid<TAB>taxid")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--prefix", required=True, help="Output prefix for renamed FASTA and species codes")
    parser.add_argument("--email", help="Email for NCBI Entrez (recommended)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for Entrez queries")
    args = parser.parse_args()

    # Configure Entrez
    Entrez.email = args.email or ""

    out_fasta = f"{args.prefix}_renamed.fasta"
    out_codes = f"{args.prefix}_speciescodes.txt"

    # Load taxid map
    taxid_map = load_taxid_map(args.taxidmap)

    # Collect taxids from FASTA
    taxids_needed = set()
    seqid_to_taxid = {}
    for rec in SeqIO.parse(args.input, "fasta"):
        tid = taxid_map.get(rec.id)
        seqid_to_taxid[rec.id] = tid or ""
        if tid:
            taxids_needed.add(tid)

    # Fetch taxid -> scientific name in parallel
    taxid_to_name = {}
    chunks = list(chunked(list(taxids_needed), 50))  # 50 per request
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(fetch_taxid_to_name, chunk): chunk for chunk in chunks}
        for future in as_completed(futures):
            try:
                taxid_to_name.update(future.result())
            except Exception as e:
                print(f"Error fetching chunk {futures[future]}: {e}")
            time.sleep(0.34)  # respect NCBI rate limits globally

    # Build species codes
    species_codes = {}
    seen_codes = {}
    connector = "|"

    with open(out_fasta, "w", encoding="utf-8") as fh_out:
        for rec in SeqIO.parse(args.input, "fasta"):
            seqid = rec.id
            gene = extract_gene(rec.description)
            tid = seqid_to_taxid.get(seqid, "")
            scname = taxid_to_name.get(tid, "")
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

            header = f"{code}{connector}{gene}{connector}{seqid}"
            fh_out.write(f">{header}\n{rec.seq}\n")

    # Write species codes file
    with open(out_codes, "w", encoding="utf-8") as fh_codes:
        fh_codes.write("Scientific name\tCode\n")
        for name, code in species_codes.items():
            fh_codes.write(f"{name}\t{code}\n")

    print(f"Sequences renamed successfully. Output saved to {out_fasta}.")
    print(f"Map to species scientifique names and species code saved to {out_codes}.")

if __name__ == "__main__":
    main()