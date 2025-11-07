
#!/usr/bin/env python3
import argparse
import os

# Set ETE3PATH before importing ete3
ete3_path = os.path.join(os.getcwd(), ".etetoolkit")
os.makedirs(ete3_path, exist_ok=True)
os.environ["ETE3PATH"] = ete3_path

from ete3 import NCBITax


# ----- Helpers -----
def load_taxid_map(path):
    """Load seqid -> taxid mapping from a tab-delimited file."""
    with open(path, "r", encoding="utf-8") as fh:
        return {parts[0]: parts[1] for parts in (line.strip().split("\t") for line in fh)}

def parse_fasta(path):
    """Simple FASTA parser: yields (header, sequence)."""
    header, seq = None, []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def extract_gene(desc):
    """Extract gene name from FASTA description."""
    parts = desc.split("[gene=", 1)
    return parts[1].split("]", 1)[0] if len(parts) > 1 else "unknowngene"

def make_code(scname):
    """Generate species code from scientific name."""
    if not scname:
        return ""
    toks = scname.strip().split()
    n = len(toks)
    lower = [t.lower() for t in toks]
    if n > 3 and 'x' in lower:  # hybrid case
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

def fetch_taxid_to_name_ete3(taxids):
    """Fetch scientific names for given taxids using ETE3 NCBITaxa."""
    ncbi = NCBITaxa()
    # Ensure taxonomy DB is present
    try:
        ncbi.get_taxid_translator([1])  # test query
    except Exception:
        print("Taxonomy database missing. Downloading...")
        ncbi.update_taxonomy_database()
    return ncbi.get_taxid_translator(taxids)

# ----- Main -----
def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers using taxid mapping and species codes (ETE3).")
    parser.add_argument("--taxidmap", required=True, help="Tab-delimited file: seqid<TAB>taxid")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--prefix", required=True, help="Output prefix for renamed FASTA and species codes")
    args = parser.parse_args()

    out_fasta = f"{args.prefix}_renamed.fasta"
    out_codes = f"{args.prefix}_speciescodes.txt"

    # Load taxid map
    taxid_map = load_taxid_map(args.taxidmap)

    # Collect taxids from FASTA
    taxids_needed = set()
    seqid_to_taxid = {}
    fasta_records = []
    for header, seq in parse_fasta(args.input):
        seqid = header.split()[0]
        tid = taxid_map.get(seqid)
        seqid_to_taxid[seqid] = tid or ""
        if tid:
            taxids_needed.add(tid)
        fasta_records.append((header, seq))

    # Fetch taxid -> scientific name using ETE3
    taxid_to_name = fetch_taxid_to_name_ete3([int(t) for t in taxids_needed if t.isdigit()])

    # Build species codes
    species_codes = {}
    seen_codes = {}
    connector = "|"

    with open(out_fasta, "w", encoding="utf-8") as fh_out:
        for header, seq in fasta_records:
            seqid = header.split()[0]
            gene = extract_gene(header)
            tid = seqid_to_taxid.get(seqid, "")
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

            new_header = f"{code}{connector}{gene}{connector}{seqid}"
            fh_out.write(f">{new_header}\n{seq}\n")

    # Write species codes file
    with open(out_codes, "w", encoding="utf-8") as fh_codes:
        fh_codes.write("Scientific name\tCode\n")
        for name, code in species_codes.items():
            fh_codes.write(f"{name}\t{code}\n")

    print(f"Sequences renamed successfully. Output saved to {out_fasta}.")
    print(f"Map to species scientific names and species code saved to {out_codes}.")

if __name__ == "__main__":
    main()