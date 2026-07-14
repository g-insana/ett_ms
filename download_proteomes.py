#!/usr/bin/env python3

import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from tqdm import tqdm


def read_proteomes_map(tsv_path):
    proteome_map = {}
    with open(tsv_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2:
                continue
            proteome_id, upid = row[0].strip(), row[1].strip()
            if proteome_id == "proteome_id" and upid == "upid":
                continue
            proteome_map[proteome_id] = upid
    return proteome_map


def download_gz_file(url, out_path):
    gz_path = str(out_path) + ".gz"
    cmd = ["wget", "-O", gz_path, url]
    subprocess.run(cmd, check=True)

    with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)

    os.remove(gz_path)


def rewrite_fasta_headers(input_path, output_path):
    def parse_header(header_line):
        header = header_line[1:].strip()
        upi = header.split()[0]

        m = re.search(r"\bSS=([^ ]+)", header)
        sourceid = ""
        if m:
            ss_value = m.group(1)
            sourceid = ss_value.split(":", 1)[-1]

        return upi, sourceid

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        seq_header = None
        seq_chunks = []

        def flush_record():
            nonlocal seq_header, seq_chunks
            if seq_header is None:
                return
            upi, sourceid = parse_header(seq_header)
            sequence = "".join(seq_chunks).replace("\n", "").replace(" ", "")
            seqlen = len(sequence)
            fout.write(f">{sourceid}|{seqlen}|{upi}\n")
            for i in range(0, len(sequence), 60):
                fout.write(sequence[i:i+60] + "\n")
            seq_header = None
            seq_chunks = []

        for line in fin:
            if line.startswith(">"):
                flush_record()
                seq_header = line
            else:
                seq_chunks.append(line.strip())

        flush_record()


def process_proteome(proteome_id, upid, base_dir):
    out_fa = base_dir / f"proteome_{proteome_id}.fa"
    if out_fa.exists() and out_fa.stat().st_size > 0:
        return

    url = f"https://rest.uniprot.org/uniparc/proteome/{upid}/stream?compressed=true&format=fasta"

    tmp_fa = base_dir / f"proteome_{proteome_id}.fa.tmp"
    download_gz_file(url, tmp_fa)
    rewrite_fasta_headers(tmp_fa, out_fa)
    tmp_fa.unlink(missing_ok=True)


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} OSCODE", file=sys.stderr)
        sys.exit(1)

    oscode = sys.argv[1]
    base_dir = Path("data") / oscode
    proteomes_dir = base_dir / "proteomes"
    proteomes_dir.mkdir(parents=True, exist_ok=True)

    tsv_path = base_dir / f"{oscode}.proteomes.tsv"
    if tsv_path.exists() and tsv_path.stat().st_size > 0:
        proteome_map = read_proteomes_map(tsv_path)
    else:
        print(f"ERROR: data for species {oscode} is missing. Have you downloaded the input data?")
        sys.exit(2)

    for proteome_id, upid in tqdm(proteome_map.items(), desc="Proteomes"):
        process_proteome(proteome_id, upid, proteomes_dir)


if __name__ == "__main__":
    main()
