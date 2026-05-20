#!/usr/bin/env python3
"""
Fetch FASTA sequences for given accessions.

Usage:
  python fetch_fasta.py 'AAN45060|502|UPI0000139077' 'AAN45060|502' [...]

Behaviour:
- Input is a single accession or a space-separated list of accessions
- Each accession is either:
    1) PID|SEQLEN|UPI (three fields) -> use UPI to fetch via uniparc rest api
    2) PID|SEQLEN (two fields) -> use PID to fetch from uniparc
- The returned FASTA header is replaced with the original accession string
- All results are written to stdout concatenated, in FASTA format
- Non-200 responses are printed to stderr and skipped
"""

import sys
import urllib.parse
import urllib.request
import ssl

API_UPI = "https://rest.uniprot.org/uniparc/{upi}.fasta"
API_PID = "https://rest.uniprot.org/uniparc/stream?format=fasta&query=(dbid:{pid})"

def fetch_url(url, timeout=30):
    ctx = ssl.create_default_context()
    req = urllib.request.Request(url, headers={"User-Agent": "down_fasta_acc/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        status = r.getcode()
        body = r.read().decode('utf-8')
        return status, body


def normalize_accession(acc_str):
    # keep the original string as the label to print in header
    parts = acc_str.split('|')
    return acc_str, parts


def build_url_from_parts(parts):
    if len(parts) >= 3 and parts[2].startswith("UPI"):
        upi = parts[2]
        return API_UPI.format(upi=urllib.parse.quote(upi, safe=''))
    else:
        # use first field as PID
        pid = parts[0]
        return API_PID.format(pid=urllib.parse.quote(pid, safe=''))


def replace_header(fasta_text, new_header_label):
    lines = fasta_text.splitlines()
    out_lines = []
    for line in lines:
        if line.startswith('>'):
            # replace header with the given accession
            out_lines.append(f'>{new_header_label}')
        else:
            out_lines.append(line)
    return "\n".join(out_lines).rstrip() + "\n"


def main(argv):
    if len(argv) < 2:
        print("Usage: python fetch_fasta.py <accession1> [accession2 ...]", file=sys.stderr)
        sys.exit(1)

    inputs = argv[1:]
    for acc in inputs:
        label, parts = normalize_accession(acc)
        url = build_url_from_parts(parts)
        try:
            status, body = fetch_url(url)
        except Exception as e:
            print(f"Error fetching {acc} from {url}: {e}", file=sys.stderr)
            continue

        if status != 200 or not body:
            print(f"Warning: non-200 response for {acc} (HTTP {status}) from {url}", file=sys.stderr)
            continue

        fasta = replace_header(body, label)
        sys.stdout.write(fasta)


if __name__ == "__main__":
    main(sys.argv)
