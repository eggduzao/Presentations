#!/usr/bin/env python3
"""
Take demo.sam, find the 3 "median" reads (by order in the file),
and replicate each N_PER_TEMPLATE times into spike.sam.
"""

N_PER_TEMPLATE = 250

hdr = []
alns = []

with open("demo.sam") as f:
    for line in f:
        if line.startswith('@'):
            hdr.append(line)
        else:
            alns.append(line.rstrip())

if not alns:
    raise ValueError("No alignments found in demo.sam!")

# Find 3 median reads
n = len(alns)
mid = n // 2
if n < 3:
    raise ValueError("Need at least 3 reads to pick medians.")
idxs = [max(0, mid - 1), mid, min(n - 1, mid + 1)]
templates = [alns[i] for i in idxs]

# Write header + duplicated median reads
with open("spike.sam", "w") as out:
    out.writelines(hdr)
    uid = 0
    for t in templates:
        fields = t.split('\t')
        for i in range(N_PER_TEMPLATE):
            uid += 1
            newfields = fields[:]  # shallow copy
            newfields[0] = f"{fields[0]}_dup{uid}"  # unique QNAME
            out.write('\t'.join(newfields) + '\n')