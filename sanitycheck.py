"""Check that all four implementations give the same result.
"""
import os
import time
import gzip
from openbabel import pybel
import nocharge

fns = [nocharge.neutralize_v1, nocharge.neutralize_v2, nocharge.neutralize_v3, nocharge.neutralize_v4, nocharge.neutralize_v5]
times = [0] * 5

testfile = os.path.join("data", "charged.sorted.smi.gz")
with gzip.open(testfile, "rt") as inp:
    for N, line in enumerate(inp):
        oldcan = None
        for i, fn in enumerate(fns):
            mol = pybel.readstring("smi", line).OBMol
            t = time.time()
            fn(mol)
            times[i] += time.time() - t
            can = pybel.Molecule(mol).write("can")
            if oldcan and can != oldcan:
                assert False # should never happen
            oldcan = can
print(times)
