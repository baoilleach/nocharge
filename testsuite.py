import unittest
from openbabel import pybel
ob = pybel.ob
from nocharge import neutralize

class Basic(unittest.TestCase):
    def testCOOH(self):
        mol = pybel.readstring("smi", "CC(=O)[O-]").OBMol
        neutralize(mol)
        self.assertEqual("CC(=O)O", pybel.Molecule(mol).write("smi").rstrip())
    def testAmino(self):
        mol = pybel.readstring("smi", "CC[NH3+]").OBMol
        neutralize(mol)
        self.assertEqual("CCN", pybel.Molecule(mol).write("smi").rstrip())
    def testLeaveUnchanged(self):
        data = [
                "CC[N+](=O)[O-]", # nitro
                "c1ccc[n+]([O-])c1", # pyridine N-oxide
                "C1=NC(=O)NC(=O)C1=[N+]=N", # diazenium
                ]
        for smi in data:
            mol = pybel.readstring("smi", smi).OBMol
            neutralize(mol)
            self.assertEqual(smi, pybel.Molecule(mol).write("smi").rstrip())

if __name__ == "__main__":
    unittest.main()
