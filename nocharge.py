from openbabel import pybel
ob = pybel.ob

def NoNegativelyChargedNbr(atom):
    for nbr in ob.OBAtomAtomIter(atom):
        if nbr.GetFormalCharge() < 0:
            return False
    return True

def NoPositivelyChargedNbr(atom):
    for nbr in ob.OBAtomAtomIter(atom):
        if nbr.GetFormalCharge() > 0:
            return False
    return True

def neutralize_v1(mol):
    """Neutralize charges of +1 or -1"""
    for atom in ob.OBMolAtomIter(mol):
        chg = atom.GetFormalCharge()
        if chg == 1:
            hcount = atom.GetImplicitHCount()
            if hcount >= 1 and NoNegativelyChargedNbr(atom):
                atom.SetFormalCharge(0)
                atom.SetImplicitHCount(hcount - 1)
        elif chg == -1: # -ve charge
            hcount = atom.GetImplicitHCount()
            if NoPositivelyChargedNbr(atom): # avoid nitro, etc.
                atom.SetFormalCharge(0)
                atom.SetImplicitHCount(hcount + 1)

pat_v2 = ob.OBSmartsPattern()
pat_v2.Init("[+1,-1]")

def neutralize_v2(mol):
    """Neutralize charges of +1 or -1"""
    pat_v2.Match(mol)
    matches = pat_v2.GetMapList()
    for match in matches:
        atom = mol.GetAtom(match[0])
        chg = atom.GetFormalCharge()
        if chg == 1: # +1 charge
            hcount = atom.GetImplicitHCount()
            if hcount >= 1 and NoNegativelyChargedNbr(atom):
                atom.SetFormalCharge(0)
                atom.SetImplicitHCount(hcount - 1)
        else: # -1 charge
            hcount = atom.GetImplicitHCount()
            if NoPositivelyChargedNbr(atom):
                atom.SetFormalCharge(0)
                atom.SetImplicitHCount(hcount + 1)

pat_v3 = ob.OBSmartsPattern()
pat_v3.Init("[+1!H0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")

def neutralize_v3(mol):
    """Neutralize charges of +1 or -1"""
    pat_v3.Match(mol)
    for match in pat_v3.GetMapList():
        atom = mol.GetAtom(match[0])
        chg = atom.GetFormalCharge()
        hcount = atom.GetImplicitHCount()
        atom.SetFormalCharge(0)
        atom.SetImplicitHCount(hcount - chg)

matchvector_v4 = ob.vectorvInt()

def neutralize_v4(mol):
    """Neutralize charges of +1 or -1"""
    match_found = pat_v3.Match(mol, matchvector_v4)
    if match_found:
        for match in matchvector_v4:
            atom = mol.GetAtom(match[0])
            chg = atom.GetFormalCharge()
            hcount = atom.GetImplicitHCount()
            atom.SetFormalCharge(0)
            atom.SetImplicitHCount(hcount - chg)


neutralize_op = ob.OBOp.FindType("neutralize")
def neutralize_v5(mol):
    """Neutralize charges of +1 or -1"""
    neutralize_op.Do(mol)

# The following is the same as v4 but I've added a boolean
# return value to indicate whether the molecule was altered.
# This is useful to avoid regenerating a SMILES string for
# those cases.
def neutralize(mol):
    """Neutralize charges of +1 or -1"""
    match_found = pat_v3.Match(mol, matchvector_v4)
    if not match_found:
        return False
    for match in matchvector_v4:
        atom = mol.GetAtom(match[0])
        chg = atom.GetFormalCharge()
        hcount = atom.GetImplicitHCount()
        atom.SetFormalCharge(0)
        atom.SetImplicitHCount(hcount - chg)
    return True

if __name__ == "__main__":
    smi = "C(=O)[O-]"
    mol = pybel.readstring("smi", smi).OBMol
    neutralize(mol)
    outsmi = pybel.Molecule(mol).write("smi",
            opt={"n": True, "nonewline": True}) # avoid whitespace after SMILES
    assert outsmi == "C(=O)O"
    print("Input:  {}\nOutput: {}".format(smi, outsmi))

