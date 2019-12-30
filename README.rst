nocharge
========

A simple Python module that attempts to neutralize all atoms with +1 or -1 charge in a molecule.

For more information, see https://baoilleach.blogspot.com/2019/12/no-charge-simple-approach-to.html.

Usage:

.. code-block:: python

	from openbabel import pybel
	from nocharge import neutralize

	for smi in ["CC(=O)[O-]", "C[N+](C)(C)C"]:
	    mol = pybel.readstring("smi", smi).OBMol
	    altered = neutralize(mol)
	    if altered:
		outsmi = pybel.Molecule(mol).write("smi", opt={"n": True, "nonewline": True})
		print("{} changed to {}".format(smi, outsmi))
	    else:
		print("{} was unaltered by neutralize".format(smi))

which gives::

        CC(=O)[O-] changed to CC(=O)O
        C[N+](C)(C)C was unaltered by neutralize

