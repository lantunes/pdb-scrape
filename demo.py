from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import gzip


if __name__ == '__main__':

    p = PDBParser()
    with gzip.open("resources/pdb1a00.ent.gz", 'rt') as f:
        structure = p.get_structure("", f)

    model = structure[0]

    dssp = DSSP(model, "resources/pdb1a00.ent.gz", dssp="/Users/luis/dssp-2.3.0/mkdssp")

    # DSSP data is accessed by a tuple - (chain id, residue id)
    # a tuple representing the data for residue #1 of chain A
    residue_A_1 = dssp[("A", 1)]

    """
    Tuple Index     Value
    0               DSSP index
    1               Amino acid
    2               Secondary structure
    3               Relative ASA
    4               Phi
    5               Psi
    6               NH–>O_1_relidx
    7               NH–>O_1_energy
    8               O–>NH_1_relidx
    9               O–>NH_1_energy
    10              NH–>O_2_relidx
    11              NH–>O_2_energy
    12              O–>NH_2_relidx
    13              O–>NH_2_energy
    
    Secondary structure:
    Code    Structure
    H       Alpha helix (4-12)
    B       Isolated beta-bridge residue
    E       Strand
    G       3-10 helix
    I       Pi helix
    T       Turn
    S       Bend
    -       None
    """

    # Amino acid symbol
    print(residue_A_1[1])

    # Phi
    print(residue_A_1[4])

    # Psi
    print(residue_A_1[5])
