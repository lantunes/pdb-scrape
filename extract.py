import gzip
import warnings
from os import listdir
from os.path import isfile, join

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.StructureBuilder import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)
warnings.simplefilter("ignore", UserWarning)


if __name__ == '__main__':

    for file in listdir("resources"):
        filename = join("resources", file)
        if isfile(filename):

            print("processing: %s" % filename)

            p = PDBParser()
            with gzip.open(filename, 'rt') as f:
                structure = p.get_structure("", f)

            model = structure[0]

            try:
                dssp = DSSP(model, filename, dssp="/Users/luis/dssp-2.3.0/mkdssp")
            except Exception as e:
                print(e)
                print()
                continue

            valid_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
            current_chain = ""
            chain = ""
            phi_psis = []
            dihedrals = []
            chains = []
            # we assume that the chains and residues are in order, i.e. A1,A2,A3,...,B1,B2,B3,...
            for key in dssp.keys():

                chain_id = key[0]
                residue_id = key[1][1]

                if not current_chain or chain_id != current_chain:
                    current_chain = chain_id
                    if chain:
                        assert len(chain) == len(phi_psis), \
                            "the length of chain '%s' does not equal the number of dihedrals: %s" % (chain, len(phi_psis))
                        chains.append(chain)
                        dihedrals.append(phi_psis)
                    chain = ""
                    phi_psis = []

                residue = dssp[(chain_id, residue_id)]
                amino_acid = residue[1]
                assert amino_acid in valid_aa, "unsupported amino acid: %s" % amino_acid
                phi = residue[4]
                assert isinstance(phi, float), "phi is not a float: %s" % phi
                psi = residue[5]
                assert isinstance(psi, float), "psi is not a float: %s" % psi

                chain += amino_acid
                phi_psis.append((phi, psi))

            if chain:
                assert len(chain) == len(phi_psis), \
                    "the length of chain '%s' does not equal the number of dihedrals: %s" % (chain, len(phi_psis))
                chains.append(chain)
                dihedrals.append(phi_psis)

            print(chains)
            print(dihedrals)
            print()
