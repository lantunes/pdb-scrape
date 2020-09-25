import gzip
import os
import warnings
from ftplib import FTP

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.StructureBuilder import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)
warnings.simplefilter("ignore", UserWarning)


if __name__ == '__main__':
    ftp = FTP("ftp.wwpdb.org")
    ftp.login()
    ftp.cwd("/pub/pdb/data/structures/all/pdb")

    filenames = []
    ftp.retrlines('NLST', callback=lambda line: filenames.append(line))
    print("files: %s" % len(filenames))

    for filename in filenames:
        print("downloading %s ..." % filename)
        with open(filename, 'wb') as fp:
            ftp.retrbinary("RETR %s" % filename, callback=fp.write)

        print("processing: %s" % filename)

        p = PDBParser()
        with gzip.open(filename, 'rt') as f:
            structure = p.get_structure("", f)

        pdb_id = structure.header["idcode"]

        assert pdb_id, "no PDB ID for %s" % filename

        model = structure[0]

        try:
            dssp = DSSP(model, filename, dssp="/Users/luis/dssp-2.3.0/mkdssp")
        except Exception as e:
            print(e)
            print()
            os.remove(filename)
            continue

        valid_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        current_chain = ""
        chain = ""
        phi_psis = []
        dihedrals = []
        chains = []
        chain_names = []
        # we assume that the chains and residues are in order, i.e. A1,A2,A3,...,B1,B2,B3,...
        for key in dssp.keys():

            chain_id = key[0]

            if not current_chain or chain_id != current_chain:
                current_chain = chain_id
                chain_names.append(chain_id)
                if chain:
                    assert len(chain) == len(phi_psis), \
                        "the length of chain '%s' does not equal the number of dihedrals: %s" % (chain, len(phi_psis))
                    chains.append(chain)
                    dihedrals.append(phi_psis)
                chain = ""
                phi_psis = []

            residue = dssp[key]
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

        with open("metadata-2020-09-25.txt", "a") as metadata_file:
            for line in zip([filename]*len(chain_names), [pdb_id]*len(chain_names), chain_names):
                metadata_file.write(",".join(line) + "\n")

        # dataset.extend(list(zip(chains, dihedrals)))

        os.remove(filename)

    ftp.quit()

    # TODO create a file that contains the filename, pdb_id and chain name of the datapoint in the dataset file
    #  first line of this file will correspond to the first line of the dataset file, etc.
    # store the date of this data scrape in the file name of the generated files
