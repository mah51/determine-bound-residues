from Bio.PDB import NeighborSearch, PDBParser, Selection

pdb_name = '3vg9' # PDB file for analysis
pdb_file = './data/{}.pdb'.format(pdb_name)
CONTACT_DISTANCE = 4.5  # Contact distance between atoms in Angstroms
EXTRA_RESIDUES = 2 # Num of extra residues on either side of CDR
h_strand_label = 'C'
l_strand_label = 'B'
ag_strand_label = 'A'
chothia_cdr_def = {"L1": (24, 34), "L2": (50, 56), "L3": (89, 97),
                   "H1": (26, 32), "H2": (52, 56), "H3": (95, 102)}

parser = PDBParser()

# read structure from file
structure = parser.get_structure(pdb_name, pdb_file)

model = structure[0]
ab_h_chain = model[h_strand_label]
ab_l_chain = model[l_strand_label]
ag_chain = model[ag_strand_label]

h_residues = ab_h_chain.get_unpacked_list()
l_residues = ab_l_chain.get_unpacked_list()
# Get CDRS


def filter_cdrs(start, end, seq):
    cdr_list = []
    for residue in seq:
        if residue.id[1] >= start and residue.id[1] <= end:
            cdr_list.append(residue)
    return cdr_list

cdrs = {}
for chotia_label in chothia_cdr_def:
    sequence = {}
    if chotia_label.startswith('L'):
        sequence = l_residues
    elif chotia_label.startswith('H'):
        sequence = h_residues
    cdrs[chotia_label] = filter_cdrs(chothia_cdr_def[chotia_label][0] - EXTRA_RESIDUES,
                                     chothia_cdr_def[chotia_label][1] + EXTRA_RESIDUES, sequence)

# Search antigen
ag_search = NeighborSearch(Selection.unfold_entities(ag_chain, 'A'))


def fill_interactions(cdr):
    list_of_res = []
    for residue in cdr:
        for atom in residue.child_list:
            result = ag_search.search(atom.coord, CONTACT_DISTANCE)
            for result_atom in result:
                if str(residue.id[1]) + residue.id[2] + residue.resname not in list_of_res:
                    list_of_res.append(str(residue.id[1]) + residue.id[2] + residue.resname)
    return list_of_res


cdr_interactions = {}

for cdr in cdrs:
    cdr_interactions[cdr] = fill_interactions(cdrs[cdr])
    print(cdr)
    print(cdr_interactions[cdr])