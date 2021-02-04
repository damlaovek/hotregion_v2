import pandas as pd
from dictionaries import DICT_ATOM
import math

def three2one(residue):
    if residue ==   "ALA": return "A"
    elif residue == "CYS": return "C"
    elif residue == "ASP": return "D"
    elif residue == "GLU": return "E"
    elif residue == "PHE": return "F"
    elif residue == "GLY": return "G"
    elif residue == "HIS": return "H"
    elif residue == "ILE": return "I"
    elif residue == "LYS": return "K"
    elif residue == "LEU": return "L"
    elif residue == "MET": return "M"
    elif residue == "ASN": return "N"
    elif residue == "PRO": return "P"
    elif residue == "GLN": return "Q"
    elif residue == "ARG": return "R"
    elif residue == "SER": return "S"
    elif residue == "THR": return "T"
    elif residue == "VAL": return "V"
    elif residue == "TRP": return "W"
    elif residue == "TYR": return "Y"
    else: return "X"
    
def write2file(interfaceID, residue_list, labels, compASAs, relCompASAs, monASAs, relMonASAs, PPs):
    res_ids = []
    res_names = []
    chain_ids= []
    comp_ASAs = []
    rel_comp_ASAs = []
    mon_ASAs = []
    rel_mon_ASAs = []
    pair_potentials = []
    info = []
    for residue in residue_list:
        res_ids.append( str(residue.id[1]) )
        res_names.append( three2one(residue.get_resname()) )
        chain_ids.append( residue.get_parent().id )
        comp_asa = "%6.2f"%(compASAs[residue])
        comp_ASAs.append(comp_asa)
        rel_comp_asa = "%6.2f"%(relCompASAs[residue])
        rel_comp_ASAs.append(rel_comp_asa)
        mon_asa = "%6.2f"%(monASAs[residue])
        mon_ASAs.append(mon_asa)
        rel_mon_asa = "%6.2f"%(relMonASAs[residue])
        rel_mon_ASAs.append(rel_mon_asa)
        pp = "%6.2f"%(PPs[residue])
        pair_potentials.append(pp)
        if (float(relCompASAs[residue]) <=20 and float(PPs[residue])>=18.0):
            info.append('HS')
        else:
            info.append('NS')
    d = {'RES_ID': res_ids, 'RES_NAME': res_names, 'CHAIN': chain_ids, 'CLUSTER': labels, 'Complex ASA':comp_ASAs, 'Relative Complex ASA': rel_comp_ASAs, 'Monomer ASA':mon_ASAs, 'Relative Monomer ASA': rel_mon_ASAs, 'PP': pair_potentials, 'Status': info}
    df = pd.DataFrame(data=d)
    df.to_csv(interfaceID+'_data.csv', sep=',', encoding='utf-8')

def distance_calculation(vector_1, vector_2):
    x1 = vector_1[0]
    y1 = vector_1[1]
    z1 = vector_1[2]

    x2 = vector_2[0]
    y2 = vector_2[1]
    z2 = vector_2[2]

    distance = math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))
    return distance

def get_com(residue):
    x_total = 0
    y_total = 0
    z_total = 0
    if(len(residue.get_resname()) >3):
        resname = residue.get_resname()[1:4]
    else:
        resname = residue.get_resname()
    num_atoms = len(DICT_ATOM[resname])
    for atom_id in DICT_ATOM[resname]:
        atom = residue[atom_id]
        x_total += atom.get_coord()[0]
        y_total += atom.get_coord()[1]
        z_total += atom.get_coord()[2]
    return [x_total/num_atoms,y_total/num_atoms,z_total/num_atoms]
    
