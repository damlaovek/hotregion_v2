from Bio import BiopythonWarning
from Bio.PDB import *
import warnings
import sys, os
from interface_extraction.extract import InterfaceExtroctor
#from biopandas.pdb import PandasPdb
from clustering.dbscan import DoubleDBSCAN
import pandas as pd
from clustering.apscan import run_apscan
import urllib.request  
from helpers import write2file, get_com
from feature_calculation.asa import ASACalculator
from feature_calculation.pp import PPCalculator
from sklearn import metrics

def parse_pdb(structure, pdb_id, chain_id_1, chain_id_2):
    os.system("chmod -R 777 PDBs")
    
    class Chain1Select(Select):
        def accept_chain(self, chain):
            if chain.id==chain_id_1:
                return 1
            else:
                return 0
                
    class Chain2Select(Select):
        def accept_chain(self, chain):
            if chain.id==chain_id_2:
                return 1
            else:
                return 0
	
    class Chain12Select(Select):
        def accept_chain(self, chain):
            if (chain.id==chain_id_1) or (chain.id==chain_id_2):
                return 1
            else:
                return 0
	
    io = PDBIO()
    io.set_structure(structure)
    io.save("PDBs/%s%s.pdb" %(pdb_id.lower(), chain_id_1.lower()), Chain1Select())
    io.save("PDBs/%s%s.pdb" %(pdb_id.lower(), chain_id_2.lower()), Chain2Select())
    io.save("PDBs/%s%s.pdb" %(pdb_id.lower(), chain_id_1.lower()+chain_id_2.lower()), Chain12Select())
    return
        
def hotregion(interface_id, counter):
    '''
        HotRegion v.2: A new method to predict hot regions in protein-protein interfaces
        Usage: python main.py <interface id>
        
        Parameters:
        -----------
        interfaceID: String
                     4 letter pdb code + chain 1 identifier + chain 2 identifier.

        Returns:
        -------
        hotRegions: CSV file
                    Output file showing which residues belong to which regions and which residues are identified as noise.
    '''
    pdb_id = interface_id[0:4].upper()
    chain_id_1 = interface_id[4].upper()
    chain_id_2 = interface_id[5].upper()
    
    path = "PDBs/%s.pdb" %(pdb_id)
    url = "http://files.rcsb.org/download/%s.pdb" %(pdb_id)
        
    if os.path.exists(path) == False:
        urllib.request.urlretrieve(url, path)

	# Set and parse protein structure
    parser = PDBParser()
    structure = parser.get_structure('myPDBStructure', "PDBs/%s.pdb" %(pdb_id))
    parse_pdb(structure, pdb_id, chain_id_1, chain_id_2)
    structure = parser.get_structure('myPDBStructure', "PDBs/%s%s.pdb" %(pdb_id.lower(), chain_id_1.lower()+chain_id_2.lower()))
    model = structure[0]
    chain_1 = model[chain_id_1]
    chain_2 = model[chain_id_2]
    residue_list = Selection.unfold_entities(structure, 'R')  
    
    # Interface Extraction
    interface_extractor = InterfaceExtroctor(structure, chain_1, chain_2)
    interface_extractor.extract_interface()
    interface_extractor.non_interface = set(residue_list) - interface_extractor.interface_list
    #interface_extractor.extract_neighbors()
    #interface_extractor.interface_list.update(interface_extractor.neighbor_list)

    # Feature Calculation
    asa_calculator = ASACalculator("PDBs", interface_id, chain_1, chain_2)
    asa_calculator.run_naccess()

    pp_calculator = PPCalculator(residue_list)
    pp_calculator.calculate_pp()

    # DataFrame Generation
    data_dict = {'X':[],'Y':[], 'Z':[]}
    for residue in list(interface_extractor.interface_list):
        try:
            com = get_com(residue)
            data_dict['X'].append(com[0])
            data_dict['Y'].append(com[1])
            data_dict['Z'].append(com[2])
            #data_dict['relCompASA'].append(asa_calculator.relative_complex_asa[residue])
            #data_dict['relMonASA'].append(asa_calculator.relative_monomer_asa[residue])
            #data_dict['compASA'].append(asa_calculator.monomer_asa[residue])
            #data_dict['monASA'].append(asa_calculator.monomer_asa[residue])
            #data_dict['PP'].append(pp_calculator.pair_potentials[residue])
        except:
            interface_extractor.interface_list.remove(residue)

    df = pd.DataFrame(data_dict)
        
    # Clustering
    '''
    db = DoubleDBSCAN(5, 3)
    labels = db.cluster(df)
    i = 0
    for residue in interface_extractor.interface_list:
       print(residue.get_resname() + " " + str(residue.id[1]) + " " + str(labels[i]))
       i += 1
    '''

    labels = run_apscan(df)
    silhouette = metrics.silhouette_score(df, labels, metric='euclidean')
    write2file(interface_id, 
        interface_extractor.interface_list,
        labels, 
        asa_calculator.complex_asa, 
        asa_calculator.relative_complex_asa,
        asa_calculator.monomer_asa, 
        asa_calculator.relative_monomer_asa, 
        pp_calculator.pair_potentials
    )
    print(f"Silhouette score of {interface_id}: {silhouette}")
    print(str(counter)+" complex run sucessfull")

if __name__ == '__main__':
    warnings.simplefilter('ignore', BiopythonWarning)
    '''
    interface_df = pd.read_csv('interfaces.csv', sep=";")
    interface_df = interface_df.dropna(axis = 0, how ='any')
    counter = 1
    for interface_id in interface_df.iloc[:, 0]:
        try:
            hotregion(interface_id, counter)
            counter += 1
        except:
            pass  
    
    if len(sys.argv) == 2:
        interface_id = sys.argv[1].upper()
        hotregion(interface_id, 1)
               
    else:
        print("Usage: python main.py <interface id>")
    '''
    #experimental_list = ["1JTGAB","1A0OAB","1A22AB","1BRSAD","1EMVAB","1IBRAB","2C5LAC","2WPTAB","4G0NAB"]
    experimental_list = ["4ZQKAB"]
    counter = 1
    for interface_id in experimental_list:
        hotregion(interface_id, counter)
        counter += 1  
