import pandas as pd
import dictionaries

class InterfaceExtroctor:
    def __init__(self, structure, chain_1, chain_2):
        self.structure = structure
        self.chain_1 = chain_1
        self.chain_2 = chain_2
        self.interface_list = set()
        self.non_interface = set()
        self.neighbor_list = set()
        
    def is_atoms_interacting(self, atom_1, atom_2):
        r1 = dictionaries.VDW_RADII_EXTENDED[atom_1.get_parent().get_resname()][atom_1.get_id()]
        r2 = dictionaries.VDW_RADII_EXTENDED[atom_2.get_parent().get_resname()][atom_2.get_id()]
        distance = abs(atom_1 - atom_2)
        return True if distance < (r1+r2+0.5) else False
            
    
    def is_residues_interacting(self, residue_1, residue_2):
        for atom_1 in residue_1:
            for atom_2 in residue_2:
                try:
                    if self.is_atoms_interacting(atom_1, atom_2):
                        return True
                except:
                    pass
        return False
    
    def is_interface_residue(self, residue_1):
        for residue_2 in self.chain_2:
            if self.is_residues_interacting(residue_1, residue_2):
                self.interface_list.add(residue_1)
                self.interface_list.add(residue_2)
        return
    
    def extract_interface(self):
        for residue_1 in self.chain_1:
            self.is_interface_residue(residue_1)

    def find_neighbors(self, residue):
        for neighbor in self.non_interface:
            if neighbor.get_resname() != "HOH":
                if neighbor.get_parent() == residue.get_parent():
                    try:
                        if (abs(residue['CA'] - neighbor['CA']) <= 6) and (abs(residue.get_id()[1]-neighbor.get_id()[1]) >= 4):
                            self.neighbor_list.add(neighbor)
                    except:
                        pass
                else:
                	try:
                		if abs(residue['CA'] - neighbor['CA']) <= 6:
                			self.neighbor_list.add(neighbor)
                	except:
                		pass
        
    def extract_neighbors(self):
        for residue in self.interface_list:
            self.find_neighbors(residue)

