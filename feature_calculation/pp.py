from dictionaries import CONTACT_POTENTIALS, PP_indices
from helpers import get_com, distance_calculation

class PPCalculator:

	def __init__(self, residue_list):
		self.residue_list = residue_list
		self.pair_potentials = dict()


	def extract_neighbors(self, residue):
		com_residue = get_com(residue)
		neighbors = []
		for neighbor in self.residue_list:
			try:
				com_neighbor = get_com(neighbor)
				if (residue != neighbor):
					distance = distance_calculation(com_residue,com_neighbor)
					if (residue.get_parent().id == neighbor.get_parent().id): #they are in the same chain
						if(abs(residue.get_id()[1]-neighbor.get_id()[1]) >= 4):
							if (distance <= 7.0):
								neighbors.append(neighbor)
					else:
						if (distance <= 7.0):
							neighbors.append(neighbor)
			except:
				pass
		return neighbors
		
	def contact_potentials(self, residue):
		neighbors = self.extract_neighbors(residue)
		sum_pp = 0.0
		for neighbor in neighbors:
			indexI = PP_indices[residue.get_resname()]
			indexJ = PP_indices[neighbor.get_resname()]
			if(CONTACT_POTENTIALS[indexI][indexJ] != 0):
				sum_pp+=CONTACT_POTENTIALS[indexI][indexJ]
			else:
				sum_pp+=CONTACT_POTENTIALS[indexJ][indexI]
		return sum_pp

	def calculate_pp(self):
		for residue in self.residue_list:
			try:
				pp = self.contact_potentials(residue)
				self.pair_potentials[residue] = abs(pp)
			except:
				pass