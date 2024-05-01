from pymoo.core.repair import Repair
import numpy as np

class SNPRepair(Repair):
	"""
	It defines the crossover process (obtain two children from two parents)
	
	Attributes:
		prob: Probability of applying the crossover operator. If it's 1, it'll 
			always be applied. It's useful to have values lower than 1 to have
			generations without crossovers and thus maintain genetic diversity
	"""

	def _do(self, problem, Z, **kwargs):

		"""
		Creates two new children from a given pair of parents
		"""
		
		# now repair each indvidiual i
		for i in range(len(Z)):
        	# the packing plan for i
			z = Z[i,0]
			for j in range(problem.dim_epi):
				while(self.buscar_snp(problem.dim_epi, z[j], z, j)):
					z[j]=(z[j]+1)%problem.loci_size	
				z.sort()
			Z[i,0] = z
		return Z
	
	def buscar_snp(self, dim_epi, snp_value, snp_sol, pos):
		
		"""
		Check that the mutation does not generate repeated SNPs. 
		"""
		for i in range(dim_epi):
			if snp_value == snp_sol[i] and pos != i:
				return True
		return False