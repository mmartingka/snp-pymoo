from pymoo.core.repair import Repair

class SNPRepair(Repair):
	
	"""
	It defines the repair method to ensure there are not repetaed SNPs groups
	"""

	def _do(self, problem, Z, **kwargs):

		"""
		Creates two new children from a given pair of parents.
		"""
		
		# Now repair each indvidiual i
		for i in range(len(Z)):
        	# The packing plan for i
			z = Z[i,0]
			for j in range(problem.dim_epi):
				while(self.buscar_snp(problem.dim_epi, z[j], z, j)):
					z[j]=(z[j]+1)%problem.loci_size	
				z.sort()
			Z[i,0] = z
		return Z
	
	def buscar_snp(self, dim_epi, snp_value, snp_sol, pos):
		
		"""
		Check for repeated SNPs.
		"""
		for i in range(dim_epi):
			if snp_value == snp_sol[i] and pos != i:
				return True
		return False