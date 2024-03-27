from pymoo.core.sampling import Sampling
import numpy as np

class SNPSampling(Sampling):
	"""
	Creates sample data from file. It chooses n_samples from snp.data (matrix
	of individuals) selecting randomly sol_size indices from snp.SNP_names. 
	They must be different and incremental.
	"""
	

	def _do(self, problem, n_samples, **kwargs):
		
		# Chosen start index
		index = 0
		# Store indices in array
		chosen_indices = []
		# Store samples in array
		samples = np.full((n_samples, 1), None, dtype=object)

		for i in range(n_samples): 
			for _ in range(problem.sol_size):
				snp_index = np.random.randint(index, problem.loci_size)
				chosen_indices.append(snp_index)
				index = snp_index + 1
			samples[i, 0] = chosen_indices
			index = 0
			 
		return samples
		
