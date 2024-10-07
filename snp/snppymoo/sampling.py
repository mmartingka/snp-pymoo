from pymoo.core.sampling import Sampling
import numpy as np

class SNPSampling(Sampling):
	"""
	Creates sample data from file. It chooses n_samples from snp.data (matrix
	of individuals) selecting randomly a number of dim_epi indices from 
	snp.SNP_names. They must be different and incremental.
	"""
	

	def _do(self, problem, n_samples, **kwargs):

		# Store samples in array
		samples = np.full((n_samples, 1), None, dtype=object)

		for i in range(n_samples): 
			# Store indices in array
			chosen_indices = []
			
			# Chosen start index
			index_init = 0
			
			# Chosen end index
			index_end = problem.loci_size - problem.dim_epi
				
			for _ in range(problem.dim_epi):
				snp_index = np.random.randint(index_init, index_end)
				chosen_indices.append(snp_index)
				
				index_init = snp_index + 1
				index_end = index_end + 1	
			
			indices = np.array(chosen_indices, dtype=object)
			samples[i, 0] = indices
		return samples
		
