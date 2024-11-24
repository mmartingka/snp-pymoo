from pymoo.core.mutation import Mutation
import random
import numpy as np

class SNPMutation(Mutation):
	"""
	It defines the mutation process
	
	Attributes:
		prob_mutation: Probability of the SNP to be mutated. If it's close to 
		1 the probability is higher, otherwise, if is close to 0, it's lower.
		range_mut: Range within the mutation of a SNP varies. If it's 10, it 
		ranges between [-10, 10]
	"""
	
	def __init__(self, prob_mutation, range_mut):
	
		""" 
		Define the mutation
		"""
		super().__init__()
		self.prob_mutation = prob_mutation # [0-100] values
		self.range_mut = range_mut 

		
	def _do(self, problem, X, **kwargs):

		"""
		For each individual, it mutates only if a random number is below the
		selected probability mutation. If so, the process is done until a good
		mutation is performed (no repeted and order SNPs).  
		"""
		for i in range(len(X)):  # Iterate through each individual
			for k in range(problem.dim_epi):  # Iterate through each SNP
				if np.random.rand() * 100 < self.prob_mutation:  # Apply mutation with probability
					X[i][k] = self.mutate_snp(X[i][k], problem)
		return X

	def mutate_snp(self, snp_value, problem):
		"""
		Mutate a value between [-range_mut, range_mut] for the snp. If it 
		exceeds the boundaries of 0 or the size of loci, then the mutation 
		value will be these respective edges.        """
		mut = 0
		while mut == 0:  # Ensure mutation is non-zero
			mut = np.random.randint(-self.range_mut, self.range_mut)

		new_value = snp_value + mut
		# Ensure the value stays within the valid range
		new_value = max(0, min(new_value, problem.loci_size - 1))
		return new_value
		