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
		self.prob_mutation = prob_mutation # [0, 1] values
		self.range_mut = range_mut 

		
	def _do(self, problem, X, **kwargs):

		"""
		For each individual, it mutates only if a random number is below the
		selected probability mutation. If so, the process is done until a good
		mutation is performed (no repeted and order SNPs).  
		"""
		
		sol_size = problem.sol_size
		
		# For each individual
		for i in range(len(X)):
			for k in range(sol_size):
				if random.randint(101) < self.prob_mutation:
					X[i, 0] = self.mutate_snp(X[i, 0], problem, k)
					while self.buscar_snp(sol_size, X[i, 0][k], X[i, 0], k):
						X[i, 0] = self.mutate_snp(X[i, 0], problem, k)
					X[i,0].sort()
		return X
		
	def mutate_snp(self, snp, problem, k):
	
		"""
		Mutate a value between [-range_mut, range_mut] for the snp. If it 
		exceeds the boundaries of 0 or the size of loci, then the mutation 
		value will be these respective edges.
		"""
		
		mut = 0
		while mut == 0:
			mut = self.range_mut - np.random.randint(self.range_mut * 2 + 1)
		snp[k] += mut
		if snp[k] < 0:
			snp[k] = 0
		elif snp[k] > problem.loci_size - 1:
			snp[k] = problem.loci_size - 1
		return snp

	def buscar_snp(self, sol_size, snp_value, snp_sol, pos):
	
		"""
		Check that the mutation does not generate repeated SNPs. 
		"""
		
		for i in range(sol_size):
			if snp_value == snp_sol[i] and pos != i:
				return True
		return False
