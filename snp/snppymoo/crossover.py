from pymoo.core.crossover import Crossover
import random
import numpy as np

class SNPCrossover(Crossover):
	"""
	It defines the crossover process (obtain two children from two parents)
	
	Attributes:
		prob: Probability of applying the crossover operator. If it's 1, it'll 
			always be applied. It's useful to have values lower than 1 to have
			generations without crossovers and thus maintain genetic diversity
	"""
	
	def __init__(self, prob=1.0):
	
		""" 
		Define the crossover: number of parents and number of offsprings
		"""
		super().__init__(2, 2)
		self.prob = prob
		
	def _do(self, problem, X, **kwargs):

		"""
		Creates two new children from a given pair of parents
		"""

		# The input has the following shape (n_parents, n_matings, n_var)
		_, n_matings, _ = X.shape

		# The output with the shape (n_offsprings, n_matings, n_var)
		# As the number of parents and offsprings are equal it keeps X's shape
		Y = np.full_like(X, None, dtype=object)

		dim_epi = problem.dim_epi

		for m in range(n_matings):
        	# Get the first and the second parent
			a, b = X[0, m, 0], X[1, m, 0]

			# prepare the offsprings
			off_a = np.full(dim_epi, None, dtype=object)
			off_b = np.full(dim_epi, None, dtype=object)

			if random.randint(0, 100) < self.prob:
				if dim_epi>2:
					pto1=random.randint(0, dim_epi-2)
					pto2=random.randint(pto1, dim_epi-1)
				else:
					pto1=0
					pto2=0
				for i in range(dim_epi):
					if i<pto1 or i>pto2:
						off_a[i]=a[i]
						off_b[i]=b[i]
					else:
						off_a[i]=b[i]
						off_b[i]=a[i]
			else:
				off_a = np.asarray(a)
				off_b = np.asarray(b)

			Y[0, m, 0] = np.sort(off_a)
			Y[1, m, 0] = np.sort(off_b)

		return Y
		
		
