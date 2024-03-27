from pymoo.core.crossover import Crossover
import random

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
		"""

		# The input has the following shape (n_parents, n_matings, n_var)
		_, n_matings, _ = X.shape

		# The output with the shape (n_offsprings, n_matings, n_var)
		# As the number of parents and offsprings are equal it keeps X's shape
		Y = np.full_like(X, None, dtype=object)

		sol_size = problem.sol_size

		for m in range(n_matings):
        	# Get the first and the second parent
			a, b = X[0, m, 0], X[1, m, 0]

			# prepare the offsprings
			off_a = np.full(problem.sol_size, None, dtype=object)
			off_b = np.full(problem.sol_size, None, dtype=object)

			if sol_size>2:
				pto1=random.randint(0, sol_size-2)
				pto2=random.randint(pto1, sol_size-1)
			else:
				pto1=0
				pto2=0

			for i in range(sol_size):
				if i<pto1 or i>pto2:
					off_a[i]=a[i]
					off_b[i]=b[i]
				else:
					off_a[i]=b[i]
					off_b[i]=a[i]

			Y[0, m, 0] = off_a.sort()
			Y[1, m, 0] = off_b.sort()
		return Y
		
		
