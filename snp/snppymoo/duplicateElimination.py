from pymoo.core.duplicate import ElementwiseDuplicateElimination
import numpy as np

class SNPDuplicateElimination(ElementwiseDuplicateElimination):
	
	"""
	Compares two SNPs group and retrives true if equals, false otherwise.
	"""

	def is_equal(self, a, b):
		return np.array_equal(a.X[0], b.X[0])
