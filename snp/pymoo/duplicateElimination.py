from pymoo.core.duplicate import ElementwiseDuplicateElimination
import numpy as np

class SNPDuplicateElimination(ElementwiseDuplicateElimination):
	
	"""
	
	"""

	def is_equal(self, a, b):
		return np.array_equal(a.X[0], b.X[0])
