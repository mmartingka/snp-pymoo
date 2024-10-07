import numpy as np
import math
from pymoo.core.problem import ElementwiseProblem

class SNPProblem(ElementwiseProblem):
	
	"""
	Define the objective functions. 
	
	Attributes:
		dim_epi: Number of SNP to be chosen
		loci_size: Number of loci (SNP) of each individual
		sample_size: Number of individuals in the sample 
		data: Values for each SNP on each individual 
	"""
	
	def __init__(self, **kwargs):
		super().__init__(n_var=1, n_obj=2, n_ieq_constr=0,  
				   elementwise_evaluation=False, **kwargs)

		values = list(kwargs.values())
		self.dim_epi = values[1]
		self.loci_size = values[2]
		self.sample_size = values[3]
		self.data = values[4]

		
	def _evaluate(self, x, out, *args, **kwargs):
		score = self.bayesian_score(x[0])
		aic = self.logistic_score(x[0])
		out["F"] = np.array([+ score, + aic], dtype=float)

		
	def bayesian_score(self, x):
	
		"""
		Calculates the measurement of the relationship between SNP nodes and 
		disease nodes
		"""
		
		comb = int(math.pow(3.0, self.dim_epi))
		observed_values = np.zeros((2, comb))
		col_sum_table = np.zeros(comb)
		score = 0.0

		for i in range(self.sample_size):
			index = 0
			cont = True
			for j in range(self.dim_epi):
				if self.data[i][x[j]] == 3:
					cont = False
					break
				else:
					power_value = int(math.pow(3.0, (self.dim_epi-1-j)))
					index = index + self.data[i][x[j]] * power_value
			
			if cont:
				observed_values[self.data[i][self.loci_size]][index] += 1
				col_sum_table[index] += 1
		
		for i in range(comb):
			score = (score + self.snp_factorial(observed_values[0][i])
			 + self.snp_factorial(observed_values[1][i])
			   - self.snp_factorial(col_sum_table[i] + 1))

		return abs(score)
		
	def logistic_score(self, x):
	
		"""
		Calculates the likelihood of the model returning higher relationship
		with the disease as lower the value of the logistic score is.
		"""

		delta = 0.001
		maxiter = 30
		aic = 0
		testsample = self.sample_size
		theta_size = self.dim_epi + 2
		newdata = np.ones((testsample, theta_size))
	
		# Update middle columns of newdata
		for h in range(1, theta_size - 1):
			newdata[:, h] = self.data[:, x[h - 1]]
		
		# Update last column of newdata  by multiplying values of middle columns
		newdata[:, theta_size - 1] = np.prod(newdata[:, 1:theta_size - 1], axis=1)

		theta = np.zeros(theta_size)
		n_iter = 0
		loss = 1
		lossnew = 0
		while loss > delta and n_iter < maxiter:
			lossold = lossnew
			
			ypre = np.dot(newdata, theta)
			ypre2=ypre
			ypre2[ypre2<-708]=-708
			ypre2[ypre2>709]=709
			pi = np.exp(ypre2) / (1 + np.exp(ypre))

			w = pi * (1 - pi)
			wz = w * ypre + self.data[:, self.loci_size] - pi
			
			xtwx = np.zeros((theta_size, theta_size))
			for i in range(theta_size):
				for j in range(theta_size):
					red = np.sum(w * newdata[:, i] * newdata[:, j])
					xtwx[i, j] = red

			if (np.linalg.det(xtwx) != 0.0):
				xtwx = np.linalg.inv(xtwx)
			xwz = np.dot(newdata.T, wz)
			theta = np.dot(xtwx, xwz)
			
			lossnew = np.sum(np.abs(theta))
			loss = np.abs(lossnew - lossold)
			n_iter += 1
		
		pre = np.abs(1 - self.data[:, self.loci_size] - pi)
		aic = -2 * np.sum(np.log(np.maximum(pre, 0.0000001)))
			
		aic += 2 * theta_size		
		return aic
		
	def snp_factorial(self, n):
		z = 0.0
		if n < 0:
			print("Illegal n, n should be a non-negative number.")
			return 0
		if n == 0:
			return 0
		for i in range(1, int(n+1.0)):
			z += math.log(i)
		return z
