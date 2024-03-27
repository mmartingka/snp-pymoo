import numpy as np
from pymoo.core.problem import ElementwiseProblem

class SNPProblem(ElementwiseProblem):
	
	"""
	
	Attributes:
			sol_size: Number of SNP to be chosen
	"""
	
	def __init__(self, sol_size, loci_size, sample_size, data):
		super().__init__(n_var=1, n_obj=2, n_ieq_constr=0)
		self.sol_size = sol_size
		self.loci_size = loci_size
		self.sample_size = sample_size
		self.data = data

		
	def _evaluate(self, x, out, *args, **kwargs):
		score = self.bayesian_score(x)
		aic = self.logistic_score(x)
		out["F"] = np.array([+ score, + aic], dtype=float)

		
	def bayesian_score(self, x):
	
		"""
		"""
		
		comb = int(np.power(3.0, self.sol_size))
		observed_values = np.zeros((2, comb))
		col_sum_table = np.zeros(comb)
		score = 0
		
		for i in range(comb):
			observed_values[0][i] = 0
			observed_values[1][i] = 0
			col_sum_table[i] = 0
		
		for i in range(self.sample_size):
			index = 0
			cont = True
			for j in range(self.sol_size):
				if self.data[i][x[j]] == 3:
					cont = False
					break
				else:
					power_value = int(np.power(3.0, (self.sol_size - 1 - j)))
					index += self.data[i][x[j]] * power_value
			
			if cont:
				observed_values[self.data[i][self.loci_size]][index] += 1
				col_sum_table[index] += 1
		
		for i in range(comb):
			score += self.snp_factorial(observed_values[0][i]) 
			+ self.snp_factorial(observed_values[1][i]) 
			- self.snp_factorial(col_sum_table[i] + 1)

		score = np.abs(score)
		return score
		
	def logistic_score(self,x):
	
		"""
		"""
		delta = 0.001
		maxiter = 30
		aic = 0
		testsample = self.data.shape[0]
		theta_size = self.sol_size + 2
		newdata = np.zeros((testsample, theta_size))
		
		# newdata inicialization
		for i in range(testsample):
			newdata[i][0] = 1
			newdata[i][theta_size - 1] = 1
			for h in range(1, theta_size - 1):
				newdata[i][h] = self.data[i][x[h - 1]]
				newdata[i][theta_size - 1] *= newdata[i][h]

		theta = np.zeros(theta_size)
		n_iter = 0
		loss = 1
		lossnew = 0
		
		while loss > delta and n_iter < maxiter:
			lossold = lossnew
			ypre = np.dot(newdata, theta)
			
			pi = np.exp(ypre) / (1 + np.exp(ypre))
			pi[ypre < -708] = 0
			pi[ypre > 709] = 1
			
			w = pi * (1 - pi)
			wz = w * ypre + self.data[:, self.loci_size] - pi
			
			xtwx = np.zeros((theta_size, theta_size))
			for i in range(theta_size):
				for j in range(theta_size):
					red = np.sum(w * newdata[:, i] * newdata[:, j])
					xtwx[i, j] = red
				
			xtwx_inv = np.linalg.inv(xtwx)
			xwz = np.dot(newdata.T, wz)
			theta = np.dot(xtwx_inv, xwz)
			
			lossnew = np.sum(np.abs(theta))
			loss = np.abs(lossnew - lossold)
			n_iter += 1
		
		for s in range(testsample):
			pre = np.abs(1 - self.data[s, self.loci_size] - pi[s])
			aic += -2 * np.log(pre)
			
		aic += 2 * theta_size		
		return aic
		
	def snp_factorial(self, n):
		z = 0
		if n < 0:
			print("Illegal n, n should be a non-negative number.")
			return 0
		if n == 0:
			return 0
		for i in range(1, int(n) + 1):
			z += np.log(i)
		return z
