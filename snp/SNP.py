import pandas as pd
import numpy as np

class SNP:

	"""
	Basic structure to store the information read from the file studied.
	
	Attributes:
		file: Data file path 
		sample_size: Number of individuals in the sample (rows[1:])
		loci_size: Number of loci (SNP) of each individual (columns[:-1]) 
		class_column: Index column where class is store 
		data: Values retrieved from the file (vector/matrix)
		SNP_names: First row of the file (names for each SNP - vector)
		number_classes:  Number of individual in each class 
						(0-No fenotype , 1-Yes fenotype)
	"""
	
	def __init__(self, file):

		"""
		Init SNP class with data read from file studied 
		"""

		self.load_data(file)
		
	def load_data(self, file):
	
		"""
		Load SNP from parameter file. 
		"""
		
		df_SNP = pd.read_csv(file, sep=',')
		self.SNP_names = np.asarray(df_SNP.columns[:-1])
		self.loci_size = len(self.SNP_names)
		self.data = df_SNP[self.SNP_names].to_numpy()
		self.sample_size = df_SNP.shape[0]
		self.class_column = len(df_SNP.columns)
		self.number_classes = np.asarray(
								df_SNP.value_counts(
									df_SNP.columns[-1:][0]))
