#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "María Martín Rodríguez"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
from pathlib import Path

from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize

from snp.SNP import SNP
from snppymoo.sampling import SNPSampling
from snppymoo.crossover import SNPCrossover
from snppymoo.mutation import SNPMutation
from snppymoo.duplicateElimination import SNPDuplicateElimination
from snppymoo.problem import SNPProblem


#"./00.1600.0.antesnp100.txt"
def main(filePath, sol_size, lambda_mutation, factor_mutation, prob):

	""" Main entry point of the app """

	# Read data from file and store in SNP object
	snp = SNP(file=filePath)
	print(snp.sampleSize)
	
	sampling = SNPSampling()
	
	crossover = SNPCrossover(prob = prob)
	
	mutation = SNPMutation(
		prob_mutation = (sol_size*(100+lambda_mutation)), 
		range_mut = (snp.loci_size*factor_mutation/100))
	
	duplication = SNPDuplicateElimination()
	
	problem = SNPProblem(
		sol_size = sol_size, 
		loci_size = snp.loci_size, 
		sample_size = snp.sample_size, 
		data = snp.data)
	
	algorithm = NSGA3(
		pop_size=11, 
		sampling=sampling, 
		crossover=crossover,
		mutation=mutation, 
		eliminate_duplicates=duplication)
	
	res = minimize(problem, algorithm, ('n_gen', 100), seed=1, verbose=False)


def parseArguments():

	""" Read arguments from command line """
	
	parser = argparse.ArgumentParser(description=
										"Multi-objective SNP program")
	
	# Required file path positional argument
	parser.add_argument("file", 
						metavar="F", 
						type=Path, 
						help="Required file path positional argument")
						
	# Required solution SNP size argument
	parser.add_argument("sol_size", 
						metavar="S", 
						type=int, 
						help="Required solution SNP size argument")
						
	# Required lamda probability mutation argument
	parser.add_argument("lambda_mutation", 
						metavar="P", 
						type=float, 
						help="Required lamda probability mutation argument")
	
	# Required factor range mutation argument
	parser.add_argument("factor_mutation", 
						metavar="R", 
						type=float, 
						help="Required factor range mutation argument")
	
	# Non required probability crossover argument
	parser.add_argument("prob", 
						metavar="C", 
						type=float, 
						required=False,
						help="Non required probability crossover argument")
							
	
	args = parser.parse_args()
	return args.file, args.sol_size, args.lambda_mutation, args.factor_mutation, args.prob

if __name__ == "__main__":

	""" This is executed when run from the command line """

	main(parseArguments())
