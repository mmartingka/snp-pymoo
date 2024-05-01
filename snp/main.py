#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "María Martín Rodríguez"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
from pathlib import Path
import random
import numpy as np

from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.util.ref_dirs import get_reference_directions

from pymoo.optimize import minimize

from SNP import SNP
from snppymoo.sampling import SNPSampling
from snppymoo.crossover import SNPCrossover
from snppymoo.mutation import SNPMutation
from snppymoo.repair import SNPRepair
from snppymoo.duplicateElimination import SNPDuplicateElimination
from snppymoo.problem import SNPProblem


def main(pop_size, max_iter, mlambda, mfactor, filePath, dim_epi, prob_cross):

	""" Main entry point of the app """

	# Read data from file and store in SNP object
	snp = SNP(file=filePath)
	
	sampling = SNPSampling()
	
	crossover = SNPCrossover(prob = prob_cross)

	mutation = SNPMutation(
		prob_mutation = ((100+mlambda)/dim_epi), 
		range_mut = (snp.loci_size*mfactor/100)
		)

	duplication = SNPDuplicateElimination()

	problem = SNPProblem(
		dim_epi = dim_epi, 
		loci_size = snp.loci_size, 
		sample_size = snp.sample_size, 
		data = snp.data)

	# create the reference directions to be used for the optimization
	ref_dirs = get_reference_directions("uniform", 2, n_points=12)

	repair = SNPRepair()

	algorithm = NSGA2(
		#ref_dirs=ref_dirs,
		pop_size=pop_size, 
		sampling=sampling, 
		crossover=crossover,
		mutation=mutation, 
		repair=repair,
		eliminate_duplicates=duplication)

	res = minimize(problem, algorithm, ('n_gen', max_iter), 
				seed=1, verbose=True)
	print(res.X)
	print(res.F)



def parseArguments():

	""" Read arguments from command line """
	
	parser = argparse.ArgumentParser(
		description="Multi-objective SNP program")
	
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	# Population size argument
	required.add_argument("pop_size", 
						type=int, 
						help="Population size argument")
	
	# Maximum number of iterations argument
	required.add_argument("max_iter", 
						type=int, 
						help="Maximum number of iterations argument")
	
	# Lamda probability mutation argument
	required.add_argument("lambda_mutation", 
						type=int, 
						help="Lamda probability mutation argument [0-100]")
	
	# Factor range mutation argument
	required.add_argument("factor_mutation", 
						type=int, 
						help="Factor range mutation argument (pe. 10)")

	# File path argument
	required.add_argument("file",
						type=Path, 
						help="File path argument (00.1600.0.antesnp100.txt)")
						
	# Solution SNP size argument
	required.add_argument("dim_epi", 
						type=int, 
						help="Solution SNP size argument (2/5/10)")
						
	# Probability crossover argument
	optional.add_argument("prob_cross", 
						type=float, 
						help="Probability crossover argument (default 1.0)",
						default=1.0)
							
	
	return parser.parse_args()


if __name__ == "__main__":

	""" This is executed when run from the command line """
	random.seed(1425236234)
	np.random.seed(1425236234)

	args = parseArguments()
	print(args)
	main(args.pop_size, args.max_iter, args.lambda_mutation, 
	  args.factor_mutation, args.file, args.dim_epi, args.prob_cross)
