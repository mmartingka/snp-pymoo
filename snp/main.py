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

from dask.distributed import Client

from pymoo.algorithms.moo.nsga3 import NSGA3
#from pymoo.algorithms.moo.age import AGEMOEA

from pymoo.util.ref_dirs import get_reference_directions
from pymoo.termination import get_termination
from pymoo.optimize import minimize
from pymoo.core.problem import DaskParallelization

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
	
	# Create sampling from readed data
	sampling = SNPSampling()
	
	# Create crossover operation given its probability
	crossover = SNPCrossover(prob = prob_cross)

	# Create mutation operation given its probability and range of mutation
	mutation = SNPMutation(prob_mutation = ((100+mlambda)/dim_epi), 
						range_mut = (snp.loci_size*mfactor/100))

	# Create duplication operation
	duplication = SNPDuplicateElimination()
	
	# Creates a client with a define number of nodes and threads per node 
	""" 
	When running multiple pymoo executions simultaneously, it is better to 
	assign a separate node to each execution to avoid shared memory usage.
	"""
	client = Client(n_workers=4, threads_per_worker=2)
	client.restart()
	print("DASK STARTED")
	
	# Initialize the thread pool and create the runner
	runner = DaskParallelization(client)

	# Define the problem by passing the starmap interface of the thread pool
	problem = SNPProblem(elementwise_runner=runner,
					  dim_epi = dim_epi, 
					  loci_size = snp.loci_size, 
					  sample_size = snp.sample_size, 
					  data = snp.data)

	# Create the reference directions to be used for the optimization
	ref_dirs = get_reference_directions("uniform", 2, n_points=12)

	# Create the repair operation 
	repair = SNPRepair()

	# Initialize algorithm with operators
	algorithm = NSGA3(ref_dirs=ref_dirs,
				   pop_size=pop_size, 
				   sampling=sampling, 
				   crossover=crossover,
				   mutation=mutation, 
				   repair=repair,
				   eliminate_duplicates=duplication)
	
	# Add termination limit
	termination = get_termination('n_gen', max_iter)

	# Execute pymoo 
	res = minimize(problem, algorithm, termination, seed=random.seed(142523623))

	formatted_list = [f"{x[0]},{x[1]}" for x in res.F]
	formatted_string = "\n".join(formatted_list)
	print(formatted_string)

	print("Elapsed time (seconds):", res.exec_time)
	client.close()
	print("DASK SHUTDOWN")


def parseArguments():

	""" Read arguments from command line """
	
	parser = argparse.ArgumentParser(description="Multi-objective SNP program")
	
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	# Population size argument
	required.add_argument("pop_size", 
						type=int, 
						help="Population size")
	
	# Maximum number of iterations argument
	required.add_argument("max_iter", 
						type=int, 
						help="Maximum number of iterations")
	
	# Lamda probability mutation argument
	required.add_argument("lambda_mutation", 
						type=int, 
						help="Lamda probability mutation [0-100]")
	
	# Factor range mutation argument
	required.add_argument("factor_mutation", 
						type=int, 
						help="Factor range mutation (pe. 10)")

	# File path argument
	required.add_argument("file",
						type=Path, 
						help="File path (00.1600.0.antesnp100.txt)")
						
	# Solution SNP size argument
	required.add_argument("dim_epi", 
						type=int, 
						help="Solution SNP size (2/5/8)")
						
	# Probability crossover argument
	optional.add_argument("prob_cross", 
						type=float, 
						help="Probability crossover (default 1.0)",
						default=1.0)
							
	return parser.parse_args()


if __name__ == "__main__":

	""" This is executed when run from the command line """
	# Initialize randomness
	random.seed(1425236234)
	np.random.seed(1425236234)

	args = parseArguments()
	main(args.pop_size, args.max_iter, args.lambda_mutation, 
	  args.factor_mutation, args.file, args.dim_epi, args.prob_cross)
