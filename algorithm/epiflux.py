import argparse
import os
#import sys
#sys.path.append("../metabolic-models/")
from codEpiflux import *

# Description:
# ------------
parser = argparse.ArgumentParser(description='Welcome to Epiflux.')

# Arguments:
# ----------
parser.add_argument("-i","--input_file", help="Name of the csv file with required information.", type=str)
parser.add_argument("-o","--output_directory", help="Name of the directory where output files will be stored.", type=str)

args = parser.parse_args()

print("Input file: ", args.input_file)

# Create an output directory:
# ---------------------------
# output directory:
resultsDir = args.output_directory
if os.path.exists(resultsDir) == False:
	path = os.path.join(resultsDir)
	os.mkdir(path)

# Run Pheflux:
# ------------
inputFile = args.input_file
getFluxes(inputFile,resultsDir) #(-i, -o)


print("Output directory:", resultsDir)
