#!/usr/bin/env python
# coding: utf-8

# # Tutorial

# ## 1. Loading Epiflux
# 
# The use Epiflux, go to the folder where the ```Epiflux.py``` is located and load it into memory as follows:

# In[3]:


from codEpiflux import *


# ## 2. Preparing the input files
# 
# To run Pheflux, four inputs are needed:
# 
# - ```inputFileName```       : Name of the csv file with required information.
# - ```resultsDir```      : Name of the directory where output files will be stored.

# ### 2. 1 Input file
# 
# The input file is a comma-separated file containing: the name of the organism, the condition, medium, gene expression file, and genome-scale metabolic network. 
# 
# More than one row can be defined in the same file. Pheflux will be run on the data provided for each line. For example, the input file ```examples/data/InputData_yeast_example.csv``` contains the following:
# 
# | **Organism** | **Condition**          | **GeneExpFile**                        | **Medium**                                       | **Network**                         | **KnownFluxes**                      |
# |--------------|------------------------|----------------------------------------|--------------------------------------------------|-------------------------------------|--------------------------------------|
# | Scerevisiae  | T1                     | ../examples/data/transcriptomes/T1.csv             | ../examples/data/mediums/Scerevisiae_Mediumfile_Kuang2.csv   | ../examples/data/gems/iMM904_ac_me_glycogen.xml | ../examples/data/knownFluxes/knownFluxes_T1.csv  |

# In[4]:


inputFile       = "./data/InputData_yeast_example.csv"  # Name of the csv file with required information


# ### 2.2 Results directory
# 
# A string describing the folder where the output files will be stored. **Important: In the "example.sh" script the output folder must be created before predictions are made.**

# In[5]:


resultsDir      = "./results" # Name of the directory where output files will be stored


# ## 3. Running Epiflux
# 
# With all input files defined, you are now ready to make fluxome estimations. For this, you need to used the ```getFluxes```:

# getFluxes(inputFile,resultsDir) 

# ## 4. What is the output
# 
# Two output files are generated:
# 
# - ```Organism_Condition_STATUS.fluxes.csv```: Comma-separated file containing the fluxome estimation. 
# - ```Organism_Condition_STATUS.Log.csv```: Comma-separated file containing various statistics.
# 
# For both files, the words ```Organism```, and ```Condition``` are extracted from the ```examples/data/InputData_example.csv``` file, and the ```STATUS``` word indicates ```IPOPT``` finalization condition.
