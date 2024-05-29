# algorithm
The use of `epiflux.py` is exemplified in the following files:
- tutorial.ipynb: Tutorial of `epiflux.py` using Jupyter notebook.
- example.sh: Example run of Epiflux.py using Linux bash script. It calls the `codEpiflux.py` file.

Both examples use the information stored in the `data` folder and save the outputs in the `results` folder. The inputs and outputs are described below.

## Inputs
Four inputs are needed:
- ```inputFileName```       : Name of the csv file with required information.
- ```resultsDir```      : Name of the directory where output files will be stored.

### Input file
The input file is a comma-separated file containing: the name of the organism, the condition, medium, gene expression file, and genome-scale metabolic network. 
More than one row can be defined in the same file. Pheflux will be run on the data provided for each line. For example, the input file ```examples/data/InputData_yeast_example.csv``` contains the following:
| **Organism** | **Condition**        | **GeneExpFile**                                          | **Medium**                                      | **Network**                          | **KnownFluxes**                      |                        
|--------------|----------------------|----------------------------------------------------------|-------------------------------------------------|--------------------------------------|--------------------------------------|
| Scerevisiae  | T1                   | ../examples/data/transcriptomes/T1.csv                               | ../examples/data/mediums/Scerevisiae_Mediumfile_Kuang2.csv  | ../examples/data/gems/iMM904_ac_me_glycogen.xml  | ../examples/data/knownFluxes/knownFluxesT1.csv  |

### Results directory
A string describing the folder where the output files will be stored. 

## Outputs
Two output files are generated:
- ```Organism_Condition_STATUS.fluxes.csv```: Comma-separated file containing the fluxome estimation. 
- ```Organism_Condition_STATUS.Log.csv```: Comma-separated file containing various statistics.
For both files, the words ```Organism```, and ```Condition``` are extracted from the ```examples/data/InputData_yeast_example.csv``` file, and the ```STATUS``` word indicates ```IPOPT``` finalization condition.
