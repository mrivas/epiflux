# Examples
The use of `epiflux.py` is exemplified in the following files:
- tutorial.ipynb: Tutorial of `epiflux.py` using Jupyter notebook.
- example.sh: Example run of Epiflux.py using Linux bash script. It calls the `epiflux-terminal.py` file.

Both examples use the information stored in the `data` folder and save the outputs in the `results` folder. The inputs and outputs are described below.

## Inputs
Four inputs are needed:
- ```inputFile```       : Name of the csv file with required information.
- ```resultsDir```      : Name of the directory where output files will be stored.
- ```verbosity```       : Verbose mode.
- ```prefix_log_file``` : Prefix log output file.
### Input file
The input file is a comma-separated file containing: the name of the organism, the condition, medium, gene expression file, and genome-scale metabolic network. 
More than one row can be defined in the same file. Pheflux will be run on the data provided for each line. For example, the input file ```examples/data/InputData_yeast_example.csv``` contains the following:
| **Organism** | **Condition**        | **GeneExpFile**                                          | **Medium**                                      | **Network**                          | **KnownFluxes**                      |                        
|--------------|----------------------|----------------------------------------------------------|-------------------------------------------------|--------------------------------------|--------------------------------------|
| Scerevisiae  | T1                   | data/transcriptomes/T1.csv                               | data/mediums/Scerevisiae_Mediumfile_Kuang2.csv  | data/gems/iMM904_ac_me_glycogen.xml  | data/knownFluxes/knownFluxesT1.csv  |

### Results directory
A string describing the folder where the output files will be stored. 
### Verbosity
Indicates whether a verbose output (```verbosity = True```) will be printed on terminal or not (```verbosity = False```).
### Prefix of Log file
Pheflux also produces a log file with summary statistics of all computations. The user-provided prefix is used to name this file. For example, if ```prefix = example``` an ```example_record_XXXX.log.csv``` file will be created, where ```XXXX``` is a random four-character tag.
## Outputs
Two output files are generated:
- ```Organism_Condition_STATUS.fluxes.csv```: Comma-separated file containing the fluxome estimation. 
- ```Organism_Condition_STATUS.Log.csv```: Comma-separated file containing various statistics.
For both files, the words ```Organism```, and ```Condition``` are extracted from the ```examples/data/InputData_yeast_example.csv``` file, and the ```STATUS``` word indicates ```IPOPT``` finalization condition.
