# Stochastic gene expression influences the selection of antibiotic resistance mutations

Lei Sun, Peter Ashcroft, Martin Ackermann, and Sebastian Bonhoeffer

## Directory structure

### aa-ms/
All manuscript related files, including the actual manuscript, figures, and supplementary material.


### bb-code/
All files related to simulations, analysis, and plotting.
This contains directories:
* `Params/`: Files which hold parameters taken from the literature;
* `R/`: Directory for all commonly used functions, including simulation, ODE, analysis, and plotting;
* `Simulate/`: Files which are used to standardise the execution and output of simulations;
* `MIC-fraction/`: Functions used to calculate the MIC-fraction of bound targets for each drug type;
* `MoleculeTrajectories/`: Functions to calculate and plot the trajectories of within-cell molecules in the presence and absence of drugs;
* `MoleculeDistributions/`: Functions to compute the distributions of within-cell molecules across multiple cells;
* `PopulationGrowth/`: Functions for calculating the populations growth rates;
* `PopulationSurvival/`: Functions to simulate populations of cells, from which we can extract survival probabilities and extinction times in response to pulsed and continuous treatments;
* `Doc/`: Generated figures and the manuscript in R-markdown format;
* `Benchmarks/`: Scripts to test the performance of the simulation functions;
* `zz-Other/`: Other scripts that don't really have a home.

## Running simulations
Simulations are executed using the function `Simulate/simulate.R`, which is ran from the command line.
It requires a configuration file to be passed as an argument, and this file must contain specific information.
To run a simulation as described in the configuration file `Simulate/config-template.R`, we run
```bash
R --vanilla --slave < Simulate/simulate.R --args Simulate/config-template.R
```
Note that all file paths are relative to the parent directory `bb-code/`.
Alternatively, we can use the shell script `simulate.sh` to run the simulations, i.e.
```bash
./simulate.sh Simulate/config-template.R
```

Configuration files must contain the following information:
* `drug`: one of "ampicillin", "cefotaxime", "ciprofloxacin", "rifampicin", "tetracycline";
* `cell.type`: one of "KO", "WT", "REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT";
* `output`: what type of output do we want from the simulations? One of "survival", "growth", "molecules";
* `total.sims`: the number of simulations to perform per parameter set
* `sims.per.node`: how many simulations are sent to a given node? If zero simulations are run locally, otherwise this runs on EULER (cluster);
* `var.params`: A data frame of parameters to scan over. These names must match those as used in the model, and can be appended with `.MULT` if we wish to simply multiply a parameter value by a factor;
* `file.id`: A unique label which identifies the dataset. This is appended with the parameter and node indices;
* `data.directory`: The location where the data will be stored (relative to the parent directory `bb-code`).

They can be easily created by copying and modifying `Simulate/config-template.R`.

Once the data is generated, it can be combined into a single file (with full details also archived) by running
```bash
./simulate.sh Simulate/config-template.R compress
```
The files in `Simulate` should not be modified.

Sometimes submitted jobs may fail.
To identify if any files are missing, we can run the `Simulate/missing-files.R` script, with the configuration file passed as an argument.
I.e.,
```bash
R --vanilla --slave < Simulate/missing-files.R --args Simulate/config-template.R
```
To re-run the missing simulations, add `run` as an extra argument to the above command:
```bash
R --vanilla --slave < Simulate/missing-files.R --args Simulate/config-template.R run
```
Again, a shell-script can be used to shorten these commands:
```bash
./missing-files.sh Simulate/config-template.R [run]
```

### Example: Computing the MIC fraction
Within the `MIC-fraction/` directory, we can sweep over the mic-fraction parameter range to find out which one gives us a net survival probability of 50%.
The configuration file for this example looks like this:
```r
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "survival"
total.sims <- 10000
sims.per.node <- 100
var.params <- expand.grid(
    drug.outsideConcentration.MULT = 1,
    MIC.fraction = seq(0.06, 0.13, by = 0.005)
)
file.id <- paste(drug, cell.type, sep = "_")
data.directory <- "MIC-fraction/mic-data/"
```

To simulate this, we run one of the following commands:
```bash
R --vanilla --slave < Simulate/simulate.R --args MIC-fraction/config-ciprofloxacin.R
```
```bash
./simulate.sh MIC-fraction/config-ciprofloxacin.R
```
Once generated, we compress the data:
```bash
./simulate.sh MIC-fraction/config-ciprofloxacin.R compress
```
The combined files (`MIC-fraction/mic-data/ciprofloxacin_WT_combined.dat`) can then be read and plotted using `MIC-fraction/plot-micFraction.R`.
