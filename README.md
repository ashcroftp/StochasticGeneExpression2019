# Stochastic gene expression influences the selection of antibiotic resistance mutations

Lei Sun, Peter Ashcroft, Martin Ackermann, and Sebastian Bonhoeffer

## Directory structure

### aa-ms/
All manuscript related files, including the actual manuscript, figures, and supplementary material.


### bb-code/
All files related to simulations, analysis, and plotting.
This contains directories:
* `Benchmarks/`: Scripts to test the performance of the simulation functions;
* `ControlNoise/`: Simulation configs, data, and script for modifying the expression noise;
* `EffluxBias/`: Simulation configs, data, and script for vary the bias of efflux pump distribution following cell division;
* `GrowthTest/`: Simulation configs, data, and script for comparing MIC measurements between growth and lineage survival data;
* `MIC-fraction/`: Functions used to calculate the MIC-fraction of bound targets for each drug type;
* `MoleculeDistributions/`: Functions to compute the distributions of within-cell molecules across multiple cells;
* `MoleculeTrajectories/`: Functions to calculate and plot the trajectories of within-cell molecules in the presence and absence of drugs;
* `Params/`: Files which hold parameters taken from the literature;
* `PopulationGrowth/`: Functions for calculating the populations growth rates;
* `PopulationSurvival/`: Functions to simulate populations of cells, from which we can extract survival probabilities and extinction times in response to pulsed and continuous treatments;
* `R/`: Directory for all commonly used functions, including simulation, ODE, analysis, and plotting;
* `Simulate/`: Files which are used to standardise the execution and output of simulations;
* `Systematic/`: Simulation configs, data, and script for the parameter grid search experiment;

## Running simulations
Simulations are executed using the function `Simulate/simulate-script.R`, which is ran from the command line.
It requires a configuration file to be passed as an argument, and this file must contain specific information.
To run a simulation as described in the configuration file `MoleculeDistributions/config-WT.R`, we run
```bash
R --vanilla --slave < Simulate/simulate-script.R --args MoleculeDistributions/config-WT.R
```
Note that all file paths are relative to the parent directory `bb-code/`.
Alternatively, we can use the shortcut shell script `simulate.sh` to run the simulations, i.e.
```bash
./simulate.sh MoleculeDistributions/config-WT.R
```

Configuration files must contain the following information:
* `drug`: one of "ciprofloxacin" or "rifampicin";
* `cell.type`: one of "KO", "WT", "REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT";
* `output`: what type of output do we want from the simulations? One of "survival", "growth", "molecules";
* `total.sims`: the number of simulations to perform per parameter set
* `sims.per.node`: how many simulations are sent to a given node? If zero, then simulations are run locally, otherwise this runs on EULER (HPC cluster at ETH Zurich);
* `queue`: which HPC queue are simulations sent to? One of "04:00" or "24:00" (not important for non-HPC sims);
* `var.params`: A data frame of parameters to scan over. These names must match those as used in the model, and can be appended with `.MULT` if we wish to simply multiply a parameter value by a factor;
* `file.id`: A unique label which identifies the dataset. This is appended with the parameter and node indices;
* `data.directory`: The location where the data will be stored (relative to the parent directory `bb-code`).

Once the data is generated, it can be combined into a single file (with full details also archived) by running
```bash
./compress.sh MoleculeDistributions/config-WT.R
```

Sometimes submitted jobs may fail.
To identify if any files are missing, we can run the `Simulate/missing-files-runs.R` script, with the configuration file passed as an argument.
I.e.,
```bash
R --vanilla --slave < Simulate/missing-files.R --args MoleculeDistributions/config-WT.R
```
To re-run the missing simulations, add `run` as an extra argument to the above command:
```bash
R --vanilla --slave < Simulate/missing-files.R --args MoleculeDistributions/config-WT.R run
```
Again, a shell-script can be used to shorten these commands:
```bash
./missing-files.sh MoleculeDistributions/config-WT.R [run]
```

### Example: Computing a lineage survival probability curve
Within the `PopulationSurvival/` directory, we can sweep over the mic-fraction parameter range to find out which one gives us a net survival probability of 90%.
The configuration file for this example could look like this:
```r
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "survival"
total.sims <- 1000
sims.per.node <- 200
queue <- "04:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    mutant.effect = c(1),
    drug.outsideConcentration.MULT = c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5))
)
file.id <- paste(drug, cell.type, sep = "_")
data.directory <- "PopulationSurvival/zz-data/"
```

To simulate this, we run:
```bash
./simulate.sh PopulationSurvival/config-ciprofloxacin_WT.R
```
Once generated, we compress the data:
```bash
./compress.sh PopulationSurvival/config-ciprofloxacin_WT.R
```
The combined files (`PopulationSurvival/zz-data/ciprofloxacin_WT_output.dat`) can then be read and plotted using `PopulationSurvival/plot-survival.R`.
