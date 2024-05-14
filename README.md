[article]: https://www.nature.com/articles/s41598-022-11522-z

# Economic-Systemic-Risk
*Quantifying firm-level economic systemic risk from supply networks*

To know about the detail of the method used, please refer to the following scientific publication
[Quantifying firm-level economic systemic risk from nation-wide supply networks][article] by C. Diem et al. </br>

## Basic instructions
1) Install Julia by running `curl -fsSL https://install.julialang.org | sh` in the command line or go to [https://julialang.org/downloads/].

2) Run the script `./instantiate.sh` from in the command line to install the necessary libraries.

3) Prepare the input CSV file with the format found in the example data/test-list.csv

4) Execute `./run.sh inputfile outputfile` to estimate ESRI with a sequential code.

5) Execute `./run-parallel.sh inputfile outputfile` to estimate ESRI with a parallelized code.

<br/><b>Good luck.</b>
