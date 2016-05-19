# SCOBI
Salmonid COmposition Bootstrap Intervals

`SCOBI` is an R package for performing compositional analyses of adults and smolts at Lower Granite Dam. Adults are analyzed using 
the function `SCOBI()`. Juveniles are analyzed using the function `SCRAPI()`. It also contains the useful functions `lgr2SCOBI()` 
and `lgr2SCRAPI()` for formatting raw data from the Lower Granite Dam trapping database (LGTrappingDB) for input into the `SCOBI()`
and `SCRAPI()` functions.

## Getting Started

To install `SCOBI` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:
```
install.packages("devtools")
library(devtools)
```
Once `devtools` is successfully installed, use the following to install SCOBI:
```
devtools::install_github("mackerman44/SCOBI")
```
If you are interested in making contributions, consider getting a GitHub account, fork this repository, modify, and send a pull request.

For further information, see:

Steinhorst, K. T. Copeland, M. W. Ackerman, W. C. Schrader, and E. C. Anderson. (*In review*) Estimates and Confidence Intervals
for Run Composition of Returning Salmonids. Fishery Bulletin.

Enjoy!
