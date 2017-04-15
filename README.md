# Sharpe Ratio Portfolio Optimization via Quadratic Programming

This repository contains a set of scripts that perform (constrained) Sharpe Ratio portfolio optimization by casting the original quasi-convex Sharpe ratio maximization problem as a convex program (i.e. a quadratic program).

### Requirements

In order to use the sharpe ration maximization scripts in this repository:
- You must be using Mac OSX or Linux
- You must have installed [`R` programming language](https://www.r-project.org/), version 3 or later
- You must have installed the [quadprog](https://cran.r-project.org/web/packages/quadprog/quadprog.pdf) library
- You must navigate to this repository after cloning it 
```
git clone https://…/sharpemax.git
cd ./sharpemax
```
before executing any of the scripts therein

### Unconstrained

Format your returns data as a `csv` that has the asset names as the column headers and the returns down each column, however do not include a date column and make sure instead of percent `(12, -3, …)` you use raw numbers `(.12, -.03, …)`. To obtain the optimal portfolio weights, simply run the following in command line

```
Rscript --vanilla ./maxSharpe.R /path/to/your/data.csv
```

and it will output `maxSharpe.csv` containing portfolio weights to the root of this repository.

### Robust Sharpe Ratio

The sharpe ratio is the ratio between the mean and variance. However, the simple arithmetic estimators for both statistics are not impervious to outliers or extreme events (which are common in financial time series) and, in the case of covariance, underestimates non-linear monotonic relationships. The provided script `maxSharpeRobust.R` functions exactly like `maxSharpe.R`,

```
Rscript --vanilla ./maxSharpeRobust.R /path/to/your/data.csv
```

however, internally it uses the truncated mean in place of the arithmetic mean and [Spearman’s rank covariance](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) in place of the arithmetic covariance. This ultimately results in a Sharpe ratio that is less statistically efficient, but significantly more robust.

### Imposing Constraints

The script `maxSharpeConstrained.R` can process asset constraints, sector constraints, and class constraints. You can call it as follows

```
Rscript --vanilla maxSharpeConstrained.R ./returns.csv ./assetConstraints.csv ./sectorConstraints.csv ./classConstraints.csv
```

where the `csv`’s follow the same format as the files contained in the included `testData` directory. The script will write a new `csv` with weights and another with performance metrics, in the same as `maxSharpe.R`. Leave constraints blank and they will be ignored. The fewer the constraints, the faster the program will typically run.

### Comparing Several Sets of Constraints

The script `maxSharpeConstrainedLoop.R` runs `maxSharpeConstrained.R` many times on a set of strategies you define. Specifically, you must fill a directory the same way as the folder `loopTestData`; containing one `returns.csv` but many different constraint files according to your strategy. They must also abide by the same naming convention as the files in the `loopTestData` folder (i.e. `assetConstraints_1.csv`, `sectorConstraints_2.csv`, and so forth for as many strategies as you'll be considering). Then, you call the script from command line as follows

```
Rscript --vanilla maxSharpeConstrainedLoop.R …/directory/with/constraint/files/ totalNumberOfStrategies
```

Where `totalNumberOfStrategies` is a positive integer specifying how many sets of constraints you have in your constraints directory. This will then run the algorithm through your strategies and output a file `sharpeLoops.csv` with strategies listed along the columns and metrics listed along the rows.

### Sortino Ratio Maximization

Additionally included is a script `maxSortinoConstrainedLoop.R` which attempts to maximize the [Sortino](http://www.investopedia.com/terms/s/sortinoratio.asp) ratio instead of the Sharpe ratio, using [Estrada’s method](http://webprofesores.iese.edu/jestrada/PDF/Research/Refereed/MSO.pdf) for estimating the semi-covariance matrix. If the semi-covariance matrix is indefinite, it is projected onto the semidefinite cone so that the underlying optimization problem remains convex.

The `maxSortinoConstrainedLoop.R` script functions just like `maxSharpeConstrainedLoop.R`

```
Rscript --vanilla maxSortinoConstrainedLoop.R …/directory/with/constraint/files/ totalNumberOfStrategies
```

### Note on Formatting Data

In all these scripts the algorithm takes returns and all constraints as raw numbers (so that if you multiplied them by 100, you get percent) and also outputs the weights as raw numbers. This way the allocation size doesn't affect how the constraints are interpreted: A minimum of .3 or 30% on asset 1 is the same irrespective of whether you're investing $10k or $1M. You can then multiply the weights by the investment amount and get the per-asset dollar allocations back.

Also, the assets/securities must have text names (column headers) or else R has some issues.
