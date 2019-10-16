

### Main Simulation

The main simulation details computing the R^2 score of mixing proportions.
This is done for sources at varying JSD values to understand the effect of
estimating the proportions when sources are similar/dissimilar.

The main simulation Rscript is in `sims/main.R`.
The most recent results were obtained using the following experiments below.

The folders `saved/jsd/10k/{jsd_value}` contains groups of sources which were
found to have the specified JSD values. These are input to the simulation script.

#### To run for JSD=0.08

```bash
Rscript sims/main.R 0.08 500 saved/jsd/10k/0080/sources_jsd_0080_009649.rds
```

#### To run for JSD=0.125

```bash
Rscript sims/main.R 0.125 500 saved/jsd/10k/0500/sources_jsd_0500_045897.rds
```

#### To run for JSD=0.5

```bash
Rscript sims/main.R 0.500 500 saved/jsd/10k/0500/sources_jsd_0500_045897.rds
```

#### To run for JSD=0.90

```bash
Rscript sims/main.R 0.90 500 saved/jsd/10k/0900/sources_jsd_0900_090180.rds
```
