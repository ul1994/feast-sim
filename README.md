
## FEAST Code and Simulation

### Main Simulation

The main simulation details computing the R^2 score of mixing proportions.
This is done for sources at varying JSD values to understand the effect of
estimating the proportions when sources are similar/dissimilar.

The main simulation Rscript is in `sims/main.R`.
The most recent results were obtained using the following experiments below.

The folders `saved/jsd/10k/{jsd_value}` contains groups of sources which were
found to have the specified JSD values. These are input to the simulation script.

#### To run for JSD=0.08

```
Rscript sims/main.R 0.08 1000 saved/jsd/10k/0080/sources_jsd_0080_009649.rds
```

#### To run for JSD=0.125

```
Rscript sims/main.R 0.125 1000 saved/jsd/10k/0125/sources_jsd_0125_014478.rds
```

#### To run for JSD=0.5

```
Rscript sims/main.R 0.500 1000 saved/jsd/10k/0500/sources_jsd_0500_051307.rds
```

#### To run for JSD=0.90

```
Rscript sims/main.R 0.90 1000 saved/jsd/10k/0900/sources_jsd_0900_090180.rds
```

### Data Preparation


#### Sources with desired JSD

Simulations may be performed for sources selected to be similar to eachother in varying amounts. The following describes the procedure used to find such groups of sources and save them in  `saved/jsd/10k/{jsd_value}` for reproducable experiments.

The procedure is:

1. Precalculate JSD between pairs of sources in the data.

2. Sort sources by their calculated JSD in increasing order.

3. Sample pairs of sources close to where the desired JSD is in the ordered list. Keep the resulting groups who have approx. the desired JSD.

To run the procedure:

1. `scripts/find_jsd_pairs.R`
2. `scripts/group_by_jsd.R`