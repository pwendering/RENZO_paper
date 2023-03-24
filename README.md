# **R**elative **ENZ**yme in **O**ptimization of metabolic models (RENZO)

## Publication

## Requirements
- Matlab (tested with 2020a and b)
- [COBRA Toolbox 3.0]() (Heirendt et al. (2019), doi: [10.1038/s41596-018-0098-2](https://doi.org/10.1038/s41596-018-0098-2)
- [GECKO Toolbox 2.0](https://github.com/SysBioChalmers/GECKO) (Sanchéz et al. (2017), doi: [10.15252/msb.20167411](https://doi.org/10.15252/msb.20167411),
  Domenzain et al. (2021), doi: [10.1101/2021.03.05.433259](https://doi.org/10.1101/2021.03.05.433259))
- [Gurobi solver](https://www.gurobi.com/) (tested with version 9.1.1)

## Steps to reproduce results
- clone repository, change into *Code* subdirectory
- add gurobi solver to Matlab path
- _FVA\_CDC48_, _FVA\_Clb2_: flux variability analysis in dark and blue light conditions
- _enrichment\_analysis_: find enriched pathways in reactions with non-overlapping flux ranges

## RENZO documentation
Run RENZO: 
```
results = RENZO(model_mt, model_wt, msrd_protein_ids, Ptot_mt, Ptot_wt,...
    protein_abundance_ratio, biomass_ratio, GAM, NCPU, bio_name, f, ratio_tol)
```


| Input parameter | Class | Description |
|-----------|-------|-------------|
| model\_mt | struct | mutant-specific ecModel |
| model\_wt | struct | wild-type-specifc ecModel |
| msrd\_protein\_ids | cellstr | Swissprot IDs of measured proteins |
| Ptot\_mt | double | total protein content (mutant)
| Ptot\_wt | double | total protein content (wild type)
| protein\_abundance\_ratio | double | ratios of protein abundances (mutant/wild type) (number of rows must correspond to the number of Swissprot IDs) |
| biomass\_ratio | double | ratio between mutant and wild type biomass |  
| GAM | struct | fields: mut, wt, contain growth-associated maintenance for mutant and wild type |
| NCPU | double | number of workers for parallel pool (optional, default: 1) | 
| bio\_name | char | name of biomass reaction in model (optional, default: biomass) |
| f | double | fraction of protein mass that is accounted for in the ecModel (optional, default: 0.5) |
| ratio_tol | double | tolerance for the ratio between predicted and measured protein ratios (optional, default: 1e-5) |


*results* is a struct with fields

| Field | Content |
| ----- | ------- |
| FVA\_min | reaction flux minima |
| FVA\_max | reaction flux maxima |
| infeasible\_min | number of infeasible minimization programs |
| infeasible\_max | number of infeasible maximization programs |
| distinct\_rxns | reactions with non-overlapping flux ranges |
| genes\_in\_distrinct\_rxns | genes associated with distinct\_rxns |
