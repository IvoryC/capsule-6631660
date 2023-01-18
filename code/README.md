# Benchmarking and analysis for SCRuB

SCRuB is a decontamination method designed to help researchers address the common issue of contamination in microbiome studies. In this capsule, we provide full code and data to replicate the analyses from the SCRuB manuscript. 

## Code Walkthrough
The entire analysis from this capsule can be executed by running the `code/run` script, which calls the necessary functions from within the `code` directory in the proper order. We structure the rest of the code directory analysis as follows:
| Folder | Description |
|--|--|
| `SCRuB` | The Rscripts that make up SCRuB. Files from This directory are used in all analyses. As we continue to improve the `SCRuB` package, we recommend that readers who want to use SCRuB access it from the most up-to date package version, available [here](https://github.com/korem-lab/SCRuB). |
| `Simulation_Functions` |  The Rscripts that outline the core structure of our simulations. See the SCRuB manuscript Methods section, and supplementary figure S1 for further details. |
| `Fig1` | Contains an Rscript that runs simulations and compares both SCRuB and alternative decontamination methods. The total runtime has been shortened in this capsule by limiting the total number of simulation trials; to recreate the original analysis, set `n_experiments <- 10` in the `finalized_simulations` function call. Outputs from the simulation are stored in the `results/Fig_1` directory. |
| `Fig2` | Contains two scripts, each decontaminating one plate from the Minich et al. experiment (10.1128/mSystems.00186-19). Outputs form both files are stored in the `results/Fig_2/` directory. |
| `Plasma` | Contains an Rscript that deontaminates the plasma dataset from Poore et al. (10.1038/s41586-020-2095-1), stores sample shannon diversity csvs in `results/Plasma/`, runs the prediction pipeline form Poore et al. on each dataset, and stores the outputs in `results/Plasma/`. |
 Fig3 | Runs the decontamination and analyses of all figures corresponding to the experiment designed by the authors of SCRuB. Contains two files, one R script to run the decontaminations, and one python script to generate the plots. |
| `Tumor` | Contains: 1) an Rscript that runs decontamination pipelines on the tumor samples from Nejman et al. (10.1126/science.aay9189), storing resulting samples and diversity metrics in `results/Tumor`; 2) the `run_cancer_preds.py` file, which reads the decontaminated cancer data, and runs the melanome immune checkpoint inhibitor prediction analysis described in our Methods. To shorten the runtime, we provide a dictionnary containing the tuned hyperparameters for each dataset (calculated via our within-center k-fold tuning approach), which our function uses by default. To re-run the full tuning process, which will increase the runtime by several hours, set `rerun_tuning=True` in the `write_melanoma_rocs_nki_test_set` function call. In either case, running this file stores csvs of ROC curves into `results/Tumor`.|
| Plot_generation | Contains a single python script, `__init__.py`, which reads in all the outputs written into `results` from the previously mentioned folders, and creates the plots presented in our main figures. All plots are stored in the `results/Plots/` directory. |

## Data
All datasets analyzed in this study are publicly available. The college dormitory dataset used in Fig 1 is available from the European Nucleotide Archive (ENA), accession ERP115809, and Qiita, study ID 12470. Negative controls used in Fig 1 are available from Qiita, study ID 12019. The well-to-well leakage dataset is available from ENA, accession ERP115213. The plasma cfDNA data is available from ENA, accessions ERP119598, ERP119596, and ERP119597; and Qiita, study IDs 12667, 12691, and 12692. The tumor microbiome dataset is available from SRA, accession PRJNA624822. The processed data was obtained from Table S2 in the origina Nejman et al. publication.

All data used in our analysis is stored in the `data` directory. 

## Environment
The environment used in this analysis involes R 3.6.3 and python 3.6. We provide a docker file which can be used to set up the environment for a successful capsule run. 

## Support
If you have any questions regarding the analyis outlined in this capsule, or run into any issues when using SCRuB, please reach out via our [issues page](https://github.com/Shenhav-and-Korem-labs/SCRuB/issues).

