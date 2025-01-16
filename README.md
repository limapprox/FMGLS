# Fully modified GLS estimation: Replication code and data
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Authors: Yicong Lin (yc.lin@vu.nl, Vrije Universiteit Amsterdam & Tinbergen Institute) and Hanno Reuvers (hannoreuvers@gmail.com, Erasmus Universiteit Rotterdam & Tinbergen Institute)

Discussion paper titled "Fully Modified Estimation in Cointegrating Polynomial Regressions: Extensions and Monte Carlo Comparison" is available at https://tinbergen.nl/discussion-paper/6215/22-093-iii-fully-modified-estimation-in-cointegrating-polynomial-regressions-extensions-and-monte-carlo-comparison.  

## Raw data
- primary balance-to-GDP ratios: "imf-dm-export-updated_pb.xls", available at https://www.imf.org/external/datamapper/pb@FPP/USA?year=2022;
- public debt-to-GDP ratios: "imf-dm-export-updated_debt.xls", available at https://www.imf.org/external/datamapper/d@FPP/USA?year=2022.

## Extracted data
- Data used for the empirical exercise: "pb_d_joint_table.csv."

## Replication code
- Run "ReplicationFRF_cubic.m" to replicate the empirical exercise of fiscal reaction functions.

## Main code
- "AndHAC.m": Andrew's HAC estimator of long-run covariance matrices;
- "fm_inference_cubic.n": Fully modified estimation and inference, including the GLS estimator proposed in this project, along with the system OLS (SOLS) and seemingly unrelated regression (SUR) methods proposed by Wagner et. al. (2020, https://www.sciencedirect.com/science/article/abs/pii/S0304407619301150);  
- "LRbiam.m": BIAM estimator of long-run covariance matrices proposed in this project;
- "multikpss_bonferroni.m": Multivariate KPSS test proposed in this project;
- 



For any questions or feedback, please feel free to contact the authors via email: 
yc.lin@vu.nl;
hannoreuvers@gmail.com.