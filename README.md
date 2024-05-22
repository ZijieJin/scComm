# scComm
scComm is a computational pipeline for inferring cell-cell communication at single-cell resolution. scComm works on R platform with any OS. If you have any questions related to scComm, please visit https://github.com/ZijieJin/scComm and post them on the Issues page or email me: jzj2035198@outlook.com

## Software Prerequisite

- R (tested on 4.2.1)
- R package: stringr, entropy, Seurat (tested on V4), pracma, RcppML, NMF, matrixStats, progress

## Recommend Configuration

- 32 GB memory or more 

- 4 CPU cores or more

## Data Requirement

- Single-cell RNA-seq gene expression data (10X or Smart-seq) with cell annotations

- (optional) custom Ligand-receptor database

- (optional) custom downstream regulation database

## Usage

`source(scComm.R)`

`res = scComm(data, anno)`

## Result Interpretation



## Commercial use

For non-academic use, please email Dr. Jin (jinzijie@bjmu.edu.cn) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
