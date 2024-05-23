# scComm
scComm is a computational pipeline for inferring cell-cell communication at single-cell resolution. scComm works on R platform with any OS. If you have any questions related to scComm, please visit https://github.com/ZijieJin/scComm and post them on the Issues page or email me: jzj2035198@outlook.com

## Software Prerequisite

- R (tested on 4.2.1)
- R package: stringr, entropy, pracma, RcppML, NMF, matrixStats, progress

## System Requirement

- 32 GB memory or more 

- 4 CPU cores or more

## Data Requirement

- Single-cell RNA-seq gene expression data (10X or Smart-seq) with cell annotations

- (optional) custom Ligand-receptor database

- (optional) custom downstream regulation database

## Usage

`source(scComm.R)`

`res = scComm(data, anno)`

After running `scComm()`, user can run `res = FindLRscoreGivenCells(scCommRes, celllist, lr_database)` to get the detailed CCI score with each L-R pair between cell type user defined. 

`scComm()` returns a list containing weights, ccires (a list), GroupCCC, cellanno, expr.

weights contains weight 1, 2, 3, and totalweight(the product of weight 1,2,3).

ccires is a list containing result below: 
- cciscore: a cell-by-cell matrix of aggregating cell-cell interaction score of each cell pair
- lrscore: cell type-cell type interaction score with each L-R pair
- sender_lr: an LR pair-by-cell matrix of aggregating cell-cell interaction score where the element (i,j) indicates the sum of interaction score of jth cell as sender cell to all other cells with L-R pair i.
- receptor_lr: an LR pair-by-cell matrix of aggregating cell-cell interaction score where the element (i,j) indicates the sum of interaction score of all other cells to jth cell as receptor cell with L-R pair i.

GroupCCC is a matrix indicating the CCI scores between one cell type and another.

cellanno is a vector indicating the annotation of each input, which is given by user.

expr is a gene-by-cell matrix indicating the expression

res$ccires$cciscore can be used to analyze the interaction patterns (aggregate by cell type to get cell-by-cell type matrix)

GroupCCC provides cell type-level cell-cell communication scores.

`FindLRscoreGivenCells()` returns a list containing lrscore and lrscoresd, which are LR pair-by-cell group matrix indicating the mean and standard variation of CCI score of each L-R pair and cell interaction group, respectively. 


## Commercial use

For non-academic use, please email Dr. Jin (jzj2035198@outlook.com) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
