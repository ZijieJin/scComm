# scComm

scComm is a computational pipeline for inferring cell-cell communication (CCC) at single-cell resolution. The main function of scComm includes:
- Calculate CCC scores between any single cell pairs to reveal the strength of intercellular communications
- Calculate CCC scores between two cell types like CellChat and CellphoneDB
- Report CCC scores between any custom cell groups without the need for rerun.

If you have any questions related to scComm, please visit https://github.com/ZijieJin/scComm and post them on the Issues page or email me: jzj2035198@outlook.com

## Software Prerequisite
scComm works on R and Python platform with any OS.
- R (tested on 4.2.3)
- R package: stringr, entropy, pracma, RcppML, NMF, matrixStats, progress
- Python (tested on 3.12)
- Python modules: numpy (1.26.4), pandas (2.2.3), torch (2.5.0), scikit-learn (1.5.2)

## System Requirement

To run scComm properly, Your computer should have:

- 32 GB memory or more 


## Data Requirement

scComm requires three inputs:

- Single-cell RNA-seq gene expression data with cell annotations. The gene expression matrix is a gene-by-cell matrix, and cell annotations are provided by a tow-column csv file, including cell names and cell types.

- (optional) custom Ligand-receptor database. User can also provide other L-R databases, which is a two-column csv file indicating ligands and receptors. 

- (optional) downstream regulation database. User can also provide other downstream regulation databaases. The file format should be the same as  `dorothea.rds`



## Usage

First, source the functions of scComm using:

`source(scComm.R)`

Then, run scComm using all-in-one function:

`res = scComm(data, anno)`

To identify significant interating cell type pairs, run the supervised contrastive learning by two steps:

Step 1: run the commend in R `MakeAugmentedData(data, res)` to generate training datasets for the  contrastive learning;

Step 2: run the commend in Python `python codes/ContrastiveLearning.py PositiveData.csv NegativeData.csv Fulldata.csv`. The result is stored in 'fulldata_prediction.csv'

## Result

`scComm()` returns a list containing weights, ccires (a list), GroupCCC, cellanno, expr.

*weights* contains weight 1, 2, 3, and totalweight(the product of weight 1,2,3).

*ccires* is a list containing main result below: 
- *cciscore*: a cell-by-cell matrix of aggregating cell-cell interaction score of each cell pair (important)
- *lrscore*: cell type-cell type interaction score with each L-R pair, which can be used to analyze the interacting L-R pairs with given cells (important)

*GroupCCC* is a matrix indicating the CCC scores between two cell types. It provides cell type-level cell-cell communication scores.

*cellanno* is a vector indicating the annotation of each input, which is given by user.

*expr* is a gene-by-cell matrix indicating the expression

`res$ccires$cciscore` can be used to analyze the interaction patterns (aggregate by cell type to get cell-by-cell type matrix)

After running `scComm()`, user can run `res2 = FindLRscoreGivenCells(res, celllist, lr_database)` to get the detailed CCI score with each L-R pair between cell groups user defined. Here, `res` is the result returned by `scComm()`, celllist is a data.frame (see the format in testdata)

`FindLRscoreGivenCells()` returns a list containing lrscore and lrscoresd, which are LR pair-by-cell group matrix indicating the mean and standard variation of CCI score of each L-R pair and cell interaction group, respectively. 

`res2$lrscore` can be used to analyze the interacting L-R pairs with given cells.


## Sample Usage on Testdata  

Download the data in the folder `testdata/`, and run 

`source('codes/scComm.R')`

`expr = readRDS('testdata/expr.rds')`

`anno = readRDS('testdata/anno.rds')`

`res = scComm(expr, anno)`

Tested on MacStudio M1Max 10 CPUs with 32GB memory, it costs less than 5 minutes.
The output is listed at `./res2.rds`

## Commercial Use

For non-academic use, please email Dr. Jin (jzj2035198@outlook.com) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
