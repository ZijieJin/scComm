library(entropy)
library(Seurat)
library(pracma)
library(RcppML)
library(NMF)
library(matrixStats)
library(progress)
library(stringr)

scCommInit <- function(expr, anno, lr_database = "./data/scriabin_LR_OmniPath.txt", tf_database = './data/dorothea.rds'){
    lr_network <<- read.table(lr_database, sep = "\t")
    tfs <<- readRDS(tf_database)
    colnames(lr_network) <- c("from", "to")
    colnames(tfs)[3:4] <- c("from", "to")
    LRpairs_str <<- str_c(lr_network$from, lr_network$to, sep = "---")
    anno <- as.character(anno)
    expr <- as.matrix(expr)
    this_data_groups <<- unique(anno)
    genes <- row.names(expr)
    meanexpr_full <<- colMeans(t(expr))
    goodgene <- rep(F, dim(lr_network)[1])
    for (l in 1:length(goodgene)) {
        if (lr_network$from[l] %in% genes & lr_network$to[l] %in% genes) {
            if (meanexpr_full[lr_network$from[l]] > 0 & meanexpr_full[lr_network$to[l]] > 0) goodgene[l] <- T
        }
    }
    lr_network <<- lr_network[goodgene, ]
    LRpairs_str <<- LRpairs_str[goodgene]
    goodgene <- rep(F, dim(tfs)[1])
    for (l in 1:length(goodgene)) {
        if (tfs$from[l] %in% genes & tfs$to[l] %in% genes) {
            goodgene[l] <- T
        }
    }
    tfs <<- tfs[goodgene & tfs$is_inhibition == 0, ]
    cellanno <<- data.frame(cellname = colnames(expr), anno = anno)
}

scComm <- function(expr, anno, lr_database = "./data/scriabin_LR_OmniPath.txt", tf_database = './data/dorothea.rds') {
    scCommInit(expr, anno, lr_database, tf_database)
    print("##### Evaluating Weights #####")
    weights <- CalcWeights(expr, anno)
    print("##### Evaluating Cell-Cell Communication Score #####")
    ccires <- CalcCCI(expr, anno, weights)
    print("##### Evaluating Group Communication #####")
    GroupCCC <- CalcGroupCCI(expr, anno, ccires$cciscore)
    GroupCCC_res <- FindHighCommGroup_Permutation_Simple(expr, anno, GroupCCC, ccires$cciscore)
    GroupCCC_all <- MakeAllGroupSig()
    CCC_sig <- FindHighCommCell(expr, anno, ccires$cciscore, GroupCCC_res$GroupCCC_sig)
    LR_sig <- FindHighCommLR(expr, anno, GroupCCC_res$GroupCCC_sig, ccires$lrscore, ccires$fakegroup_LR_dist)
    CCC_sig_all <- FindHighCommCell(expr, anno, ccires$cciscore, GroupCCC_all)
    LR_sig_all <- FindHighCommLR(expr, anno, GroupCCC_all, ccires$lrscore, ccires$fakegroup_LR_dist)
    print("##### Result Saved #####")
    return(list(
        weights = weights, ccires = ccires, GroupCCC = GroupCCC, GroupCCC_res = GroupCCC_res,
        LR_sig = LR_sig, LR_sig_all = LR_sig_all, CCC_sig = CCC_sig,
        CCC_sig_all = CCC_sig_all, cellanno = cellanno, expr = as.sparse(expr)
    ))
}

CalcWeights <- function(expr, anno, set_weight = F) {
    this_data_groups <- unique(anno)
    weight1 <- rep(0, dim(lr_network)[1])
    weight2 <- weight3 <- totalweight <- array(0, c(length(this_data_groups), length(this_data_groups), dim(lr_network)[1]))
    rownames(weight2) <- this_data_groups
    colnames(weight2) <- this_data_groups
    rownames(weight3) <- this_data_groups
    colnames(weight3) <- this_data_groups
    if (set_weight == T) {
        weight1 <- rep(1, dim(lr_network)[1])
        weight2 <- weight3 <- totalweight <- array(1, c(length(this_data_groups), length(this_data_groups), dim(lr_network)[1]))
        return(list(totalweight = totalweight, weight1 = weight1, weight2 = weight2, weight3 = weight3))
    }
    # Calc the entropy of l and r
    meanexpr <- matrix(0, nrow = nrow(expr), ncol = length(this_data_groups))
    rownames(meanexpr) <- rownames(expr)
    for (k in 1:length(this_data_groups)) {
        if (sum(anno == this_data_groups[k]) >= 2) {
            meanexpr[, k] <- colMeans(t(expr[, anno == this_data_groups[k]]))
        } else {
            meanexpr[, k] <- expr[, anno == this_data_groups[k]]
        }
    }
    print("##### Calculating Weight 1")
    progress <- progress_bar$new(total = dim(lr_network)[1], clear = F)
    for (l in 1:dim(lr_network)[1]) {
        groupexpr_l <- meanexpr[lr_network$from[l], ] + 0.001
        groupexpr_r <- meanexpr[lr_network$to[l], ] + 0.001
        groupexpr_l[is.nan(groupexpr_l)] <- 0.001
        groupexpr_r[is.nan(groupexpr_r)] <- 0.001
        ent_l <- entropy::entropy(groupexpr_l)
        ent_r <- entropy::entropy(groupexpr_r)
        weight1[l] <- 1 / (1 + ent_l + ent_r)
        progress$tick()
    }
    weight1[is.na(weight1)] <- 0
    weight1 <- weight1 - min(weight1[weight1 > 0])
    weight1[weight1 < 0] <- 0
    weight1[weight1 > 1] <- 1

    # Calc the competition score
    print("##### Calculating Weight 2")
    progress <- progress_bar$new(total = dim(lr_network)[1], clear = F)
    eta <- 0.1
    for (l in 1:dim(lr_network)[1]) {
        alll <- which(lr_network$to == lr_network$to[l])
        allr <- which(lr_network$from == lr_network$from[l])
        w1mat <- w2mat <- matrix(0, nrow = length(this_data_groups), ncol = length(this_data_groups))
        for (i in 1:length(this_data_groups)) {
            maxl <- max(meanexpr[lr_network$from[alll], i])
            w1mat[i, ] <- (exp(eta * meanexpr[lr_network$from[l], i]) - 1) / (exp(eta * maxl) - 1)
            maxr <- max(meanexpr[lr_network$to[allr], i])
            w2mat[, i] <- (exp(eta * meanexpr[lr_network$to[l], i]) - 1) / (exp(eta * maxr) - 1)
        }
        weight2[, , l] <- sqrt(w1mat^2 + w2mat^2) / sqrt(2)
        totalweight[, , l] <- weight1[l] * weight2[, , l]
        progress$tick()
    }
    weight2[is.nan(weight2)] <- 0
    totalweight[is.nan(totalweight)] <- 0

    # Calc the downstream score
    print("##### Calculating Weight 3")
    progress <- progress_bar$new(total = dim(tfs)[1], clear = F)
    tfcor <- matrix(0, nrow = dim(tfs)[1], ncol = length(this_data_groups))
    for (l in 1:dim(tfs)[1]) {
        if (meanexpr_full[tfs$from[l]] == 0 | meanexpr_full[tfs$to[l]] == 0) next()
        thisv1 <- expr[tfs$from[l], ]
        thisv2 <- expr[tfs$to[l], ]
        for (k in 1:length(this_data_groups)) {
            tfcor[l, k] <- cor(thisv1[anno == this_data_groups[k]], thisv2[anno == this_data_groups[k]])
        }
        progress$tick()
    }
    tfcor[is.na(tfcor)] <- 0

    progress <- progress_bar$new(total = dim(lr_network)[1], clear = F)
    for (l in 1:dim(lr_network)[1]) {
        lcor <- rcor <- rep(0, length(this_data_groups))
        for (k in 1:length(this_data_groups)) {
            lcor[k] <- mean(tfcor[which(tfs$to == lr_network$from[l]), k])
            rcor[k] <- mean(tfcor[which(tfs$to == lr_network$to[l]), k])
            rcor[is.nan(rcor)] <- 0
            lcor[is.nan(lcor)] <- 0
        }
        weight3[, , l] <- (matrix(rep(lcor, length(this_data_groups)), ncol = length(this_data_groups)) +
            t(matrix(rep(rcor, length(this_data_groups)), ncol = length(this_data_groups)))) / 2
        progress$tick()
    }
    weight3 <- sigmoid(weight3, a = 10)
    totalweight <- totalweight * weight3
    return(list(totalweight = totalweight, weight1 = weight1, weight2 = weight2, weight3 = weight3))
}

CalcCCI <- function(expr, anno, weights) {
    fakeanno <- sample(anno, length(anno))
    this_data_groups <- unique(anno)
    totalweight <- weights$totalweight
    weight1 <- weights$weight1
    cciscore <- matrix(0, nrow = dim(expr)[2], ncol = dim(expr)[2])
    lrscore <- fakelrscore <- array(0, c(length(this_data_groups), length(this_data_groups), dim(lr_network)[1]))
    sender_lr <- receptor_lr <- matrix(0, nrow = dim(lr_network)[1], ncol = dim(expr)[2])
    fakegroup_LR_dist <- matrix(0, nrow = dim(lr_network)[1], ncol = 2)
    rownames(cciscore) <- colnames(cciscore) <- cellanno$cellname
    colnames(sender_lr) <- colnames(receptor_lr) <- cellanno$cellname
    rownames(sender_lr) <- rownames(receptor_lr) <- LRpairs_str
    dimnames(lrscore)[[1]] <- dimnames(lrscore)[[2]] <- this_data_groups
    dimnames(lrscore)[[3]] <- str_c(lr_network$from, lr_network$to, sep='---')
    progress <- progress_bar$new(total = length(lr_network$from), clear = F)
    exprsd = rowSds(as.matrix(expr))
    names(exprsd) = rownames(expr)
    for (l in 1:length(lr_network$from)) {
        progress$tick()
        thissd <- sqrt(exprsd[lr_network$from[l]] * exprsd[lr_network$to[l]])
        if (thissd < 0.1) {
            next() 
        }
        exprlij <- sqrt(expr[lr_network$from[l], ] %*% t(expr[lr_network$to[l], ]))
        exprlij_nonzero <- exprlij[exprlij > 0]
        nonzeroprop <- length(exprlij_nonzero) / length(exprlij)
        if (nonzeroprop <= 0.0001) {
            next()
        }
        thismean <- mean(exprlij_nonzero) * nonzeroprop
        thissd <- sqrt((sd(exprlij_nonzero)^2 + (1 - nonzeroprop) * mean(exprlij_nonzero)^2) * nonzeroprop)
        exprlij <- (exprlij - thismean) / thissd
        exprlij[exprlij < 0] <- 0
        exprlij[exprlij > 5] <- 5
        cciaddmat <- exprlij * totalweight[anno, anno, l]
        fakecciaddmat <- exprlij * totalweight[fakeanno, fakeanno, l]
        sender_lr[l, ] <- colSums(t(cciaddmat))
        receptor_lr[l, ] <- colSums(cciaddmat)
        cciscore <- cciscore + cciaddmat
        for (i in 1:length(this_data_groups)) {
            cccmat_s <- cciaddmat[which(anno == this_data_groups[i]), ]
            fakecccmat_s <- fakecciaddmat[which(fakeanno == this_data_groups[i]), ]
            for (j in 1:length(this_data_groups)) {
                if (sum(anno == this_data_groups[i]) == 1) {
                    subcccmat <- cccmat_s[which(anno == this_data_groups[j])]
                    fakesubcccmat <- fakecccmat_s[which(fakeanno == this_data_groups[j])]
                } else {
                    subcccmat <- cccmat_s[, which(anno == this_data_groups[j])]
                    fakesubcccmat <- fakecccmat_s[, which(fakeanno == this_data_groups[j])]
                }
                subcccmat <- subcccmat[subcccmat > 0]
                fakesubcccmat <- fakesubcccmat[fakesubcccmat > 0]
                if (length(subcccmat) > 0) {
                    lrscore[i, j, l] <- mean(subcccmat)
                }
                if (length(fakesubcccmat) > 0) {
                    fakelrscore[i, j, l] <- mean(fakesubcccmat)
                }
            }
        }
        nullmean <- mean(fakelrscore[, , l])
        nullsd <- sd(fakelrscore[, , l])
        fakegroup_LR_dist[l, ] <- c(nullmean, nullsd)
    }
    return(list(cciscore = cciscore, lrscore = lrscore, sender_lr = sender_lr, receptor_lr = receptor_lr, fakegroup_LR_dist = fakegroup_LR_dist))
}

CalcGroupCCI <- function(expr, anno, cciscore) {
    GroupCCC <- matrix(0, nrow = length(this_data_groups), ncol = length(this_data_groups))
    rownames(GroupCCC) <- this_data_groups
    colnames(GroupCCC) <- this_data_groups
    for (i in 1:nrow(GroupCCC)) {
        for (j in 1:ncol(GroupCCC)) {
            subcccmat <- cciscore[which(anno == this_data_groups[i]), which(anno == this_data_groups[j])]
            if (length(subcccmat) > 0) GroupCCC[i, j] <- mean(as.vector(subcccmat))
        }
    }
    rownames(GroupCCC) <- colnames(GroupCCC) <- this_data_groups
    return(GroupCCC)
}

FindHighCommGroup_Permutation_Simple <- function(expr, anno, GroupCCC, cciscore) {
    GroupCCC_sig <- data.frame(g1 = 0, g2 = 0)
    nullgroupccc <- c()
    for (i in 1:5) {
        fakeanno <- sample(anno, length(anno))
        fakeGroupCCC <- CalcGroupCCI(expr, fakeanno, cciscore)
        nullgroupccc <- c(nullgroupccc, as.vector(fakeGroupCCC))
    }
    fakemean <- mean(nullgroupccc)
    fakesd <- sd(nullgroupccc)
    GroupCCC_Zscore <- (GroupCCC - fakemean) / fakesd
    GroupCCC_pv <- 1 - pnorm(GroupCCC_Zscore)
    GroupCCC_padj <- matrix(p.adjust(GroupCCC_pv), nrow = nrow(GroupCCC_pv))
    for (i in 1:nrow(GroupCCC)) {
        for (j in 1:ncol(GroupCCC)) {
            if (sum(is.na(GroupCCC_padj)) > 0) {
                next()
            }
            if (GroupCCC_padj[i, j] < 0.05) {
                GroupCCC_sig <- rbind(GroupCCC_sig, data.frame(g1 = this_data_groups[i], g2 = this_data_groups[j]))
            }
        }
    }
    GroupCCC_sig <- GroupCCC_sig[-1, ]
    return(list(GroupCCC_sig = GroupCCC_sig, GroupCCC_Zscore = GroupCCC_Zscore, GroupCCC_padj = GroupCCC_padj))
}

FindHighCommCell <- function(expr, anno, cciscore, GroupCCC_sig) {
    cciscore_sig <- matrix(0, nrow = nrow(cciscore), ncol = ncol(cciscore))
    cciscore_v <- as.vector(cciscore)
    quant <- quantile(cciscore_v, c(0.05, 0.95))
    cciscore_sub <- cciscore_v[cciscore_v >= quant[1] & cciscore_v <= quant[2]]
    cciscoremean <- mean(cciscore_sub)
    cciscoresd <- sd(cciscore_sub)
    cciscore_pv <- 1 - pnorm(cciscore, cciscoremean, cciscoresd)
    cciscore_padj <- matrix(p.adjust(cciscore_pv), nrow = nrow(cciscore_pv))
    for (i in 1:nrow(GroupCCC_sig)) {
        cciscore_sig[which(anno == GroupCCC_sig[i, 1]), which(anno == GroupCCC_sig[i, 2])] <- 1
    }
    cciscore_sig <- cciscore_sig * (cciscore_padj < 0.05)
    rownames(cciscore_sig) <- colnames(cciscore_sig) <- cellanno$cellname
    return(cciscore_sig)
}

FindHighCommLR <- function(expr, anno, GroupCCC_sig, lrscore, group_LR_dist) {
    LR_sig <- data.frame(g1 = 0, g2 = 0, l = 0, r = 0)
    if (nrow(GroupCCC_sig) == 0) {
        LR_sig <- LR_sig[-1, ]
        return(LR_sig)
    }
    lrpv <- array(1, c(length(this_data_groups), length(this_data_groups), dim(lr_network)[1]))
    for (l in 1:dim(lr_network)[1]) {
        if (group_LR_dist[l, 2] < 0.001) {
            next()
        }
        for (i in 1:length(this_data_groups)) {
            for (j in 1:length(this_data_groups)) {
                lrpv[i, j, l] <- 1 - pnorm(lrscore[i, j, l], mean = group_LR_dist[l, 1], sd = group_LR_dist[l, 2])
            }
        }
    }
    lrpv[, , l] <- matrix(p.adjust(lrpv[, , l]), nrow = length(this_data_groups))
    for (k in 1:nrow(GroupCCC_sig)) {
        for (l in 1:nrow(lr_network)) {
            if (lrpv[which(this_data_groups == GroupCCC_sig$g1[k]), which(this_data_groups == GroupCCC_sig$g2[k]), l] < 0.05) {
                LR_sig <- rbind(LR_sig, data.frame(
                    g1 = GroupCCC_sig$g1[k], g2 = GroupCCC_sig$g2[k],
                    l = lr_network$from[l], r = lr_network$to[l]
                ))
            }
        }
    }
    LR_sig <- LR_sig[-1, ]
    return(LR_sig)
}

MakeAllGroupSig <- function() {
    GroupCCC_all <- data.frame(g1 = 0, g2 = 0)
    for (i in this_data_groups) {
        for (j in this_data_groups) {
            GroupCCC_all <- rbind(GroupCCC_all, data.frame(g1 = i, g2 = j))
        }
    }
    GroupCCC_all <- GroupCCC_all[-1, ]
    return(GroupCCC_all)
}

FindLRscoreGivenCells <- function(expr, scCommRes, celllist, lr_database = "scriabin_LR_OmniPath.txt"){
    scCommInit(expr, anno, lr_database)
    totalweight <- scCommRes$weights$totalweight
    
    celllist = celllist[celllist$CellName %in% cellanno$cellname,]
    senderCells = unique(celllist$CellName[celllist$Tag == 'Sender'])
    receptorCells = unique(celllist$CellName[celllist$Tag == 'Receptor'])
    Groups = unique(celllist$Group)
    if (length(Groups) == 1){
        print('Only One Group! Cannot compare! Exit....')
        return()
    }
    if (length(senderCells) == 0){
      print('No Sender Cells! Exit....')
      return()
    }
    if (length(receptorCells) == 0){
      print('No Receptor Cells! Exit....')
      return()
    }
    
    cciscore <- matrix(0, nrow = length(senderCells), ncol = length(receptorCells))
    lrscore <- lrscoresd <- matrix(0, nrow = dim(lr_network)[1], ncol = length(Groups))
    rownames(cciscore) <- senderCells
    colnames(cciscore) <- receptorCells
    colnames(lrscore) <- colnames(lrscoresd) <- Groups
    rownames(lrscore) <- rownames(lrscoresd) <- str_c(lr_network$from, lr_network$to, sep='---')
    
    expr = expr[, c(senderCells, receptorCells)]
    annorow = cellanno$anno[cellanno$cellname %in% senderCells]
    annocol = cellanno$anno[cellanno$cellname %in% receptorCells]
    
    progress <- progress_bar$new(total = length(lr_network$from), clear = F)
    for (l in 1:length(lr_network$from)) {
        progress$tick()
        exprlij <- sqrt(expr[lr_network$from[l], ] %*% t(expr[lr_network$to[l], ]))
        exprlij_nonzero <- exprlij[exprlij > 0]
        if (length(exprlij_nonzero) == 0) next()
        nonzeroprop <- length(exprlij_nonzero) / length(exprlij)
        thismean <- mean(exprlij_nonzero) * nonzeroprop
        thissd <- sqrt((sd(exprlij_nonzero)^2 + (1 - nonzeroprop) * mean(exprlij_nonzero)^2) * nonzeroprop)
        exprlij <- (exprlij - thismean) / thissd
        exprlij[exprlij < 0] <- 0
        exprlij[exprlij > 5] <- 5
        exprlij <- exprlij * totalweight[c(annorow, annocol), c(annorow, annocol), l]
        colnames(exprlij) = rownames(exprlij) = c(senderCells, receptorCells)
        for (g in Groups){
            thisSender = celllist$CellName[celllist$Group == g & celllist$Tag == 'Sender']
            thisReceptor = celllist$CellName[celllist$Group == g & celllist$Tag == 'Receptor']
            if (length(thisSender) == 0 | length(thisReceptor) == 0) next()
            subexprlij = exprlij[thisSender, thisReceptor]
            lrscore[l, g] = mean(subexprlij)
            lrscoresd[l, g] = sd(subexprlij)
        }
    }
    return(list(lrscore = lrscore, lrscoresd = lrscoresd))
}

