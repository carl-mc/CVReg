

#' Estimate covoting coefficients for categorical cleavage variables from multinomial model
#'
#' @param data Individual-level data.
#' @param out Main outcome of interest, must be either integer, string, or factor
#' @param cov Character vector of covariates. First entry must denote the main
#'    group variable of interest for which covoting coefficient is computed.
#'    The first entry must be either an integer, string, or factor.
#' @param method \code{"ols"} or \code{"glm"}. Defaults to OLS, thus using \link[fixest]{feols}
#'   for estimation. If \code{"glm"}, it estimates a multinomial logistic regression
#'   via  \link[nnet]{multinom}.
#' @param ... arguments passed to  \link[fixest]{feols} or  \link[nnet]{multinom}.
#'
#' @returns A vector with the main covoting coefficient \code{"beta"}, as well as the
#'   average predicted covoting probabilities within and between groups,
#'   based on mean values on all covariates.
#' @export
#'
#' @import fixest
#' @import nnet
#' @importFrom stats formula
#' @importFrom stats na.omit
#' @importFrom stats predict
#' @examples
#' # Load some data
#' data(synthvote)
#'
#' # Covoting coefficient from multinomial choice model
#' est_covote_indlevel(data = synthvote, out = "prtvt", cov = c("group", "income"), method = "ols")
est_covote_indlevel <- function(data, out, cov, method = "ols", ...){
    # To DataFrame if data.table
    if("data.table" %in% class(data)){
        data <- data.frame(data)
    }

    # Check outcome
    if(!is.factor(data[, out]) & !is.integer(data[, out]) & !is.character(data[, out])){
        stop("out must be either integer, character, or factor. ")
    }

    # Check main cleavage
    if(!is.factor(data[, cov[1]]) & !is.integer(data[, cov[1]]) & !is.character(data[, cov[1]])){
        stop("The first entry in cov must be either integer, character, or factor. ")
    }

    ## Subset data
    data <- na.omit(data[, c(out, cov)])

    ## Outcomes
    out.val <- sort(unique(data[, out]))

    ## Prep control variables
    data.cov <- data.frame(matrix(NA, nrow = nrow(data), ncol = 0))
    for(v in cov[-1]){
        if(is.character(data[,v]) | is.factor(data[,v])){
            vals <- na.omit(unique(data[,v]))[-1]
            for(i in vals){
                data.cov <- cbind(data.cov,
                                  as.numeric(data[,v] == i))
                colnames(data.cov)[ncol(data.cov)] <- paste0(v, ".", which(vals == i))
                data.cov[,ncol(data.cov)] <- data.cov[,ncol(data.cov)] - mean(data.cov[,ncol(data.cov)])
            }
        } else {
            data.cov <- cbind(data.cov,
                              data[,v])
            colnames(data.cov)[ncol(data.cov)] <- paste0(v)
            data.cov[,ncol(data.cov)] <- data.cov[,ncol(data.cov)] - mean(data.cov[,ncol(data.cov)])
        }
    }

    # Data for regression
    data.m <- cbind(data[,!colnames(data) %in% colnames(data.cov)],
                    data.cov)

    # Make group factor
    data.m[, cov[1]] <- as.factor(data.m[, cov[1]])

    # Data for predicted probs
    grps <- sort(unique(data[,cov[1]]))
    pred.df <- data.m[seq_along(grps),]
    pred.df[, colnames(data.cov)] <- 0
    pred.df[, cov[1]] <- factor(grps,
                                levels = attr(data.m[, cov[1]] , "levels"),
                                ordered = T)

    # Multinomial logit // NNET
    if(method == "glm"){
        ## Estimate
        mlog <- nnet::multinom(formula(paste("factor(", out, ") ~ 0 + ",
                                             paste(c(paste0("factor(", cov[1], ")"), colnames(data.cov)), collapse = " + "))),
                               data = data.m, ...)

        ## Voting probabilities
        pmat <- predict(mlog, newdata = pred.df, type = "probs")
        if(length(unique(data.m[,out])) <= 2){
            pmat <- cbind(pmat, 1-pmat)
        }
        rownames(pmat) <- grps

    } else if(method == "ols"){
        ## Estimate
        m_ls <- fixest::feols(formula(paste("c(", paste(paste0("I(", out, "=='", out.val, "')"), collapse = ","), ") ~ 0 + ",
                                            paste(c(paste0("factor(", cov[1], ")"), colnames(data.cov)), collapse = " + "))),
                              data = data.m, ...)



        ## Voting probabilities
        pmat <- do.call(cbind, lapply(m_ls, function(m){
            predict(m, newdata = pred.df, type = "response")
        }))
        rownames(pmat) <- grps

    } else {
        stop("Invalid method")
    }

    ## Group comparisons
    covote.mat <- do.call(rbind, lapply(grps, function(i){
        sapply(grps, function(j){
            sum(pmat[i,] * pmat[j,])
        })
    }))

    ## Weights
    weights <- matrix(table(data[,cov[1]])[grps], ncol = 1) %*% t(table(data[,cov[1]])[grps])
    diag(weights) <- diag(weights) - table(data[,cov[1]])[grps]

    ## Covoting among (non)group members
    coeth <- sum(diag(covote.mat) * diag(weights) / sum(diag(weights)))
    nocoeth <- sum(covote.mat[upper.tri(covote.mat)] * weights[upper.tri(weights)] / sum(weights[upper.tri(weights)]))
    beta <- coeth - nocoeth

    ## Return
    return(c(beta = beta, covoting_within_grp = coeth, covoting_between_grp = nocoeth))
}




