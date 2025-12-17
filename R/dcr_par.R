
##############################################################
### Compute cluster–robust standard errors for dyadic data ###
##############################################################

#' Fast Export function for parallel processing
#'
#' @param cl
#' @param varlist
#' @param envir
clusterExport_fast <- function (cl, varlist, envir = .GlobalEnv){
    require(foreach)
    registerDoParallel(cl)
    tmp.file <- tempfile()
    save(list = varlist, file = tmp.file, envir = envir)
    clusterExport(cl, list("tmp.file"), envir = environment())
    load <- foreach(i = seq_along(cl), .noexport = c("tmp.file")) %dopar% {
        load(tmp.file, envir = .GlobalEnv)
        TRUE
    }
    clusterEvalQ(cl = cl, expr = {
        rm(tmp.file)
        TRUE
    })
    unlink(tmp.file)
    return(TRUE)
}


#' Helper function to create dyad IDs
#'
#' @param data
#' @param dyad_mem1
#' @param dyad_mem2
#' @param directed
#'
#' @returns
#'
#' @examples
create_dyadid <- function (data, dyad_mem1, dyad_mem2, directed = F){
    dvars <- c(dyad_mem1, dyad_mem2)
    sorteddyads <- as.matrix(data[, ..dvars])
    if (directed != T) {
        sorteddyads[sorteddyads[,1] > sorteddyads[,2], c(1,2)] <-
            sorteddyads[sorteddyads[,1] > sorteddyads[,2], c(2,1)]
    }
    return(as.integer(as.factor(paste0(sorteddyads[,1], ".", sorteddyads[,2]))))
}

#' This function computes DCR variance estimates and standard errors
#' (DCRSEs) for model parameters.
#'
#' @param model A model object.
#' @param dyad_mem1 A string for the name of first dyad member identifier variable.
#' @param dyad_mem2 A string for the name of second dyad member identifier variable.
#' @param data A `data.frame` object containing dyadic data.
#' @param posdef A logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite.
#' @param ncore
#' @param type
#' @param dofcorr A logical value indicating whether or not to apply a small sample correction to the final DCRSE estimates. This correction is equivalent to multiplying DCRSEs by sqrt(N/(N-1)), where N is the number of unique dyad members in the analytic sample. Correspondingly, when computing p-values, the test statistic should be compared to a t-distribution with DOF = N - 1.
#'
#' @return A list containing DCRSEs and the DCR variance-covariance matrix for model parameters, as well as the number of unique dyad members in the analytic sample.
#' @export
dcr_par <- function(model, dyad_mem1, dyad_mem2, data, ncore = 1, type = "sandwich", posdef = FALSE, dofcorr = FALSE) {

    # Make Dyad-IDs
    data[, dcr_dyad_id := create_dyadid(data, dyad_mem1, dyad_mem2, directed = F)]

    # More house keeping
    gp.tag <- -99
    unique.dyad.mem <- sort(unique(c(data[[dyad_mem1]], data[[dyad_mem2]])))
    N_dyad <- length(unique.dyad.mem)
    d.vars <- c(dyad_mem1, dyad_mem2)
    dyad.by.obs <- data[, ..d.vars]

    # Cluster
    if(ncore > 1){
        require(foreach)
        cl <- parallel::makeCluster(getOption("cl.cores", ncore))
        registerDoParallel(cl)
        exp.vec <- c("unique.dyad.mem", "dyad.by.obs", "gp.tag", "model", "dyad_mem1", "dyad_mem2", "type")
        clusterExport_fast(cl, varlist = exp.vec, envir = environment())
    }

    ## Compute
    cov.mat.ls <- foreach(i = 1:N_dyad, .noexport = exp.vec, .packages = c("data.table")) %dopar% {
        dyad.mem.i <- unique.dyad.mem[i]
        dyad.with.i <-  (unlist(dyad.by.obs[,..dyad_mem1]) == dyad.mem.i | unlist(dyad.by.obs[,..dyad_mem2]) == dyad.mem.i)
        dyad.category <- dyad.with.i * gp.tag + (1 - dyad.with.i) * 1:nrow(dyad.by.obs)

        if(type == "sandwich"){
            sandwich::vcovCL(
                model,
                dyad.category,
                type = "HC0", multi0 = TRUE, cadjust = FALSE, fix = FALSE)
        } else {
            multiwayvcov::cluster.vcov(model, dyad.category, df_correction = FALSE)
        }
    }

    ## Stop cluster
    if(ncore > 1){
        stopCluster(cl)
    }


    ## Sum
    cov.mat.sum <- Reduce('+', cov.mat.ls)

    # substract HC and repeated variance estimator for repeated dyads
    if(type == "sandwich"){
        cov.mat.sum.intermed <- cov.mat.sum - sandwich::vcovCL(
            model,
            data[["dcr_dyad_id"]],
            type = "HC0", multi0 = TRUE, cadjust = FALSE, fix = FALSE)
        V.hat <- cov.mat.sum.intermed - (N_dyad - 2) *
            sandwich::vcovCL(model, 1:nrow(data), type = "HC0",
                             multi0 = TRUE, cadjust = FALSE, fix = FALSE)
    } else {
        cov.mat.sum.intermed <- cov.mat.sum - multiwayvcov::cluster.vcov(model, data[["dcr_dyad_id"]], df_correction = FALSE)
        V.hat <- cov.mat.sum.intermed - (N_dyad - 2) * multiwayvcov::cluster.vcov(model, 1:nrow(data), df_correction = FALSE)
    }


    # force posdef
    param.names <- colnames(V.hat)
    if (sum(diag(V.hat) < 0) >= 1 & posdef == TRUE) {
        decomp <- eigen(V.hat, symmetric = TRUE)
        pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
        V.hat <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
    }

    # DOF correction
    if (dofcorr == TRUE) {
        coef.vcov <- (N_dyad / (N_dyad - 1)) * V.hat
    } else {
        coef.vcov <- V.hat
    }

    # Compute SEs
    coef.se <- sqrt(diag(coef.vcov))

    # Finalize output
    outputlst <- list(coef.se, coef.vcov, N_dyad)
    names(outputlst) <- c("dcrse", "dcrvcov", "N_udm")
    names(outputlst$dcrse) <- param.names

    # Return
    return(outputlst)

}



###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
# https://github.com/jscarlson/dcr/blob/master/R/dyad_funcs.R
