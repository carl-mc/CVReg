

##################################
# MAIN ESTIMATION FUNCTION
##################################

#' Estimate Covoting Regression Model
#'
#' @param fml Regression formula, following the notation used by \code{fixest}
#' @param data Dataset, either a data.table or data.frame.
#' @param fixef Character vector with fixed effects.
#' @param cluster Character vector of clusters to cluster standard errors on.
#' @param weights Variable name recording weights of observations to be used in regression.
#' @param cvr.det.var Variable name of the main explanatory variable
#' to compute conditional covoting and shared identity quantities for.
#' @param main.int.var Variable name for main explanatory variable to construct
#' cross-cleavage interactions. Defaults to NULL (no interactions added to the regression).
#' Will be expanded using c("a", "b") prefixes and can include CVReg suffix.
#' @param int.vars Variables to be interacted with each other and the main explanatory variable.
#' Will be expanded using c("a", "b") prefixes and can include CVReg suffix.
#' @param collin.tol Tolerance for collinearity in fixest estimation.
#' Increase to avoid multi-collinearity issues when including many interaction terms.
#' @param int_lasso Logical flag on whether to reduce set of interaction variables through a
#' one-step LASSO which selects at most max.int.vars interaction variables with the highest predictive
#' power over the main outcome of interest.
#' @param max.int.vars Maximum number of interaction terms selected through the LASSO procedure.
#' @param lean Return a "lean" fixest model with all chunky objects dropped.
#' @param model Defaults to "ols", but can be "glm" to estimate a logistic CVReg model.
#' @param dcr Logical flag to compute cluster-robust standard errors for network data (Samii, 2015).
#' @param ncore Number of threads to use for DCR standard error computation. Defaults to 1.
#' @param ... Additional arguments to pass to \link[fixest]{feols} and  \link[fixest]{feglm}
#'
#' @returns A \code{fixest} model object with some additional entries:
#'  \itemize{
#'      \item{\code{"Y_mean"}:} Average of outcome variable
#'      \item{\code{"N_nodes"}:} Number of unique nodes. Only returned if identifiers "a.ID", "b.ID" in the data.
#'      \item{\code{"sum.ls"}:} Results from computing conditional covoting and shared group membership quantities, if specified.
#'      \item{\code{"dcr"}:} Cluster-robust standard errors and covariance matrix for network data, if specified.
#' }
#' @export
#'
#' @importFrom glmnet cv.glmnet
#' @import fixest
#' @import data.table
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats na.omit
#'
#' @examples
#' # Load some data
#' data(synthvote)
#'
#' # Transform to pairs
#' pair.dt <- data_to_pairs(df = synthvote, by.var = c("unit"), structure = "full")
#'
#' # Process pair level data to create pairwise similarity/difference measures
#' pair.dt <- process_pairs(df = pair.dt,
#'                          vars =  c("prtvt.eq", "group.eq", "income.sim"))
#'
#' # Estimate model
#' m <- estimate_cvreg(fml = prtvt.eq ~ group.eq + income.sim, fixef = c("unit"),
#'                     cluster = c("a.ID", "b.ID"),
#'                     data = pair.dt, cvr.det.var = "group.eq")
#'
#' summary(m)
#' unlist(m$sum.ls[c(1,2)])
#' do.call(cbind, m$sum.ls[-c(1,2)])
estimate_cvreg <- function(fml,
                         data, fixef, cluster, weights = NULL,
                         cvr.det.var = NULL,
                         main.int.var = NULL, int.vars = c(),
                         collin.tol = 1e-10, int_lasso = TRUE,
                         max.int.vars = 100,
                         lean = TRUE,
                         model = "ols", dcr = FALSE, ncore = 1, ...){

    ## Check Data format
    if(!"data.table" %in% class(data)){
        data <- data.table(data)
    }

    ## Check cluster variables and drop NA
    if(!is.null(cluster)){
        if(all(paste0(c("a", "b"),".", rep(cluster, each = 2)) %in% colnames(data))){
            cluster <- paste0(c("a", "b"),".", rep(cluster, each = 2))
        }
        na.omit(data, cols = cluster)
    }

    ## Fixed effects vars
    if(length(fixef) > 0){
        if(all(paste0(c("a", "b"),".", rep(fixef, each = 2)) %in% colnames(data))){
            fixef <- paste0(c("a", "b"),".", rep(fixef, each = 2))
        }
        na.omit(data, cols = fixef)
    }

    # Weights
    if(!is.null(weights)){
        na.omit(data, cols = weiqhts)
    } else {
        weight.vec <- rep(1, nrow(data))
    }

    # Make interactions
    if(!is.null(main.int.var)){
        if(length(main.int.var) > 1){
            stop("Only one main interaction variable allowed. ")
        }

        ## Make interactions
        data <- cvr_cleav_int(df = data,
                             main.var = main.int.var,
                             int.vars = int.vars )
        add.int.vars <- colnames(data)[grepl("int__", colnames(data))]

        ## Lasso Selection
        if(int_lasso){
            print("LASSO Interaction Selection")

            ## Main covariates (fixed)
            cov <- all.vars(fml)[-1]

            ## Ridge regression for penalties
            explvars <- c(cov, add.int.vars)
            modelr <- cv.glmnet(x = as.matrix(data[, ..explvars]),
                                y = data[[all.vars(fml)[1]]],
                                alpha = 0)
            coefr <- as.matrix(coef(modelr, s = modelr$lambda.min))
            w.r <- 1/(abs(coefr[-1,]))
            w.r[1:length(cov)] <- 0

            ## Fit Lasso
            cvfit <- cv.glmnet(x = as.matrix(data[, ..explvars]),
                               y = data[[all.vars(fml)[1]]], alpha = 1, dfmax = max.int.vars + length(cov),
                               penalty.factor = w.r)

            ## Get coefs
            lasso_coef <- coef(cvfit, s = "lambda.min")

            ## Adjust add.int.vars
            add.int.vars <- rownames(lasso_coef)[as.vector(lasso_coef) != 0]
            add.int.vars <- add.int.vars[!add.int.vars %in% c(cov, "(Intercept)")]
        }

        # Adjust formula
        f.chr <- as.character(fml)[-1]
        cov.f <- strsplit(f.chr[2], "|", fixed = T)[[1]]
        cov.f[1] <- paste(cov.f[1], "+" , paste(add.int.vars, collapse ="+") )
        fml <- as.formula(paste(f.chr[1], "~", paste(cov.f, collapse = "|")))
    } else {
        add.int.vars <- c()
    }


    # Estimate
    if(model == "ols"){
        m <- feols(fml,
                   data = data,
                   cluster = cluster,
                   fixef = fixef,
                   weights = weight.vec,
                   lean = FALSE, collin.tol = collin.tol, ...)
    } else if(model == "logit"){
        m <- feglm(m.spec,
                   data = data,
                   family = "logit",
                   cluster = cluster,
                   fixef = fixef,
                   weights = weight.vec,
                   lean = FALSE, collin.tol = collin.tol,  ...)
    }


    # Post process
    if("fixest_multi" %in% class(m)){
        m <- lapply(m, function(i){
            post_process_cvr(df = data, m = i, lean = lean, weights = weights,
                             fe = fixef, cluster= cluster, cvr.det.var = cvr.det.var, dcr = dcr, weight.vec = weight.vec, ncore = ncore)
        })
    } else {
        m <- post_process_cvr(df = data, m = m, lean = lean, weights = weights,
                              fe = fixef, cluster= cluster, cvr.det.var = cvr.det.var, dcr = dcr, weight.vec = weight.vec, ncore = ncore)
    }


    # Return
    return(m)
}


#' Post-processing function for estimate_cvreg
#'
#' Helper function used to process the internal results of \code{estimate_cvreg()}.
#'
#' @param df The data used to fit model \code{m}
#' @param m Model as returned by \code{fixest} estimator
#' @param lean Logical flag to strip the model of all large objects.
#' @param weights Weight variable name
#' @param fe Names of fixed effects
#' @param cluster Names of cluster variables used for standard error clustering
#' @param cvr.det.var Name of main cleavage measure to estimate within and between group covoting
#'     and shared group membership among covoters and non-covoters for. Only works for binary measures (i.e., variable.eq)
#' @param dcr Logical flag to compute dyadic clustered standard errors (Samii, 2015)
#' @param weight.vec Vector of weights
#' @param ncore Number of CPUs used for DCR standard error clustering.
#'
#' @importFrom msm deltamethod
#' @import data.table
#'
#' @returns A modified model object based on \code{m}.
post_process_cvr <- function(df, m, lean, weights, fe, cluster, cvr.det.var, dcr, weight.vec = NULL, ncore = 1){

    # Record numbers

    ## Mean outcome
    m$Y_mean <- mean(m$fitted.values)

    ## Respondents
    if(all(c("a.ID", "b.ID") %in% colnames(df))){
        ind_vec <- unlist(df[, c("a.ID", "b.ID")])
        m$N_nodes <- length(unique(ind_vec))
    }



    # Standard errors
    if(dcr){
        print("SE clustering")

        ## Warning if weights
        if(!is.null(weights)){
            warning("Weights are not implemented with DCR.")
        }

        ## Make vars
        if(!is.null(cluster) & !all(cluster == c("a.ID", "b.ID"))){
            # Put them on the same factor levels
            fac <- unique(unlist(df[,..cluster]))
            df[, (cluster) := lapply(.SD, function(x){as.integer(factor(x, levels = fac, ordered = T))}), .SDcols = cluster]
        } else if(is.null(cluster)) {
            cluster = c("a.ID", "b.ID")
        }

        ## Run
        m$dcr <- dcr_par(model = m,
                         dyad_mem1 = cluster[1], dyad_mem2 = cluster[2],
                         data = df, ncore = ncore,
                         type = "sandwich", posdef = TRUE)
    }

    # Record group details
    if(!is.null(cvr.det.var)){
        if(!is.null(weight.vec)){
            weight.vec <- rep(1, nrow(df))
        }

        # Get correct vcv
        if(dcr){
            vcov <- m$dcr$dcrvcov
        } else {
            vcov <- m$cov.scaled
        }

        # Check subsetting
        if(!is.null(m$obs_selection) & length(m$obs_selection) > 0){
            df <- df[m$obs_selection$obsRemoved,]
            weight.vec <- weight.vec[m$obs_selection$obsRemoved]
        }

        # Covoting measures
        id <- which(names(m$coefficients) == cvr.det.var)
        grpeq.mean = sum(df[, ..cvr.det.var] * weight.vec / sum(weight.vec), na.rm = T)
        tot.cov.rate <- sum(m$fitted.values * weight.vec / sum(weight.vec))

        m$sum.ls <- list(
            covote_rate = tot.cov.rate,
            grpeq_mean = grpeq.mean
        )

        # Covariate means
        vars <- names(m$coefficients)
        if(length(fe) == 0){
            vars <- vars[-1]
        }
        nogrp.mean <- sapply(vars,
                             function(v){
                                 mean(df[[v]])
                             })
        names(nogrp.mean) <- vars
        if(length(fe) == 0){
            nogrp.mean <- c("(Intercept)" = 1,
                            nogrp.mean)
            add = 0
        } else {
            add = tot.cov.rate
        }
        all.mean <- grp.mean <- nogrp.mean
        grp.mean[cvr.det.var] <- 1
        nogrp.mean[cvr.det.var] <- 0

        stopifnot(length(grp.mean) == length(m$coefficients))

        grp.mean.mat <- matrix(grp.mean, ncol = 1) %*% grp.mean
        nogrp.mean.mat <- matrix(nogrp.mean, ncol = 1) %*% nogrp.mean


        ### Covoting within and between groups
        m$sum.ls$covoting_within_grp = c(est = sum(m$coefficients * grp.mean) + add,
                                         se = sum(grp.mean.mat * vcov)^.5)

        m$sum.ls$covoting_between_grp = c(est = sum(m$coefficients * nogrp.mean) + add,
                                          se = sum(nogrp.mean.mat * vcov)^.5)

        m$sum.ls$covoting_ratio <- c(est = m$sum.ls$covoting_within_grp[1] / m$sum.ls$covoting_between_grp[1],
                                     se = msm::deltamethod(g = as.formula(paste("~", "(", add, "+", paste(paste0(grp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ") /",
                                                                                "(", add, "+", paste(paste0(nogrp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ") ")),
                                                           mean = c(m$coefficients),  cov = vcov, ses = TRUE))

        ## Shared group membership among covoters and non-covoters
        m$sum.ls$sharedgrp_covoters <- c(est = grpeq.mean *
                                             m$sum.ls$covoting_within_grp[1] / tot.cov.rate,
                                         se = msm::deltamethod(g = as.formula(paste("~", grpeq.mean, "* (", add , "+", paste(paste0(grp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ") /",
                                                                                    "(", tot.cov.rate, ") ")),
                                                               mean = c(m$coefficients),  cov = vcov, ses = TRUE) )
        m$sum.ls$sharedgrp_noncovoters <- c(est = grpeq.mean *
                                                (1 - m$sum.ls$covoting_within_grp[1]) / (1-tot.cov.rate),
                                            se = msm::deltamethod(g = as.formula(paste("~", grpeq.mean, "* (1-(", add , "+", paste(paste0(grp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ")) /",
                                                                                       "(1-(", tot.cov.rate, ")) ")),
                                                                  mean = c(m$coefficients),  cov = vcov, ses = TRUE))

        m$sum.ls$sharedgrp_ratio <- c(est = m$sum.ls$sharedgrp_covoters[1] / m$sum.ls$sharedgrp_noncovoters[1],
                                      se = msm::deltamethod(g = as.formula(paste("~",
                                                                                 "((", add , "+", paste(paste0(grp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ") / (",
                                                                                 tot.cov.rate, ")) /",
                                                                                 "((1 - (",add , "+",  paste(paste0(nogrp.mean, "*x", seq_along(m$coefficients)), collapse = " + "), ")) / (1 - (",
                                                                                 tot.cov.rate, "))) ")),
                                                            mean = c(m$coefficients),  cov = vcov, ses = TRUE))

        # Fix names
        for(i in seq_along(m$sum.ls)){
            if(length(m$sum.ls[[i]]) == 2)
                names(m$sum.ls[[i]]) <- c("est", "se")
        }

    } else {
        m$sum.ls <- NULL
    }


    # Make lean
    if(lean){
        m$call_env <- NULL
        m$fixef_id <- NULL
        m$obs_selection <- NULL
        m$weights <- NULL
        m$residuals <- NULL
        m$fitted.values <- NULL
        m$sumFE <- NULL
        m$scores <- NULL
        m$call <- NULL
    }

    # Return model
    return(m)
}


#' Make interactions between cleavages for fully interacted models
#'
#' @param df Pairwise dataset
#' @param main.var Main variable of interest
#' @param int.vars Control variables for which to build interactions.
#'    Can include either a CVR suffix or not. Must be in \code{df} for a and b.
#'
#' @import data.table
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#'
#' @returns An augmented dataset
cvr_cleav_int <- function(df,
                           main.var,
                           int.vars){
    # Check if data.table
    if(!"data.table" %in% class(df)){
        df <- data.table(df)
    }

    # Force factors
    for(v in c(main.var, int.vars[int.vars != main.var])){
        if(grepl("\\.eq|\\.uneq", v)){
            v_strip <- clean_cvr_suffix(v)
            if(!is.factor(df[,paste0("a.", v_strip)])){
                var <- paste0("a.", v_strip)
                df[, (var) := lapply(.SD, as.factor), .SDcols = var]
            }
            if(!is.factor(df[,paste0("b.", v_strip)])){
                var <- paste0("b.", v_strip)
                df[, (var) := lapply(.SD, as.factor), .SDcols = var]
            }
        }
    }

    # Clean variable names
    main.var = clean_cvr_suffix(main.var)
    int.vars = clean_cvr_suffix(int.vars)
    int.vars <- int.vars[int.vars != main.var]

    # Prepare and standardize controls
    std.ls <- lapply(c("a.","b."), function(d){

        # One-hot encode
        var.df <- model.matrix(as.formula(paste0("~ 0 + ", paste(paste0(d, int.vars), collapse = " + "))),data = df)

        # To data.table
        var.df <- data.table(var.df)

        # Center around 0
        cols_to_std <- colnames(var.df)
        var.df[, (cols_to_std) := lapply(.SD, scale, scale = FALSE), .SDcols = cols_to_std]

        # Rename
        colnames(var.df) <- substr(colnames(var.df), 3, nchar(colnames(var.df)))

        # Return
        return(var.df)
    })

    # Single
    single.df <- std.ls[[1]] + std.ls[[2]]
    colnames(single.df) <- paste0("int__cov", 1:ncol(single.df))

    # Covariate Interactions
    varseq <- 1:ncol(std.ls[[1]])
    cols <- colnames(std.ls[[1]])
    covint.df <- do.call(cbind, lapply(varseq, function(i){
        col.i <- cols[i]
        col.j <- cols[varseq >= i]
        d <- std.ls[[2]][,lapply(.SD, function(x){x * std.ls[[1]][[col.i]]}), .SDcols = col.j] +
            std.ls[[1]][,lapply(.SD, function(x){x * std.ls[[2]][[col.i]]}), .SDcols = col.j]
        colnames(d) <- paste0("int__c", i, "_c", varseq[varseq >= i])
        d
    }))

    # Main treatment: Prepare
    main.treat.ids <- df[[(paste0(main.var, ".eq"))]] == 1

    main.ls <- lapply(c("a.","b."), function(d){
        # One-hot encode
        var.df <- model.matrix(as.formula(paste0("~ 0 + ", paste(paste0(d, main.var), collapse = " + "))), data = df)

        # To data.table
        var.df <- data.table(var.df)

        # Rename
        colnames(var.df) <- substr(colnames(var.df), 3, nchar(colnames(var.df)))

        # Return
        return(var.df)

    })

    # Standardize
    main.ls.std <- lapply(seq_along(main.ls), function(i){
        m <- copy(main.ls[[i]])
        cols_to_std <- colnames(m)
        m[, (cols_to_std) := lapply(.SD, function(x){as.vector(scale(x, scale = FALSE))}), .SDcols = cols_to_std, by = main.treat.ids]
        m
    })

    # Main treatment: Single // Is included in Self-interactions
    main.single.df <- main.ls.std[[1]] + main.ls.std[[2]]
    colnames(main.single.df) <- paste0("int__", "main", 1:ncol(main.single.df))

    # Main treatment: Self-Interactions

    ## Prepare
    a.main <- paste0("a.", main.var)
    b.main <- paste0("b.", main.var)
    cell.mat <- cbind(as.numeric(unlist(df[, ..a.main])),
                      as.numeric(unlist(df[, ..b.main])))
    cell.mat[cell.mat[,1] > cell.mat[,2], c(1,2)] <- cell.mat[cell.mat[,1] > cell.mat[,2], c(2,1)]
    main.cells <- paste0("int__", main.var, cell.mat[,1], ".", cell.mat[,2])
    uni.vals <- sort(unique(as.vector(cell.mat)))
    main.vals <- sort(unique(main.cells))
    same.grp <- paste0("int__", main.var, uni.vals, ".", uni.vals)

    ## Execute
    for(j in seq_along(main.vals)[-1]){
        main.vals.j <- main.vals[j]
        covint.df[,(main.vals.j) := as.numeric(main.cells == main.vals[j])]
        if(main.vals[j] %in% same.grp){
            mean <- mean(covint.df[[main.vals.j]][main.treat.ids])
            covint.df[,(main.vals.j) := main.treat.ids * (covint.df[[main.vals.j]] - mean)]
        } else {
            mean <- mean(covint.df[[main.vals.j]][!main.treat.ids])
            covint.df[,(main.vals.j) := (!main.treat.ids) * (covint.df[[main.vals.j]] - mean)]
        }
    }

    # Interactions Main with Covariates
    m.varseq <- 1:ncol(main.ls[[1]])
    c.varseq <- 1:ncol(std.ls[[1]])
    covint.df <- cbind(covint.df, do.call(cbind, lapply(m.varseq, function(i){
        col.i <- colnames(main.ls[[1]])[i]
        d <- std.ls[[2]][,lapply(.SD, function(x){x * main.ls[[1]][[col.i]]})] +
            std.ls[[1]][,lapply(.SD, function(x){x * main.ls[[2]][[col.i]]})]
        colnames(d) <- paste0("int__m", i, "_c", c.varseq)
        d
    })))

    # Add to data
    df <- cbind(df, single.df, main.single.df, covint.df)

    ## Return
    return(df)
}
