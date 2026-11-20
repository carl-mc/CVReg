
#' Synthetic voting data
#'
#' Synthetically generated data for code examples.
#'
#' @examples
#' data(synthvote)
#' summary(synthvote)
#'
#' @name synthvote
#' @docType data
#' @keywords data


#' Transform dataset into pairwise observations.
#'
#' Function to transform a dataset into a set of pairs of the observations in it.
#' Its main use is to transform survey (or other) data into pairwise observations
#' to fit a covoting regression model. For further processing of the pairwise
#' data, see \code{process_pairs_default()}.
#'
#' @param df Dataset to transform. Can be \code{data.frame} or \code{data.table}.
#' @param vars Variables to keep for each observation after pairwise transformation. Defaults to NULL (all variables).
#' @param na.rm Drop observations with missing values on vars.
#' @param by.var Pairwise transformation only within clusters defined by by.vars.
#' @param structure Structure of pairs.
#' Can also by one of:
#' \itemize{
#'      \item{\code{"full"}:} Default. Each observation with every other observation, producing all undirected dyads.
#'      \item{\code{"dyads"}:} Each observation gets only paired with one other observation.
#'      \item{\code{"circle"}:} Observations get arranged in \code{degree/2} intersecting circles.
#'      \item{\code{"delaunay"}:} Pairs based on a
#'      [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of \code{struct.vars}.
#'      Results in a planar graph in those dimensions.
#'      \item{\code{"random"}:} Pairs observations completely at random such that each observation
#'      gets paired on average with \code{degree} other observations.
#'      \item{\code{"randomD"}:} Pairs observations based on a
#'      [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of
#'      \code{ndim} random variables drawn from the same uniform distribution \code{runif(0,1)}
#'      \item{\code{"all"}:} All directed pairs (without self-loops) are created.
#'      }
#' @param struct.vars Variables used for Delaunay triangulation when \code{structure = "delaunay"}
#' @param degree Number of edges per vertex, used when  \code{structure = "circle"} or \code{"random"}
#' @param ndim Number of dimensions for Delaunay triangulation of random variables, used when  \code{structure = "randomD"}
#' @param return_df Logical flag. If \code{TRUE}, returns a \code{data.frame} rather than \code{data.table}.
#'
#' @returns A \code{data.table} or a \code{data.frame}, if \code{return_df = TRUE}
#' in which each row is a pair of observations from the original dataset \code{df}. All variables \code{vars}
#' are kept using prefixes \code{c("a.", "b.")} for the first and second original observation, respectively,
#' Variables in \code{by.vars} are kept without prefix. Two new columns are created: \code{a.ID} and \code{b.ID}
#' record the original row of observation a and b in each pair as integer.
#'
#' @export
#'
#' @import data.table
#' @importFrom stats na.omit
#'
#' @examples
#' # Load some data
#' data(synthvote)
#'
#' # Transform to pairs
#' pair.dt <- data_to_pairs(df = synthvote, by.var = c("unit"), structure = "full")
#'
#' # Inspect
#' colnames(pair.dt) ## notice a. and b. prefix
#' nrow(pair.dt) ## 500 -> 62250 pairwise obs
data_to_pairs <- function(df,
                          vars = NULL,
                          na.rm = FALSE,
                          by.var = NULL,
                          structure = "full",
                          struct.vars = c(),
                          degree = NULL,
                          ndim = NULL, return_df = FALSE){

    # Ensure data.table
    if(!"data.table" %in% class(df)){
        df <- data.table(df)
    } else {
        df <- copy(df)
    }

    # Variables

    ## Select all if none selected
    if(is.null(vars)){
        vars <- colnames(df)
    }

    ## Error
    if(!all(vars %in% colnames(df))){
        stop(paste("Some variables are not in the data: ",
                   paste(vars[!vars %in% colnames(df)], collapse = ",")))
    }
    vars <- vars[vars %in% colnames(df)]

    ## Drop variables
    keep.vars <- unique(c(by.var, vars))
    df <- df[, ..keep.vars]

    # Make ID
    df$ID <- 1:nrow(df)

    # Drop missing data

    ## From all data
    if(na.rm){
        org.num <- nrow(df)
        na.omit(df, cols = unique(c(by.var, struct.vars, vars)))
        new.num <- nrow(df)
        if(new.num < org.num){
            message(paste0("Dropping ", org.num - new.num,
                           " observations with missing values (",
                           round(100*(org.num - new.num)/org.num,),
                           "%) from the data"))
        }

    ## Or only by vars
    } else if(!is.null(c(by.var, struct.vars))){
        org.num <- nrow(df)
        na.omit(df, cols = c(by.var, struct.vars))
        new.num <- nrow(df)
        if(new.num < org.num){
            message(paste0("Dropping ",
                           org.num - new.num,
                           " observations with missing values on by and structure variables (",
                           round(100*(org.num - new.num)/org.num,),
                           "%) from the data"))
        }
    }

    # Structure
    if(structure %in% c("full", "all")){
        # Simple full, undirected graph via data.table


        ## Duplicate data
        b.df <- copy(df)
        setnames(b.df, c(vars, "ID"), paste0("b.", c(vars, "ID")))
        setnames(df, c(vars, "ID"), paste0("a.", c(vars, "ID")))

        ## By?
        if(is.null(by.var)){
            by.var <- c("TEMP")
            df$TEMP <- b.df$TEMP <- 1L
        } else {
            for(v in by.var[by.var %in% vars]){
                setnames(b.df, paste0("b.",v), v)
                setnames(df, paste0("a.",v), v)
            }
        }

        ## Join
        df <- df[b.df,
                 on = by.var,
                 allow.cartesian = TRUE]
      
        ## Subset if undirected
        if(structure != "all"){
          df <- df[a.ID < b.ID]
        }
      
        ## Remove temp
        if("TEMP" %in% colnames(df)){
            df[, TEMP := NULL]
        }


    } else {
        # More Complex architectures via igraph

        ## Split data
        keep.vars <- c(by.var, struct.vars, "ID")
        if(is.null(by.var)){
            df.ls <- list(df[,..keep.vars])
        } else {
            df.ls <- split(df[,..keep.vars], df[, ..by.var])
            by.df <- do.call(rbind, lapply(df.ls, function(d){d[1, ..by.var]}))
        }


        ## Make graphs
        data.ls <- lapply(seq_along(df.ls), function(i){
            make_pairs_raw(df = df.ls[[i]],
                             structure = structure,
                             struct.vars = struct.vars,
                             degree = degree,
                             ndim = ndim)
        })

        # Delete NULLs
        del.vec <- sapply(data.ls, is.null)
        data.ls[del.vec] <- NULL
        by.df <- by.df[!del.vec,]

        # Combine edge data
        edf <- do.call(rbind, data.ls)

        # Add by.vars
        if(!is.null(by.var)){
            edf <- cbind(by.df[rep(1:nrow(by.df), sapply(data.ls, nrow)), ], edf)
        }

        # Merge with main data

        ## Duplicate data
        b.df <- copy(df)
        colnames(b.df)[colnames(b.df) %in% c(vars, "ID")] <- paste0("b.", colnames(b.df)[colnames(b.df) %in% c(vars, "ID")])
        colnames(df)[colnames(df) %in% c(vars, "ID")] <- paste0("a.", colnames(df)[colnames(df) %in% c(vars, "ID")])

        ## By?
        if(is.null(by.var)){

        } else {
            setnames(b.df, paste0("b.",by.var[by.var %in% vars]), by.var[by.var %in% vars])
            setnames(df, paste0("a.",by.var[by.var %in% vars]), by.var[by.var %in% vars])
        }

        ## Join
        edf <- edf[df,
                   on = c(by.var, "a.ID"),
                   nomatch = NULL]
        df <- edf[b.df,
                   on = c(by.var, "b.ID"),
                  nomatch = NULL]

        # Clean up
        rm("edf")

    }

    # Clean up
    rm("b.df")

    # Transform to data.frame?
    if(return_df){
        df <- data.frame(df)
    }

    # Return
    return(df)
}


#' Make a raw list of pairs of observations in dataset
#'
#' Takes observations in dataset \code{df} and creates a list of pairs or graph of a given structure.
#' If \code{return_graph = FALSE}, the result is a \code{data.table} in which, for each pair, a.ID and b.ID
#' indicate the row in \code{df} of the first and second entry, respectively.
#'
#'
#' @param df Data to based pair creation on. Either a \code{data.table} or \code{data.frame} object
#' @param structure Structure of pairs.
#' Can also by one of:
#' \itemize{
#'      \item{\code{"full"}:} Default. Each observation with every other observation, producing all undirected dyads.
#'      \item{\code{"dyads"}:} Each observation gets only paired with one other observation.
#'      \item{\code{"circle"}:} Observations get arranged in \code{degree/2} intersecting circles. The resulting object contains a
#'         variable "dim" which denotes the circle an edge is part of.
#'      \item{\code{"delaunay"}:} Pairs based on a
#'      [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of \code{struct.vars}.
#'      Results in a planar graph in those dimensions.
#'      \item{\code{"random"}:} Pairs observations completely at random such that each observation
#'      gets paired on average with \code{degree} other observations.
#'      \item{\code{"randomD"}:} Pairs observations based on a
#'      [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of
#'      \code{ndim} random variables drawn from the same uniform distribution \code{runif(0,1)}
#'      \item{\code{"all"}:} All directed pairs (without self-loops) are created.
#'      }
#' @param struct.vars Variables used for Delaunay triangulation when \code{structure = "delaunay"}
#' @param degree Number of edges per vertex, used when  \code{structure = "circle"} or \code{"random"}
#' @param ndim Number of dimensions for Delaunay triangulation of random variables, used when  \code{structure = "randomD"}
#' @param return_graph Logical flag of whether to return an \code{igraph} object (TRUE) or a \code{data.table} (FALSE). Defaults to FALSE.
#'
#' @returns An \code{igraph} object (\code{return_graph = TRUE}) or a \code{data.table} (\code{return_graph = FALSE})
#' @export
#' @import igraph
#' @import data.table
#' @importFrom stats na.omit
#'
#' @examples
#'
#' # Load some data
#' data(synthvote)
#'
#' # Transform to pairs
#' pairs.raw <- make_pairs_raw(df = synthvote, structure = "circle", return_graph = FALSE)
#'
#' # Inspect
#' head(pairs.raw)
make_pairs_raw <- function(df,
                             structure = "random",
                             struct.vars = c(),
                             degree = 6,
                             ndim = NULL,
                             return_graph = FALSE){

    # Make data.table
    if(!"data.table" %in% class(df)){
        df <- data.table(df)
    }

    # Check if more than one row
    if(nrow(df) < 2){
        return(NULL)
    }

    # Add ID column if does not exist
    if(!"ID" %in% colnames(df)){
        df[, ID := seq_len(nrow(df))]
    }

    # Check for duplicates
    stopifnot(anyDuplicated(df[, ID]) == 0)

    # Randomize order
    df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE),]

    # Variables for structure
    if(structure == "ndim"){
        na.omit(df, cols = struct.vars)
    }

    ## Return NULL if no data is left
    if(nrow(df) < 2){
        return(NULL)
    }

    # Graph structures
    if (structure == "dyads"){
        ## Even number of rows
        if(nrow(df) %% 2 == 1){
            df <- df[-(sample(1:nrow(df), 1)),]

        }

        ## Make graph
        graph <- dyad_graph(N = nrow(df), return_edgemat = TRUE)


    } else if(structure == "circle"){
        if(is.null(degree)){
            stop("Please specify degree")
        }

        ## Adjust degree
        if(degree > nrow(df)){
            degree <- 2 * (ceiling((nrow(df) - 1)/2))
        }

        ## Make random structure graph
        graph <- circle_graph(N = nrow(df),
                               degree = degree, return_edgemat = TRUE)


    } else if(structure == "random"){
        if(is.null(degree)){
            stop("Please specify degree")
        }

        ## Make random structure graph
        graph <- random_graph(N = nrow(df),
                              degree = degree, return_edgemat = TRUE)

    } else if(structure == "randomD"){
        if(is.null(ndim)){
            stop("Please specify dimensionality")
        }

        ## Make N-dimensional random structure graph
        graph <- random_ndim_graph(N = nrow(df),
                                   Ndim = ndim)

    } else if(structure == "all"){
        graph <- expand.grid(1:nrow(df), 1:nrow(df))
        graph <- graph[graph[,1] != graph[,2],]

    } else if(structure == "full"){
        graph <- expand.grid(1:nrow(df), 1:nrow(df))
        graph <- graph[graph[,1] < graph[,2],]

    } else if(structure == "delaunay"){
        if(is.null(struct.vars)){
            stop("Please specify structure variables")
        }

        ## Make graph
        graph <- make_ndim_graph(data = data.frame(df[, ..struct.vars]))

    }

    # Finalize edge level dataset
    if("igraph" %in% class(graph)){
        edf <- igraph::as_data_frame(graph, what = "edges")
        edf <- data.table(edf)
    } else {
        edf <- data.table(graph)
        colnames(edf)[1:2] <- c("from", "to")
    }


    ## Set to and from to original IDs
    edf[, a.ID := df$ID[as.numeric(from)]]
    edf[, b.ID := df$ID[as.numeric(to)]]
    edf[, c("from", "to") := NULL]

    # Return
    if(return_graph){
        # Make a graph object
        e.vrs <- c(from = "a.ID", to = "b.ID")
        setnames(edf, e.vrs, names(e.vrs))
        graph <- igraph::graph_from_data_frame(data.frame(edf))

        # Return
        return(graph)
    } else {
        return(edf)
    }

}




clean_cvr_suffix <- function(v){
    suf <- get_suffix(v)
    sapply(seq_along(v), function(i){
        gsub(paste0(".", suf[i]), "", v[i], fixed = T)
    })
}


#' Get CVR suffixed from variable names
#'
#' @param v Variable names
#'
#' @returns Just the CVR suffix
#'
#' @importFrom utils tail
get_suffix <- function(v){
    s <- strsplit(v, ".", fixed = T)
    sapply(s, function(i){
        if(length(i) == 1){
            NA
        } else {
            tail(i, 1)
        }
    })
}


cvr_suffix <- function(){
    suffix.ls <- list(c(suffix = ".eq", description = "Tests if a and b are equal"),
                      c(suffix = ".uneq", description = "Tests if a and b are unequal"),
                      c(suffix = ".diff", description = "Computes absolute difference between a and b"),
                      c(suffix = ".sim", description = "Computes -1 * absolute difference between a and b"),
                      c(suffix = ".max", description = "Computes max of a and b"),
                      c(suffix = ".min", description = "Computes min of a and b"),
                      c(suffix = ".mean", description = "Coputes mean of a and b"))
    data.frame(do.call(rbind, suffix.ls))
}


#' Process pair data based on common comparison functions
#'
#' @param df Pair dataset, with variables of observations a and b denoted by prefix "a." and "b." respectively.
#'      Best use a dataset as returned by \link[CVReg]{data_to_pairs}. Can be either a \code{data.table} or \code{data.frame}.
#' @param vars Variables to create based on comparison functions indicated by the following suffixes:
#' \itemize{
#'      \item{\code{".eq"}:} Binary equality indicator. E.g, \code{group.eq = as.numeric(a.group == b.group)}
#'      \item{\code{".uneq"}:}  Binary inequality indicator. E.g,  \code{group.uneq = as.numeric(a.group != b.group)}
#'      \item{\code{".diff"}:} Absolute difference. E.g,  \code{income.diff = abs(a.income - b.income)}
#'      \item{\code{".sim"}:}  Similarity as -1 x (absolute difference). E.g,  \code{income.sim = -1 * abs(a.income - b.income)}
#'      \item{\code{".max"}:} Pairwise maximum. E.g,  \code{income.max = max(a.income, b.income)}
#'      \item{\code{".min"}:} Pairwise minimum. E.g,  \code{income.min = min(a.income, b.income)}
#'      \item{\code{".mean"}:} Pairwise mean. E.g,  \code{income.mean = mean(a.income, b.income)}
#' }
#' @param fun.ls Named list of functions to process \code{df} with. The input and sole argument of each
#'      function must be the pair-wise dataset \code{df} and the function must return the augmented \code{df}
#'      (e.g., the \code{df} with new variables added).
#'      Note that per default, \code{df} is passed as a \code{data.table}. If
#'      the functions are based on the simpler (but slower) \code{data.frame} notation
#'      and include operations that are inadmissible in a \code{data.table} environment,
#'      set flag \code{fun_data_frame = TRUE}.
#'      Can be used to generate custom-made variables. E.g., do create the Cartesian distance between a and b
#'      on the basis of variables x and y use:
#'      \code{fun.ls = list(function(df){df[, dist :=((a.x - b.x)^2 + (a.y - b.y)^2)^.5]})}
#' @param fun_data_frame Logical flag on whether fun.ls is based on
#'      \code{data.frame} (TRUE) or \code{data.table} (FALSE, the default) notation.
#'
#' @returns A \code{data.table} object
#' @export
#' @importFrom  matrixStats rowMins
#' @importFrom  matrixStats rowMaxs
#' @importFrom  matrixStats rowMeans2
#' @import data.table
#' @examples
#' # Load some data
#' data(synthvote)
#'
#' # Transform to pairs
#' pair.dt <- data_to_pairs(df = synthvote, by.var = c("unit"), structure = "full")
#'
#' # Process pair level data to create pairwise similarity/difference measures
#' pair.vars <- c("prtvt.eq", "group.eq", "income.sim")
#' pair.dt <- process_pairs(df = pair.dt,
#'                          vars = pair.vars,
#'                          fun.ls = list(function(df){
#'                             df[, dist :=((a.x - b.x)^2 + (a.y - b.y)^2)^.5]
#'                          }))
#'
#' # Inspect
#' subset(pair.dt,,c(pair.vars, "dist"))
process_pairs <- function(df, vars = c(), fun.ls = list(), fun_data_frame = FALSE){
    # Process
    for(v in vars){
        # Process
        suf <- get_suffix(v)
        var <- gsub(paste0(".", suf), "", v, fixed = T)
        a.var <- paste0("a", ".", var)
        b.var <- paste0("b", ".", var)

        if(suf == "eq"){
            df[, (v) := as.numeric(.SD[[1]] == .SD[[2]]), .SDcols = c(a.var, b.var)]

        } else if(suf == "uneq"){
            df[, (v) := as.numeric(.SD[[1]] != .SD[[2]]), .SDcols = c(a.var, b.var)]

        } else if(suf == "diff"){
            if(!is.numeric(df[[a.var]])){warning(paste("Transforming", a.var, "to numeric."))}
            if(!is.numeric(df[[b.var]])){warning(paste("Transforming", b.var, "to numeric."))}
            df[, (v) := abs(as.numeric(.SD[[1]]) - as.numeric(.SD[[2]])), .SDcols = c(a.var, b.var)]

        } else if(suf == "sim"){
            if(!is.numeric(df[[a.var]])){warning(paste("Transforming", a.var, "to numeric."))}
            if(!is.numeric(df[[b.var]])){warning(paste("Transforming", b.var, "to numeric."))}
            df[, (v) := -1 * abs(as.numeric(.SD[[1]]) - as.numeric(.SD[[2]])), .SDcols = c(a.var, b.var)]

        } else if(suf == "max"){
            if(!is.numeric(df[[a.var]])){warning(paste("Transforming", a.var, "to numeric."))}
            if(!is.numeric(df[[b.var]])){warning(paste("Transforming", b.var, "to numeric."))}
            df[, (v) := rowMaxs(as.matrix(as.numeric(.SD))), .SDcols = c(a.var, b.var)]

        } else if(suf == "min"){
            if(!is.numeric(df[[a.var]])){warning(paste("Transforming", a.var, "to numeric."))}
            if(!is.numeric(df[[b.var]])){warning(paste("Transforming", b.var, "to numeric."))}
            df[, (v) := rowMins(as.matrix(as.numeric(.SD))), .SDcols = c(a.var, b.var)]

        }  else if(suf == "mean"){
            if(!is.numeric(df[[a.var]])){warning(paste("Transforming", a.var, "to numeric."))}
            if(!is.numeric(df[[b.var]])){warning(paste("Transforming", b.var, "to numeric."))}
            df[, (v) := rowMeans2(as.matrix(as.numeric(.SD))), .SDcols = c(a.var, b.var)]

        } else {
            stop(paste("suffix", suf, "not in default repertoire. "))
        }
    }

    # Generate additional variables based on fun.ls

    ## Make data.frame if necessary
    if(fun_data_frame){
        df <- data.frame(df)
    }

    ## Process
    for(f in seq_along(fun.ls)){
        df <- fun.ls[[f]](df)
    }

    ## Retransform to data.table
    if(fun_data_frame){
        df <- data.table(df)
    }



    # Return
    return(df)
}


