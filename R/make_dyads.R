


#' Atomistic Dyad Graph
#'
#' @param N Number of vertices
#' @param return_edgemat Return edge matrix rather than an igraph object.
#'
#' @returns Either an edge matrix or an igraph object
#' @export
#' @import igraph
#'
#' @examples
#' g <- dyad_graph(10, return_edgemat = FALSE)
#' plot(g)
dyad_graph <- function(N, return_edgemat = FALSE){
    ## Check
    if(N %% 2 == 1){
        stop("N needs to be an even number.")
    }

    ## Random order
    order <- sample(1:N, size = N, replace = F)

    ## Edges
    edge.mat <- cbind(order[1:(N/2)], order[(N/2+1):N])

    ## Return
    if(return_edgemat){
        colnames(edge.mat) <- c("from", "to")
        return(edge.mat)
    } else {
        g <- graph_from_edgelist(edge.mat[,1:2], directed = FALSE)
        return(g)
    }
}


#' Circle with supra-local structure
#'
#' @param N Number of vertices
#' @param degree Edge degree (number of edges per vertex)
#' @param return_edgemat Return edge matrix rather than an igraph object.
#'
#' @returns Either an edge matrix or an igraph object
#' @export
#' @import igraph
#'
#' @examples
#' g <- circle_graph(10, degree = 2, return_edgemat = FALSE)
#' plot(g)
circle_graph <- function(N, degree = 2, return_edgemat = FALSE){
    ## Check
    stopifnot(degree <= N)
    stopifnot(degree %% 2 == 0)


    ## Fully connected variant
    if(degree >= N-1){
        edge.mat <- expand.grid(seq_len(N), seq_len(N))
        edge.mat <- edge.mat[edge.mat[,1] < edge.mat[,2],]
        g <- graph_from_edgelist(as.matrix(edge.mat[,1:2]), directed = FALSE)
        E(g)$dim <- -99
        return(g)
    }

    ## Make circles
    pot.targets <- rep(list(1:N), N)
    pot.targets <- lapply(seq_along(pot.targets), function(i){
        base.t <- i + c(-1,0,1)
        if(i == 1){
            base.t <- c(N, base.t)
        } else if(i == N){
            base.t <- c(1, base.t)
        }
        pot.targets[[i]][!pot.targets[[i]] %in% base.t]
    })
    edge.mat <- cbind(1:N, c(N, 1:(N-1)), rep(1, N))
    c <- 1
    rep = FALSE
    target <- ifelse(degree > (N-1)/2, floor((N-1)/2) - degree/2, degree/2)
    while(c < target){
        c = c + 1
        # print(c)

        i = 0
        n = 1
        v = c()
        temp.edges <- NULL
        temp.targets <- pot.targets
        rep = FALSE
        while(i < N){
            i = i + 1
            v <- c(v, n)
            if(i == N){ ## Back to start
                v <- v[-1]
            }

            pool <- temp.targets[[n]][!temp.targets[[n]] %in% v]
            if(length(pool) > 1){
                m = sample(pool, size = 1)
            } else if(length(pool) == 0){
                rep = T
                c = c-1
                break
            } else {
                m <- pool
            }

            ## Add Edge
            temp.edges <- rbind(temp.edges,
                                c(n, m))

            ## Delete from potential set
            temp.targets[[n]] <- temp.targets[[n]][temp.targets[[n]] != m]
            temp.targets[[m]] <- temp.targets[[m]][temp.targets[[m]] != n]

            ## Move to next round
            n <- m
            #
            # print(i)
            # print(n)
            # print(pot.targets[[7]])
        }

        if(!rep){
            edge.mat <- rbind(edge.mat,
                              cbind(temp.edges, rep(c, N)))
            pot.targets <- temp.targets
        }
    }

    # Subtract edges if we are in deletion mode
    if(target < degree/2){
        all.edge.mat <- as.matrix(expand.grid(seq_len(N), seq_len(N)))
        all.edge.mat <- all.edge.mat[all.edge.mat[,1] < all.edge.mat[,2],]
        edge.mat[edge.mat[,1] > edge.mat[,2],c(1,2)] <- edge.mat[edge.mat[,1] > edge.mat[,2], c(2,1)]
        sub.edge.mat <- rbind(all.edge.mat, edge.mat[,c(1,2)])
        sub.edge.mat <- sub.edge.mat[!duplicated(sub.edge.mat) &
                                         !duplicated(sub.edge.mat, fromLast = T),]

        # Drop random dyads if necessary
        if(N %% 2 == 0){
            for(i in 1:N){
                if(sum(sub.edge.mat == i) <= degree){next}
                t <- c(sub.edge.mat[sub.edge.mat[,1] == i, 2],
                       sub.edge.mat[sub.edge.mat[,2] == i, 1])
                t <- t[sapply(t, function(j){sum(sub.edge.mat == j) > degree})]
                if(length(t) == 0){next}
                tt <- sample(t, size = 1)
                sub.edge.mat <- sub.edge.mat[!(sub.edge.mat[,1] == i &
                                                   sub.edge.mat[,2] == tt), ]
                sub.edge.mat <- sub.edge.mat[!(sub.edge.mat[,2] == i &
                                                   sub.edge.mat[,1] == tt), ]
            }
        }
        edge.mat <- sub.edge.mat
        edge.mat <- cbind(edge.mat, rep(NA, nrow(edge.mat)))
    }


    ## Return
    if(return_edgemat){
        colnames(edge.mat) <- c("from", "to", "dim")
        return(edge.mat)
    } else {
        g <- graph_from_edgelist(edge.mat[,1:2], directed = FALSE)
        E(g)$dim <- edge.mat[,3]
        return(g)
    }
}


#' Random graph
#'
#' @param N Number of vertices
#' @param degree Average edge degree (avg. number of edges per vertex)
#' @param return_edgemat Return edge matrix rather than an igraph object.
#'
#' @returns Either an edge matrix or an igraph object
#' @export
#' @import igraph
#'
#' @examples
#' g <- random_graph(10, degree = 2, return_edgemat = FALSE)
#' plot(g)
random_graph <- function(N, degree = 4, return_edgemat = FALSE){
    ## Degree
    eN <- N * degree / 2
    eN <- min(c(eN, (N*N - N)/2))

    ## Edge list
    edge.mat <- expand.grid(1:N, 1:N)
    edge.mat <- edge.mat[edge.mat[,1] < edge.mat[,2],]
    keep.e <- sample(seq_len(nrow(edge.mat)),
                     eN,
                     replace = F)
    edge.mat <- edge.mat[keep.e,]
    edge.mat <- as.matrix(edge.mat)

    ## Return
    if(return_edgemat){
        colnames(edge.mat) <- c("from", "to")
        return(edge.mat)
    } else {
        g <- make_graph(as.vector(t(edge.mat)), n = N, directed = FALSE)
        return(g)
    }
}



#' Planar graph via Delauney tesselation of Ndim random variables
#'
#' @param N Number of vertices
#' @param Ndim Number of dimensions (random variables)
#' @param return_edgemat Return edge matrix rather than an igraph object.
#'
#' @returns Either an edge matrix or an igraph object
#' @export
#' @import igraph
#' @importFrom stats runif
#'
#' @examples
#' g <- random_ndim_graph(10, Ndim = 2)
#' plot(g)
random_ndim_graph <- function(N, Ndim, return_edgemat = FALSE){
    ## Draw random variables
    data <- do.call(cbind, lapply(seq_len(Ndim), function(d){
        runif(N)
    }))
    colnames(data) <- paste0("rdim_", 1:Ndim)

    ## Make graph
    g <- make_ndim_graph(data)

    ## Return
    if(return_edgemat){
        return(igraph::as_data_frame(g, what = "edges"))
    } else {
        return(g)
    }
}

#' Planar graph via Delauney tesselation of data
#'
#' @param data Data to use for Delauney tesselation.
#' @param  return_edgemat Return edge matrix rather than an igraph object.
#'
#' @returns Either an edge matrix or an igraph object
#' @export
#' @import igraph
#' @importFrom geometry delaunayn
#'
#' @examples
#' g <- make_ndim_graph(data.frame(x = runif(10), y = runif(10)))
#' plot(g)
make_ndim_graph <- function(data, return_edgemat = FALSE){

    ## Name rows
    rownames(data) <- as.character(1:nrow(data))

    ## Make numeric
    data.num <- apply(data, 2, function(x){
        if(is.character(x)){
            as.numeric(as.factor(x))
        } else if(is.factor(x)){
            as.numeric(x)
        } else {
            x
        }
    })
    colnames(data.num) <- colnames(data)

    ## Delete redundant dimensions
    data.num <- data.num[,apply(data.num, 2, function(x){length(unique(x)) > 1})]

    ## Delaunay Triangulation
    tc <- delaunayn(data.num)

    ## Edges
    emat <- do.call(rbind, apply(tc, 1, function(x){
        expand.grid(x, x)
    }))
    emat <- rbind(as.matrix(emat), cbind(emat[,2], emat[,1]))
    emat <- emat[emat[,1] < emat[,2],]
    emat <- unique(emat)
    emat <- emat[order(emat[,1]),]

    ## Graph
    graph <- graph_from_data_frame(d = emat,
                                   directed = FALSE)

    ## ID Variable
    V(graph)$ID <- as.numeric(V(graph)$name)
    graph <- igraph::permute(graph, V(graph)$ID)

    ## Vertex attributes
    for(c in seq_len(ncol(data))){
        vertex_attr(graph, colnames(data)[c]) <-
            data[,c]
    }


    ## Return
    if(return_edgemat){
        return(igraph::as_data_frame(graph, what = "edges"))
    } else {
        return(graph)
    }
}
