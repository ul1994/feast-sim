
########################################################################
#' find_jsd_pairs.R
#'
#' This code calculates the JSD between many pairs of sources
#' in the dataset.
#'
#' @param limitN Number of taxa can be limited for small-scale tests
#' @param subsample Reduce the # of sources considered for faster results
#' @param data Sources x Taxa raw EMP data
#'
########################################################################


source('metrics.R')
source('sim.R')

limitN <- 10000
subsample <- 100 # Every Nth will be considered
data <- readRDS('liat.rds')

# Count of sources to be considered
nn <- floor((nrow(data)-1)*nrow(data)/2 / subsample)

# entries will be (ith source, jth source, jsd value)
mat <- matrix(,nrow=nn, ncol=3)

########################################################################
#' Collect all (i, j) indicies to be considered.
#' Iterating over (i, j)s can help later if jsd computation
#' is to be parallelized
########################################################################

inds <- matrix(,nrow=nn, ncol=2)
it <-1
ij_ind <- 1
for (ii in 1:(nrow(data)-1)) {
	for (jj in (ii+1):nrow(data)) {
		if (it %% subsample == 0) {
			if (ij_ind <= nn) {
    			inds[ij_ind, ] <- c(ii, jj)
				ij_ind <- ij_ind + 1
			}
		}
		it <- it + 1
	}
}


########################################################################
#' Iterate over each (i, j) and compute their jsd.
#' The result is saved into mat as (i, j, jsd value)
########################################################################

for(i in 1:nrow(inds)) {
	v1 <- data[inds[i,1],][1:limitN]
	v2 <- data[inds[i,2],][1:limitN]
	v1 <- v1 / sum(v1)
	v2 <- v2 / sum(v2)

	mat[i,] <- c(
		inds[i,1],
		inds[i,2],
		jsd(v1, v2)
	)
	if (i %% 100 == 0) print(paste(i, nrow(inds)))
}

########################################################################
#' Remember to save mat somewhere
########################################################################

saveRDS(mat, 'jsd_10k.rds')
