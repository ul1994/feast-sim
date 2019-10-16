
source('metrics.R')
source('sim.R')

limitN <- 10000
subsample <- 100 # between every Nth

data <- readRDS('liat.rds')

nn <- floor((nrow(data)-1)*nrow(data)/2 / subsample)
# pairs <- matrix(,nrow=nn, ncol=nn)
mat <- matrix(,nrow=nn, ncol=3)


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

# system.time({
# 	# 13.927   1.179  28.696
# 	smt <- foreach(i=1:10000) %do% {
# 	pairs[inds[i,1], inds[i,2]] <- jsd(
# 		data[inds[i,1],][1:limitN],
# 		data[inds[i,2],][1:limitN])
# 	}
# })

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
