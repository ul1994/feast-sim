
library('permute')
source('metrics.R')

# saved <- readRDS('jsd.rda')
batch<- '10k'
saved <- readRDS('jsd_10k.rds')
data <- readRDS('liat.rds')
scores <- saved[,3]

targets <- c(0.08, 0.125, 0.50, 0.90)
thresh <- 0.005
nSources <- 21

# this has some NAs...
scores[is.na(scores)] <- -1
sample_range <- 400
# targets <- c(0.125, 0.25, 0.5, 0.75)
# targets <- c(0.25)
track_index <- matrix(0, nrow=length(targets), ncol=2)
ord <- order(scores)

# assuming scores are spread nicely, find ranges to sample from
starting <- 1
for (ii in 1:length(ord)) {
	ind <- ord[ii]
	if (scores[ind] == -1) {
		starting <- starting + 1
		next
	}

	ent <- saved[ind,]
	for (ti in 1:length(targets)) {
		target <- targets[ti]
		prev <- track_index[ti, 2]
		if (abs(scores[ind] - target) < abs(prev - target)) {
			track_index[ti, 1] <- ii
			prev <- track_index[ti, 2] <- scores[ind]
		}
	}
}

print(track_index)
print(length(ord))
# print(saved[track_index[1,1]-1,])

collect_rows <- function(inds, limitN=10000) {
	mat <- matrix(, nrow=length(inds), ncol=limitN)
	for (ii in 1:length(inds)) {
		mat[ii,] <- data[inds[ii],1:limitN]
	}
	return(mat)
}

for (si in 1:3) {
	# sample around the indicies found
	for (ti in 1:length(targets)) {
		target <- targets[ti]
		median <- track_index[ti, 1]
		check <- shuffle(-sample_range:sample_range)

		group <- c()
		running_avg <- -1
		for (offset in check) {
			index <- ord[median + offset]
			entry <- saved[index,]

			test_group <- group
			for (candidate in entry[1:2]) {
				candidate <- as.integer(candidate)
				if (!(candidate %in% test_group)) {
					test_group <- c(test_group, candidate)
				}
			}

			if (min(rowSums(collect_rows(test_group))) < 250) {
				next
			}
			test_avg <- jsdavg(collect_rows(test_group))
			# print(test_avg)
			if (running_avg == -1 || abs(test_avg - running_avg) < thresh) {
				# acceptable jsd
				# print(test_avg)
				group <- test_group
				running_avg <- test_avg
			}

			if (length(group) >= nSources) break
		}

		group <- group[1:nSources]
		final_avg <- jsdavg(collect_rows(group))
		print(paste('Obtained:', target, length(group), final_avg))
		saveRDS(
			collect_rows(group),
			file=sprintf(
				'saved/jsd/%s/%04d/sources_jsd_%04d_%06d.rds',
				batch,
				as.integer(target * 1000),
				as.integer(target * 1000),
				as.integer(final_avg * 100000)))
	}

}
