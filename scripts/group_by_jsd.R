
########################################################################
#' group_by_jsd.R
#'
#' Taking the output of `find_jsd_pairs.R` finds groups of sources which
#' have approximately a desired pairwise avg JSD
#'
#' @param batch For bookkeeping.
#'  Umbrella folder where everything will be saved inside
#' @param saved Saved output of `find_jsd_pairs.R`
#' @param data Raw EMP data
#' @param targets Desired JSDs to find groups for
#' @param thresh Amount of avg JSD change at which sampled pairs
#'  will be rejected from a group
#' @param sample_range Max range to sample before and after where
#'  the desired JSD appears after sorting
#' @param nSources Number of sources to be collected per group
#'
########################################################################

library('permute')
source('metrics.R')

batch<- '10k'
saved <- readRDS('jsd_10k.rds')
data <- readRDS('liat.rds')
targets <- c(0.08, 0.125, 0.50, 0.90)
thresh <- 0.005
sample_range <- 400
nSources <- 21

scores <- saved[,3] # Precomputed JSD is in 3rd col
scores[is.na(scores)] <- -1 # this has some NAs...
ord <- order(scores)

########################################################################
#' For each desired JSD, save where they appear in the sorted list
########################################################################

track_index <- matrix(0, nrow=length(targets), ncol=2)
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

# Helper function to construct sources matrix given their inds in EMP
collect_rows <- function(inds, limitN=10000) {
	mat <- matrix(, nrow=length(inds), ncol=limitN)
	for (ii in 1:length(inds)) {
		mat[ii,] <- data[inds[ii],1:limitN]
	}
	return(mat)
}

########################################################################
#' Try 3 times per desired JSD to find a group of sources
#' which approx. has that average JSD.
#'
#' Groups are saved to `saved/jsd/{batch}/{desired JSD}/${saved sources mat}.rds`
#'
#' NOTE: nSources should be K+1 to include an unknown source
#'
########################################################################


for (ti in 1:length(targets)) {
	target <- targets[ti]
	print(paste('Target', target))
	median <- track_index[ti, 1]

	for (si in 1:10) {
		# sample before and after where the JSD appearss
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
