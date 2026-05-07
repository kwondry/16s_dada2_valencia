# Custom error learning script for platforms with binned quality scores
# (NextSeq, NovaSeq, iSeq). Replaces the default loessErrfun with a version
# that enforces monotonicity — error rates must not increase with quality.
#
# Based on community workaround from:
# https://github.com/benjjneb/dada2/issues/1307

library(dada2)
library(ggplot2)

# Modified loess error function: uses log10(tot) weights, span=2, and
# clamps any fitted value that dips below Q40 back up to the Q40 value.
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))

  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)

        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)

        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }

  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  # Enforce monotonicity: no fitted rate should be lower than the Q40 value
  estorig <- est
  est <- as.data.frame(est)
  q40_col <- ncol(est)  # last column = highest quality
  q40_vals <- est[[q40_col]]
  est <- as.matrix(est)
  for (i in seq_len(nrow(est))) {
    est[i, est[i, ] < q40_vals[i]] <- q40_vals[i]
  }
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)

  # Reconstruct the full 16x16 error matrix including self-transitions
  err <- rbind(
    1 - colSums(est[1:3, ]), est[1:3, ],
    est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
    est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
    est[10:12, ], 1 - colSums(est[10:12, ])
  )
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}

# --- Run learnErrors with the modified function ---

fls <- unlist(snakemake@input)
threads <- snakemake@threads

# Pass any extra params from the Snakemake rule (e.g., randomize)
extra <- snakemake@params
extra <- extra[names(extra) != "" & !is.na(names(extra))]

args <- c(
  list(
    fls = fls,
    errorEstimationFunction = loessErrfun_mod,
    multithread = threads
  ),
  extra
)

err <- do.call(learnErrors, args)

# Save error model
saveRDS(err, snakemake@output[["err"]])

# Save diagnostic plot
png(snakemake@output[["plot"]], width = 1200, height = 900, res = 150)
print(plotErrors(err, nominalQ = TRUE))
dev.off()
