#' Information--resolution profile of an observed-entry hierarchy
#'
#' Replays the chronological cuts of a declared CAPIVARA hierarchy. The
#' representation, eligible footprint and measured entries remain fixed across
#' cuts. Empirical distortion is weighted observed-entry mean squared error;
#' it is not Shannon information or independent predictive evidence.
#'
#' Optional independent validation compares each discovery reconstruction with
#' the finest evaluated reconstruction. The difference of squared losses has
#' an expectation equal to the difference of latent squared losses under the
#' declared unbiased measurement model. Marginal variances give a conservative
#' standard-deviation bound even under unknown spatial/spectral correlation.
#' Simultaneous one-sided bounds use Gaussian or Cantelli inequalities with
#' Bonferroni correction. A selection is conditional on these assumptions,
#' the reference resolution, the feature weights and the scientific tolerance.
#' @param segmentation Result of segment_large(..., representation=...,
#'   return_details=TRUE). The original hierarchy and its prepared representation
#'   are required; legacy sparse partitions are not accepted.
#' @param candidate_k Integer candidate region counts. Counts below unresolved
#'   roots or above the eligible pixel count are reported as unavailable.
#'   The coarsest feasible count is always included.
#' @param weights Optional nonnegative vector of weights on selected native
#'   feature coordinates; at least one observed entry must receive positive weight.
#'   The weights must be fixed without seeing validation observations.
#' @param validation Optional list with aligned feature and marginal variance
#'   matrices (features, variance), representation_id, selected_channel_indices,
#'   source, independent=TRUE, unbiased=TRUE, and error_model ("gaussian" or
#'   "variance_only"). These are caller-supplied assumptions, not verified by
#'   the package. Validation must be finite on the entire fixed scoring domain.
#'   For nonlinear features, the caller must establish the target mean and
#'   marginal variances in feature space, not pass raw flux variances.
#' @param tolerance Nonnegative tolerated excess loss in weighted squared feature
#'   units relative to the finest evaluated cut. NULL disables selection.
#' @param alpha Family-wise one-sided error probability.
#' @return List containing a profile data frame, native label maps, available
#'   and unavailable counts, and a conditional selection or explicit abstention.
#' @export
partition_resolution <- function(segmentation, candidate_k=seq_len(32L),
                                 weights=NULL, validation=NULL,
                                 tolerance=NULL, alpha=0.05) {
  h <- segmentation$hierarchy
  p <- segmentation$prepared_representation
  if (!inherits(h,"capivara_observed_hierarchy") ||
      !inherits(p,"capivara_prepared_representation") ||
      !identical(segmentation$representation_id,p$representation_id))
    stop("Require declared hierarchy with return_details=TRUE.",call.=FALSE)
  final <- segmentation$analysis_support$final_analysis_mask
  if (!is.logical(final) || anyNA(final) || !identical(dim(final),dim(segmentation$cluster_map)))
    stop("Invalid frozen analysis support.")
  original <- which(final)
  n <- length(original)
  m <- nrow(h$children)
  if (n<1L || is.null(m) || m!=length(h$cost) || length(h$labels)!=n ||
      !is.matrix(h$children) || ncol(h$children)!=2L ||
      n-m!=length(h$roots) || any(!is.finite(h$cost)|h$cost<0))
    stop("Invalid chronological hierarchy.")
  if (!is.numeric(candidate_k) || !length(candidate_k) ||
      any(!is.finite(candidate_k)|candidate_k<1|candidate_k!=floor(candidate_k)|
            candidate_k>.Machine$integer.max))
    stop("candidate_k must contain positive integers.")
  if (!is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("alpha must lie in (0,1).")
  if (!is.null(tolerance) && (!is.numeric(tolerance)||length(tolerance)!=1L||
      !is.finite(tolerance)||tolerance<0)) stop("tolerance must be finite and nonnegative.")
  candidate_k <- sort(unique(as.integer(candidate_k)))
  root_k <- n-m
  available <- sort(unique(c(root_k,candidate_k[candidate_k>=root_k & candidate_k<=n])))
  unavailable <- setdiff(candidate_k,available)
  if (length(available)>256L) stop("Request at most 256 feasible cuts per profile.")
  X <- p$features[original,,drop=FALSE]
  q <- p$sample_validity[original,,drop=FALSE]
  if (!is.numeric(X)||!is.logical(q)||anyNA(q)||!identical(dim(X),dim(q))||
      any(!is.finite(X[q]))) stop("Invalid prepared measured features.")
  nc <- ncol(X)
  if (is.null(weights)) weights <- rep(1,nc)
  if (!is.numeric(weights)||length(weights)!=nc||any(!is.finite(weights)|weights<0))
    stop("weights must be nonnegative and aligned to selected feature coordinates.")
  # The target domain is fixed by discovery, not by candidate-specific coverage
  # or the realised validation fluxes. Zero-weight entries remain excluded.
  domain <- q & matrix(rep(weights>0,each=n),nrow=n)
  W <- sum(colSums(domain)*weights)
  if (!is.finite(W)||W<=0) stop("No positive weighted observed entries.")
  X[!q] <- 0
  # Replay full chronological history once; preserve raw inversions and roots.
  group <- seq_len(n)
  members <- vector("list", n+m)
  for (i in seq_len(n)) members[[i]] <- i
  cuts <- vector("list",length(available)); names(cuts) <- as.character(available)
  capture <- function(k) {
    lab <- match(group,unique(group))
    stopifnot(length(unique(lab))==k)
    cuts[[as.character(k)]] <<- lab
  }
  if (n %in% available) capture(n)
  if (m>0L) for (t in seq_len(m)) {
    child <- h$children[t,]
    node <- n+t
    if (anyNA(child)||any(child<1L|child>=node)||child[1]==child[2]||
        any(lengths(members[child])==0L)) stop("Invalid merge children.")
    idx <- c(members[[child[1]]],members[[child[2]]])
    group[idx] <- node
    members[child] <- list(integer(),integer())
    members[[node]] <- idx
    if ((n-t) %in% available) capture(n-t)
  }
  # Check the original reported cut independently of label numbering.
  if (h$actual_k %in% available &&
      !identical(as.integer(cuts[[as.character(h$actual_k)]]),as.integer(h$labels)))
    stop("Chronological replay does not reproduce the original labels.")
  if (any(lengths(cuts)!=n)) stop("Missing feasible hierarchy cut.")
  # No imputation: zero is used only in sums with an explicit contributor count.
  predict <- function(lab) {
    cnt <- rowsum(q+0,lab,reorder=TRUE)
    mu <- rowsum(X,lab,reorder=TRUE)/pmax(cnt,1)
    # On the scoring domain each mean has at least its own measured contributor.
    mu[lab,,drop=FALSE]
  }
  weighted_sum <- function(a) sum(colSums(a)*weights)
  profile <- data.frame(k=available, observed_mse=NA_real_,
                        observed_sse=NA_real_, observed_entries=sum(domain),
                        compression_ratio=n/available,
                        internal_boundary_edges=NA_integer_,
                        excess_validation_loss=NA_real_, upper_excess_loss=NA_real_)
  edges <- .four_neighbour_edges(final)
  reference_k <- max(available)
  have_validation <- !is.null(validation)
  if (have_validation) {
    if (!is.list(validation) || !isTRUE(validation$independent) ||
        !isTRUE(validation$unbiased) || !is.character(validation$source) ||
        length(validation$source)!=1L || is.na(validation$source) || !nzchar(validation$source) ||
        !identical(validation$representation_id,p$representation_id) ||
        !identical(as.integer(validation$selected_channel_indices),as.integer(p$selected_channel_indices)) ||
        !is.character(validation$error_model) || length(validation$error_model)!=1L ||
        !validation$error_model %in% c("gaussian","variance_only") ||
        !is.matrix(validation$features) || !is.numeric(validation$features) ||
        !is.matrix(validation$variance) || !is.numeric(validation$variance) ||
        !identical(dim(validation$features),dim(p$features)) ||
        !identical(dim(validation$variance),dim(p$features)))
      stop("Require independently measured, aligned, unbiased validation with marginal variance and provenance.")
    Y <- validation$features[original,,drop=FALSE]
    V <- validation$variance[original,,drop=FALSE]
    if (any(!is.finite(Y[domain]))||any(!is.finite(V[domain])|V[domain]<0))
      stop("Validation must cover the entire fixed scoring domain; no silent deletion.")
    Y[!domain] <- 0
    V[!domain] <- 0
    sd <- sqrt(V)
    base <- predict(cuts[[as.character(reference_k)]])
    # A one-sided union bound is simultaneous over every non-reference candidate.
    tail <- alpha/max(1,length(available)-1L)
    critical <- if (validation$error_model=="gaussian")
      stats::qnorm(tail,lower.tail=FALSE) else sqrt((1-tail)/tail)
  }
  for (i in seq_along(available)) {
    lab <- cuts[[i]]
    fit <- predict(lab)
    resid <- X-fit
    resid[!domain] <- 0
    sse <- weighted_sum(resid^2)
    if (!is.finite(sse)) stop("Nonfinite reconstruction loss.")
    profile$observed_sse[i] <- sse
    profile$observed_mse[i] <- sse/W
    profile$internal_boundary_edges[i] <- if(nrow(edges))sum(lab[edges[,1]]!=lab[edges[,2]]) else 0L
    if (have_validation) {
      # (Y-m)^2 - (Y-base)^2 = (base-m)*(2Y-m-base).
      # The contrast is linear in validation noise. No noise-floor subtraction
      # or independence between voxels is required for its expectation.
      a <- base-fit
      a[!domain] <- 0
      delta <- weighted_sum(a*(2*Y-fit-base))/W
      # By Cauchy--Schwarz, |Cov(e_i,e_j)| <= sd_i*sd_j, so this
      # is a worst-case standard deviation over all compatible covariances.
      bound_sd <- 2*weighted_sum(abs(a)*sd)/W
      upper <- delta+critical*bound_sd
      if (!is.finite(delta)||!is.finite(upper)) stop("Nonfinite validation contrast.")
      profile$excess_validation_loss[i] <- delta
      profile$upper_excess_loss[i] <- upper
    }
  }
  # Same-entry, same-weight squared reconstruction error cannot increase
  # under refinement. Check the exact unweighted SSE against Ward increments.
  if (all(weights==1)) {
    expected <- c(0,cumsum(h$cost))[n-available+1L]
    if (!isTRUE(all.equal(profile$observed_sse,expected,tolerance=1e-7,
                         check.attributes=FALSE)))
      stop("Observed-entry loss does not reproduce the frozen merge objective.")
  }
  coarse_loss <- profile$observed_mse[1]
  profile$relative_distortion <- if (coarse_loss>0)
    profile$observed_mse/coarse_loss else rep(0,nrow(profile))
  selected <- NA_integer_
  reason <- if (!have_validation) "NO_INDEPENDENT_VALIDATION" else
    if (is.null(tolerance)) "NO_SCIENTIFIC_TOLERANCE" else "NO_ADMISSIBLE_CANDIDATE"
  if (have_validation && !is.null(tolerance)) {
    ok <- profile$upper_excess_loss<=tolerance
    if (any(ok)) {
      selected <- min(profile$k[ok])
      reason <- if(selected==reference_k) "REFERENCE_ONLY_NO_SUPPORTED_COMPRESSION" else
        "CONDITIONAL_EXCESS_LOSS_BOUND"
    }
  }
  native_maps <- lapply(cuts,function(lab) {
    map <- matrix(NA_integer_,nrow(final),ncol(final))
    map[original] <- lab
    map
  })
  list(profile=profile,labels=native_maps,requested_k=candidate_k,
       unavailable_k=unavailable,coarsest_feasible_k=root_k,
       reference_k=reference_k,eligible_pixels=n,scoring_weight=W,
       selected_k=selected,selection_status=reason,tolerance=tolerance,alpha=alpha,
       representation_id=p$representation_id,representation=p$representation$type,
       coverage_qc=segmentation$coverage_qc,
       provenance=list(empirical="In-sample observed-entry weighted squared distortion; not information in bits",
         selection="Independent validation excess loss relative to finest evaluated cut; conditional on caller declarations",
         covariance="Worst-case covariance bound from marginal variances; no voxel independence assumption",
         validation_source=if(have_validation)validation$source else NULL,
         error_model=if(have_validation)validation$error_model else NULL,
         feature_weights=weights,
         hierarchy="Chronological cuts; no monotonicized heights or forced merges",
         physical_populations="Not inferred"))
}
