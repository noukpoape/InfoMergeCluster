# conversion of the data to the internal long format (id, t, y), by column
# position:
#   - long data.frame (or matrix): 1st column = subject identifier,
#     second-to-last = time, last = response;
#   - list of data.frames or matrices (one element per subject):
#     second-to-last column = time, last = response; the identifier is the
#     element name if all names are non-empty and distinct, otherwise its rank.
.im_num <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  suppressWarnings(as.numeric(v))
}

#   n_cov > 0 (columns identified by position):
#     - long data.frame: 1st = identifier, 2nd = time, 3rd = response,
#       columns 4..(3 + n_cov) = covariates;
#     - list: 1st column = time, 2nd = response, columns 3..(2 + n_cov) =
#       covariates.
#   The covariates are returned as columns z1, ..., z<n_cov>, with their type.
.im_as_long <- function(data, n_cov = 0L) {
  zn <- if (n_cov > 0) paste0("z", seq_len(n_cov)) else character(0)
  grab <- function(e, cols) {
    cv <- e[cols]
    names(cv) <- zn
    cv
  }
  if (is.data.frame(data) || is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    nc <- ncol(data)
    if (nc < 3)
      stop("A long data.frame must have at least 3 columns: identifier (1st), ",
           "time (second-to-last) and response (last).")
    id <- data[[1]]
    if (is.factor(id)) id <- as.character(id)
    if (n_cov > 0) {
      if (nc < 3L + n_cov)
        stop(sprintf("With n_cov = %d, a long data.frame must have at least %d columns: identifier, time, response and the covariates.", n_cov, 3L + n_cov))
      out <- data.frame(id = id, t = .im_num(data[[2]]), y = .im_num(data[[3]]),
                        stringsAsFactors = FALSE)
      return(cbind(out, grab(data, 3L + seq_len(n_cov))))
    }
    return(data.frame(id = id, t = .im_num(data[[nc - 1]]), y = .im_num(data[[nc]]),
                      stringsAsFactors = FALSE))
  }
  if (is.list(data)) {
    if (length(data) == 0) stop("The list 'data' is empty.")
    ok <- vapply(data, function(e) is.data.frame(e) || is.matrix(e), logical(1))
    if (!all(ok))
      stop("Each element of the list must be a data.frame or a matrix ",
           "(offending elements: ", paste(which(!ok), collapse = ", "), ").")
    nm <- names(data)
    ids <- if (!is.null(nm) && all(nzchar(nm)) && !anyDuplicated(nm)) nm else seq_along(data)
    parts <- lapply(seq_along(data), function(i) {
      e <- as.data.frame(data[[i]], stringsAsFactors = FALSE)
      nc <- ncol(e)
      if (nc < 2)
        stop(sprintf("Element %d must have at least 2 columns: time (second-to-last) and response (last).", i))
      if (nrow(e) == 0) return(NULL)
      if (n_cov > 0) {
        if (nc < 2L + n_cov)
          stop(sprintf("With n_cov = %d, element %d must have at least %d columns: time, response and the covariates.", n_cov, i, 2L + n_cov))
        out_i <- data.frame(id = ids[i], t = .im_num(e[[1]]), y = .im_num(e[[2]]),
                            stringsAsFactors = FALSE)
        return(cbind(out_i, grab(e, 2L + seq_len(n_cov))))
      }
      data.frame(id = ids[i], t = .im_num(e[[nc - 1]]), y = .im_num(e[[nc]]),
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, parts)
    if (is.null(out)) stop("No measurements in 'data'.")
    return(out)
  }
  stop("'data' must be a data.frame, a list of data.frames or a list of matrices.")
}

# numeric design of the covariates (dummy coding, centering) on the kept rows
.im_covariate_design <- function(cov_raw, ok) {
  if (is.null(cov_raw)) return(list(X = NULL, names = character(0), center = numeric(0)))
  cv <- cov_raw[ok, , drop = FALSE]
  cv[] <- lapply(cv, function(v) if (is.character(v)) factor(v) else v)
  cv <- droplevels(cv)
  X <- stats::model.matrix(~ ., data = cv)[, -1, drop = FALSE]
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  keep <- apply(X, 2, function(v) stats::var(v) > 0)
  if (any(!keep))
    warning("Covariate column(s) without variation removed: ",
            paste(colnames(X)[!keep], collapse = ", "), call. = FALSE)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) return(list(X = NULL, names = character(0), center = numeric(0)))
  center <- colMeans(X)
  X <- sweep(X, 2, center)
  nm <- paste0("cov_", colnames(X))
  colnames(X) <- nm
  list(X = as.data.frame(X), names = nm, center = stats::setNames(center, nm))
}
