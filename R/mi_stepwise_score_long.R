# ==============================================================================
# Stepwise Score-Based Measurement Invariance for Repeated Cross-Sectional Data
# ==============================================================================
#
# This version implements the classic sequential strategy using lavTestScore:
#   1) Identify equality constraints that are most violated under metric invariance
#   2) Free the strongest one (largest score statistic / smallest p-value)
#   3) Refit the model
#   4) Re-check fit
#   5) Repeat until the metric model is acceptable
#
# This is done separately for each time point (across groups) and each group
# (across time points), and the resulting non-invariances are gathered into the
# final omnibus model across all G x T cells.
#
# NOTE: This implementation follows the general/global strategy, i.e., we use
# group.partial = c("factor =~ item", ...), not cell-specific labels.
# ==============================================================================

#' Stepwise measurement invariance for longitudinal multigroup data
#'
#' Identifies loading non-invariance across groups within each time point and
#' across time points within each group, then fits a final measurement model
#' across all group-time cells. This score-based implementation is the
#' measurement-invariance procedure used by [lmmgsem()].
#'
#' @param data Data frame containing the observed variables and grouping and
#'   time identifiers.
#' @param model Measurement model in lavaan syntax.
#' @param group_var Name of the grouping column in `data`.
#' @param time_var Name of the time column in `data`.
#' @param method Selection rule for freeing a loading:
#'   `"group_violations"` or `"test_statistic"`.
#' @param alpha Significance level for the score tests.
#' @param d_cfi_threshold Maximum configural-to-metric CFI difference.
#' @param max_iter Maximum number of iterations per group/time slice.
#' @param fit_indices Fit indices to return for the configural and metric
#'   models.
#' @param specific_noninv If `TRUE`, allow non-invariance for specific
#'   group-time pairs; otherwise free loadings globally.
#' @param ... Additional arguments passed to lavaan's CFA fitting functions.
#'
#' @details
#' Every group must have observations at each of at least two time points.
#' Missing group or time identifiers are rejected rather than silently dropped.
#'
#' @returns A list containing the final lavaan measurement model, latent
#'   covariance matrices for each group-time cell, invariance results, and fit
#'   indices.
#' @export
compute_mi_long_stepwise_score <- function(data,
                                           model,
                                           group_var,
                                           time_var,
                                           method = c("group_violations", "test_statistic"),
                                           alpha = 0.05,
                                           d_cfi_threshold = 0.010,
                                           max_iter = 20,
                                           fit_indices = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                                           specific_noninv = FALSE,
                                           ...) {

  method <- match.arg(method)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.character(model) || length(model) != 1L || !nzchar(model)) {
    stop("'model' must be a single non-empty lavaan model string.",
         call. = FALSE)
  }
  if (!is.character(group_var) || length(group_var) != 1L ||
      !group_var %in% names(data) ||
      !is.character(time_var) || length(time_var) != 1L ||
      !time_var %in% names(data)) {
    stop("'group_var' and 'time_var' must name columns in 'data'.",
         call. = FALSE)
  }
  if (identical(group_var, time_var)) {
    stop("'group_var' and 'time_var' must name different columns.",
         call. = FALSE)
  }
  if (anyNA(data[[group_var]]) || anyNA(data[[time_var]])) {
    stop("Group and time identifiers cannot contain missing values.",
         call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(d_cfi_threshold) != 1L ||
      !is.numeric(d_cfi_threshold) ||
      !is.finite(d_cfi_threshold) || d_cfi_threshold < 0) {
    stop("'d_cfi_threshold' must be a non-negative finite number.",
         call. = FALSE)
  }
  if (length(max_iter) != 1L || !is.numeric(max_iter) ||
      !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop("'max_iter' must be a positive integer.", call. = FALSE)
  }
  if (!is.character(fit_indices) || length(fit_indices) == 0L) {
    stop("'fit_indices' must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!is.logical(specific_noninv) || length(specific_noninv) != 1L ||
      is.na(specific_noninv)) {
    stop("'specific_noninv' must be TRUE or FALSE.", call. = FALSE)
  }

  # Normalize levels so covariance matrices follow the documented cell order.
  data[[group_var]] <- factor(
    as.character(data[[group_var]]),
    levels = unique(as.character(data[[group_var]]))
  )
  data[[time_var]] <- factor(
    as.character(data[[time_var]]),
    levels = unique(as.character(data[[time_var]]))
  )
  group_time_var <- ".lmmgsem_group_time"
  while (group_time_var %in% names(data)) {
    group_time_var <- paste0(group_time_var, "_")
  }
  data[[group_time_var]] <- interaction(
    data[[group_var]], data[[time_var]], sep = ".Time", drop = TRUE
  )

  unique_times  <- unique(as.character(data[[time_var]]))
  unique_groups <- unique(as.character(data[[group_var]]))
  all_cells     <- levels(data[[group_time_var]])

  if (length(unique_groups) == 0L || length(unique_times) < 2L) {
    stop("At least one group and two time points are required.",
         call. = FALSE)
  }
  if (any(table(data[[group_var]], data[[time_var]]) == 0L)) {
    stop("Every group must have observations at every time point.",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------
  build_metric_model <- function(sub_data,
                                 split_var,
                                 partial_constraints = character(0),
                                 ...) {
    metric_syn <- semTools::measEq.syntax(
      configural.model = model,
      data = sub_data,
      group = split_var,
      group.equal = "loadings",
      group.partial = if (length(partial_constraints) > 0) partial_constraints else NULL,
      std.lv = TRUE
    )

    metric_fit <- lavaan::cfa(
      model = as.character(metric_syn),
      data = sub_data,
      group = split_var,
      ...
    )

    return(metric_fit)
  }

  extract_score_pairs <- function(fit, fixed_val, mode = c("wave", "group"), alpha = 0.05) {
    mode <- match.arg(mode)

    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(list(global = character(0), pairs = NULL))
    }

    ts <- tryCatch(lavaan::lavTestScore(fit), error = function(e) NULL)
    if (is.null(ts) || is.null(ts$uni)) {
      return(list(global = character(0), pairs = NULL))
    }

    uni_tab <- ts$uni
    uni_tab <- uni_tab[uni_tab$op == "==", , drop = FALSE]
    if (nrow(uni_tab) == 0) {
      return(list(global = character(0), pairs = NULL))
    }

    uni_tab <- uni_tab[!is.na(uni_tab$p.value), , drop = FALSE]
    if (nrow(uni_tab) == 0) {
      return(list(global = character(0), pairs = NULL))
    }

    corrected_alpha <- alpha / nrow(uni_tab)
    sig_tests <- uni_tab[uni_tab$p.value < corrected_alpha, , drop = FALSE]
    if (nrow(sig_tests) == 0) {
      return(list(global = character(0), pairs = NULL))
    }

    pt <- lavaan::parTable(fit)
    grp_labels <- lavaan::lavInspect(fit, "group.label")

    glob_out <- character(0)
    pair_out <- list()

    for (i in seq_len(nrow(sig_tests))) {
      p1 <- sig_tests$lhs[i]
      p2 <- sig_tests$rhs[i]

      row1 <- pt[pt$plabel == p1 | pt$label == p1, , drop = FALSE]
      row2 <- pt[pt$plabel == p2 | pt$label == p2, , drop = FALSE]

      if (nrow(row1) > 0 && nrow(row2) > 0 && row1$op[1] == "=~") {
        l_fac <- row1$lhs[1]
        l_item <- row1$rhs[1]
        glob_out <- c(glob_out, paste(l_fac, "=~", l_item))

        g1_name <- grp_labels[row1$group[1]]
        g2_name <- grp_labels[row2$group[1]]

        if (mode == "wave") {
          c1 <- paste0(g1_name, ".Time", fixed_val)
          c2 <- paste0(g2_name, ".Time", fixed_val)
        } else {
          c1 <- paste0(fixed_val, ".Time", g1_name)
          c2 <- paste0(fixed_val, ".Time", g2_name)
        }

        pair_out[[length(pair_out) + 1]] <- data.frame(
          lhs = l_fac,
          rhs = l_item,
          cell1 = c1,
          cell2 = c2,
          stringsAsFactors = FALSE
        )
      }
    }

    p_df <- if (length(pair_out) > 0) do.call(rbind, pair_out) else NULL
    return(list(global = unique(glob_out), pairs = p_df))
  }

  build_specific_noninv_syntax <- function(base_model, cell_labels, flagged_pairs) {
    parsed_model <- lavaan::lavaanify(base_model)
    load_rows <- parsed_model[parsed_model$op == "=~", , drop = FALSE]
    syntax_lines <- character(0)

    for (fac in unique(load_rows$lhs)) {
      fac_items <- load_rows$rhs[load_rows$lhs == fac]
      fac_elements <- character(0)

      for (i in seq_along(fac_items)) {
        itm <- fac_items[i]

        if (i == 1) {
          fac_elements <- c(fac_elements, itm)
          next
        }

        item_flags <- flagged_pairs[flagged_pairs$lhs == fac & flagged_pairs$rhs == itm, , drop = FALSE]

        lbl_vec <- character(length(cell_labels))
        for (c_idx in seq_along(cell_labels)) {
          c_name <- cell_labels[c_idx]
          in_violation <- any(item_flags$cell1 == c_name | item_flags$cell2 == c_name)

          if (in_violation) {
            lbl_vec[c_idx] <- paste0("l_", fac, "_", itm, "_", gsub("[^[:alnum:]]", "_", c_name))
          } else {
            lbl_vec[c_idx] <- paste0("l_", fac, "_", itm, "_common")
          }
        }

        lbl_string <- paste0("c(", paste(lbl_vec, collapse = ", "), ")")
        fac_elements <- c(fac_elements, paste0(lbl_string, "*", itm))
      }

      syntax_lines <- c(syntax_lines, paste0(fac, " =~ ", paste(fac_elements, collapse = " + ")))
    }

    paste(syntax_lines, collapse = "\n")
  }

  build_specific_metric_model <- function(sub_data, split_var, pair_flags, ...) {
    if (nrow(pair_flags) == 0) {
      return(build_metric_model(sub_data, split_var, character(0), ...))
    }

    cell_labels <- unique(as.character(sub_data[[split_var]]))
    custom_syntax <- build_specific_noninv_syntax(model, cell_labels, pair_flags)

    lavaan::cfa(
      model = custom_syntax,
      data = sub_data,
      group = split_var,
      ...
    )
  }

  # Extract the most violated invariance constraint based on how often it is
  # implicated across the group-vs-reference comparisons evaluated by lavTestScore.
  candidate_from_score <- function(fit, alpha = 0.05, selection_method = "group_violations") {
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(NULL)
    }

    ts <- tryCatch(lavaan::lavTestScore(fit), error = function(e) NULL)
    if (is.null(ts) || is.null(ts$uni)) {
      return(NULL)
    }

    uni_tab <- ts$uni                                      # Extract the univariate tests
    uni_tab <- uni_tab[uni_tab$op == "==", , drop = FALSE] # Extract only the equality constraints
    if (nrow(uni_tab) == 0) {
      return(NULL)
    }

    uni_tab <- uni_tab[!is.na(uni_tab$p.value), , drop = FALSE]
    if (nrow(uni_tab) == 0) {
      return(NULL)
    }

    # Bonferroni-adjusted significance threshold for deciding whether a given
    # comparison is a violation.
    alpha_adj <- alpha / nrow(uni_tab)
    uni_tab$constraint_key <- vapply(
      seq_len(nrow(uni_tab)),
      function(i) constraint_key_from_score_row(fit, uni_tab$lhs[i], uni_tab$rhs[i]),
      character(1)
    )
    uni_tab$violated <- as.logical(uni_tab$p.value < alpha_adj)

    if (!any(uni_tab$violated)) {
      return(NULL)
    }

    if (selection_method == "test_statistic") {
      # Filter to Bonferroni violations, sort descending by univariate X2 score
      viol_tab <- uni_tab[uni_tab$violated, , drop = FALSE]
      viol_tab <- viol_tab[order(-viol_tab$X2, viol_tab$p.value), , drop = FALSE]
      return(viol_tab[1, , drop = FALSE])
    } else {
      # group_violations: Pick item with the most pairwise violations across slices
      agg <- aggregate(
        violated ~ constraint_key,
        data = uni_tab,
        FUN = sum,
        na.action = na.omit
      )
      names(agg)[2] <- "n_violations"

      agg <- agg[order(-agg$n_violations), , drop = FALSE]
      chosen_key <- agg$constraint_key[1]

      candidate_rows <- uni_tab[uni_tab$constraint_key == chosen_key, , drop = FALSE]
      candidate_rows <- candidate_rows[order(-candidate_rows$X2, candidate_rows$p.value), , drop = FALSE]

      return(candidate_rows[1, , drop = FALSE])
    }
  }

  parameter_name_from_constraint <- function(fit, lhs_label, rhs_label) {
    pt <- lavaan::parTable(fit)

    # lavTestScore equality rows are often parameter labels, so we map back to
    # a loading expression such as "eta1 =~ y1" using the partable.
    row1 <- pt[pt$plabel == lhs_label | pt$label == lhs_label, , drop = FALSE]
    row2 <- pt[pt$plabel == rhs_label | pt$label == rhs_label, , drop = FALSE] # In case we add the group-specific version later

    if (nrow(row1) == 0 || nrow(row2) == 0) {
      return(NULL)
    }

    # Keep the row that corresponds to a loading constraint
    idx <- which(row1$op == "=~")
    if (length(idx) == 0) {
      return(NULL)
    }

    r1 <- row1[idx[1], ]
    return(paste(r1$lhs, "=~", r1$rhs))
  }

  constraint_key_from_score_row <- function(fit, lhs_label, rhs_label) {
    # lavTestScore often returns group-specific labels (e.g., a parameter label that
    # differs by group) even when they represent the same underlying loading.
    # We want the stable semantic key: the item loading itself, not the group-specific
    # label used in that contrast.
    pt <- lavaan::parTable(fit)

    row1 <- pt[pt$label == lhs_label | pt$plabel == lhs_label, , drop = FALSE]
    row2 <- pt[pt$label == rhs_label | pt$plabel == rhs_label, , drop = FALSE]

    if (nrow(row1) == 0 || nrow(row2) == 0) {
      return(NULL)
    }

    load_rows <- rbind(row1[row1$op == "=~", , drop = FALSE], row2[row2$op == "=~", , drop = FALSE])
    if (nrow(load_rows) == 0) {
      return(NULL)
    }

    r1 <- load_rows[1, ]
    paste(r1$lhs, "=~", r1$rhs)
  }

  stepwise_slice <- function(sub_data,
                             split_var,
                             label,
                             selection_method,
                             specific_noninv = FALSE,
                             ...) {
    config_syn <- semTools::measEq.syntax(
      configural.model = model,
      data = sub_data,
      group = split_var,
      std.lv = TRUE
    )
    config_fit <- lavaan::cfa(
      model = as.character(config_syn),
      data = sub_data,
      group = split_var,
      ...
    )

    config_cfi <- lavaan::fitMeasures(config_fit, "cfi")
    if (is.null(config_cfi) || is.na(config_cfi)) {
      config_cfi <- NA_real_
    }

    partial_constraints <- character(0)
    history <- list()
    pair_flags <- data.frame(
      lhs = character(0),
      rhs = character(0),
      cell1 = character(0),
      cell2 = character(0),
      stringsAsFactors = FALSE
    )

    if (specific_noninv) {
      current_fit <- build_specific_metric_model(sub_data, split_var, pair_flags, ...)
    } else {
      current_fit <- build_metric_model(sub_data, split_var, partial_constraints, ...)
    }

    iter <- 0
    delta_cfi <- Inf

    while (iter < max_iter && (is.na(delta_cfi) || delta_cfi > d_cfi_threshold)) {
      iter <- iter + 1

      metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
      if (is.null(metric_cfi) || is.na(metric_cfi)) {
        break
      }

      candidate <- candidate_from_score(current_fit, alpha = alpha, selection_method = method)
      if (is.null(candidate)) {
        break
      }

      if (specific_noninv) {
        if (!"lhs" %in% names(candidate) || !"rhs" %in% names(candidate)) {
          break
        }

        pt <- lavaan::parTable(current_fit)
        row1 <- pt[pt$plabel == candidate$lhs | pt$label == candidate$lhs, , drop = FALSE]
        row2 <- pt[pt$plabel == candidate$rhs | pt$label == candidate$rhs, , drop = FALSE]

        if (nrow(row1) == 0 || nrow(row2) == 0) {
          break
        }

        row1_load <- row1[row1$op == "=~", , drop = FALSE]
        row2_load <- row2[row2$op == "=~", , drop = FALSE]
        if (nrow(row1_load) == 0 || nrow(row2_load) == 0) {
          break
        }

        group1 <- as.character(unique(row1_load$group))
        group2 <- as.character(unique(row2_load$group))
        candidate_key <- paste0(row1_load$lhs[1], "=~", row1_load$rhs[1], "|", group1, "-", group2)

        if (length(group1) == 0 || length(group2) == 0 || group1 == group2) {
          break
        }

        same_pair <- any(pair_flags$lhs == row1_load$lhs[1] & pair_flags$rhs == row1_load$rhs[1] &
                           ((pair_flags$cell1 == group1 & pair_flags$cell2 == group2) |
                              (pair_flags$cell1 == group2 & pair_flags$cell2 == group1)))
        if (same_pair) {
          break
        }

        pair_flags <- rbind(
          pair_flags,
          data.frame(
            lhs = row1_load$lhs[1],
            rhs = row1_load$rhs[1],
            cell1 = group1,
            cell2 = group2,
            stringsAsFactors = FALSE
          )
        )

        history[[length(history) + 1]] <- list(
          iter = iter,
          parameter = paste(row1_load$lhs[1], "=~", row1_load$rhs[1]),
          groups = c(group1, group2),
          stat = candidate$X2,
          p_value = candidate$p.value
        )

        current_fit <- build_specific_metric_model(sub_data, split_var, pair_flags, ...)
      } else {
        param_name <- parameter_name_from_constraint(current_fit, candidate$lhs, candidate$rhs)
        if (is.null(param_name) || param_name %in% partial_constraints) {
          break
        }

        delta_cfi <- if (is.na(config_cfi)) NA_real_ else config_cfi - metric_cfi

        partial_constraints <- c(partial_constraints, param_name)
        history[[length(history) + 1]] <- list(
          iter = iter,
          parameter = param_name,
          stat = candidate$X2,
          p_value = candidate$p.value
        )

        current_fit <- build_metric_model(sub_data, split_var, partial_constraints, ...)
      }

      metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
      delta_cfi <- if (is.na(config_cfi)) NA_real_ else config_cfi - metric_cfi
    }

    final_metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
    final_delta_cfi <- if (is.na(config_cfi) || is.na(final_metric_cfi)) NA_real_ else config_cfi - final_metric_cfi

    return(list(
      label = label,
      config_fit = config_fit,
      metric_fit = current_fit,
      partial_constraints = partial_constraints,
      pair_flags = pair_flags,
      history = history,
      config_cfi = config_cfi,
      metric_cfi = final_metric_cfi,
      delta_cfi = final_delta_cfi,
      acceptable = if (is.na(final_delta_cfi)) FALSE else final_delta_cfi <= d_cfi_threshold
    ))
  }

  # Helper to rewrite syntax with an invariant marker variable per factor
  adjust_marker_variables <- function(base_model, non_invariant_loadings) {
    pt <- lavaan::lavParseModelString(base_model, as.data.frame = TRUE)
    load_pt <- pt[pt$op == "=~", , drop = FALSE]
    factors <- unique(load_pt$lhs)

    clean_key <- function(lhs, rhs) gsub("\\s+", "", paste(lhs, "=~", rhs))
    freed_keys <- gsub("\\s+", "", non_invariant_loadings)

    new_syntax_lines <- character(0)

    for (fac in factors) {
      fac_items <- load_pt$rhs[load_pt$lhs == fac]
      fac_keys  <- clean_key(fac, fac_items)

      # Determine which items for this factor remained invariant
      is_invariant <- !(fac_keys %in% freed_keys)

      if (any(is_invariant)) {
        # Select the first invariant item as the new marker
        marker_idx <- which(is_invariant)[1]
      } else {
        # Fallback: if no invariant item exists on this factor, keep the first item
        marker_idx <- 1
      }

      marker_item <- fac_items[marker_idx]
      other_items <- fac_items[-marker_idx]

      # Reconstruct factor definition: marker fixed to 1, others estimated freely (NA*)
      if (length(other_items) > 0) {
        line <- sprintf("%s =~ 1*%s + %s", fac, marker_item, paste0(other_items, collapse = " + "))
      } else {
        line <- sprintf("%s =~ 1*%s", fac, marker_item)
      }
      new_syntax_lines <- c(new_syntax_lines, line)
    }

    # Retain any regressions, covariances, or intercepts from the original syntax
    # other_pt <- pt[pt$op != "=~", , drop = FALSE]
    # if (nrow(other_pt) > 0) {
    #   other_lines <- vapply(seq_len(nrow(other_pt)), function(i) {
    #     paste(other_pt$lhs[i], other_pt$op[i], other_pt$rhs[i])
    #   }, character(1))
    #   new_syntax_lines <- c(new_syntax_lines, other_lines)
    # }

    paste(new_syntax_lines, collapse = "\n")
  }

  # ---------------------------------------------------------------------------
  # Step 1: Across groups within each time point
  # ---------------------------------------------------------------------------
  flagged_global <- character(0)
  flagged_pairs <- data.frame(
    lhs = character(0),
    rhs = character(0),
    cell1 = character(0),
    cell2 = character(0),
    stringsAsFactors = FALSE
  )

  message("Executing Step 1: cross-sectional metric invariance, freeing the strongest violations sequentially...")
  cross_sectional_results <- vector("list", length(unique_times))
  names(cross_sectional_results) <- paste0("Time_", unique_times)

  pb_1 <- utils::txtProgressBar(min = 0, max = length(unique_times), style = 3)
  for (t in seq_along(unique_times)) {
    sub_data <- data[data[[time_var]] == unique_times[t], , drop = FALSE]
    result <- stepwise_slice(sub_data, split_var = group_var, label = unique_times[t], selection_method = method, specific_noninv = specific_noninv, ...)

    if (specific_noninv) {
      if (nrow(result$pair_flags) > 0) {
        pair_details <- result$pair_flags
        pair_details$cell1 <- paste0(pair_details$cell1, ".Time", unique_times[t])
        pair_details$cell2 <- paste0(pair_details$cell2, ".Time", unique_times[t])
        flagged_pairs <- rbind(flagged_pairs, pair_details)
      }
    }

    cross_sectional_results[[t]] <- result
    utils::setTxtProgressBar(pb_1, t)
  }
  close(pb_1)

  # ---------------------------------------------------------------------------
  # Step 2: Across time points within each group
  # ---------------------------------------------------------------------------
  message("Executing Step 2: longitudinal metric invariance, freeing the strongest violations sequentially...")
  longitudinal_results <- vector("list", length(unique_groups))
  names(longitudinal_results) <- paste0("Group_", unique_groups)

  pb_2 <- utils::txtProgressBar(min = 0, max = length(unique_groups), style = 3)
  for (g in seq_along(unique_groups)) {
    sub_data <- data[data[[group_var]] == unique_groups[g], , drop = FALSE]
    result <- stepwise_slice(sub_data, split_var = time_var, label = unique_groups[g], selection_method = method, specific_noninv = specific_noninv, ...)

    if (specific_noninv) {
      if (nrow(result$pair_flags) > 0) {
        pair_details <- result$pair_flags
        pair_details$cell1 <- paste0(unique_groups[g], ".Time", pair_details$cell1)
        pair_details$cell2 <- paste0(unique_groups[g], ".Time", pair_details$cell2)
        flagged_pairs <- rbind(flagged_pairs, pair_details)
      }
    }

    longitudinal_results[[g]] <- result
    utils::setTxtProgressBar(pb_2, g)
  }
  close(pb_2)

  # ---------------------------------------------------------------------------
  # Collect the final partial invariance set (global/general strategy)
  # ---------------------------------------------------------------------------
  final_partials <- character(0)

  for (res in cross_sectional_results) {
    final_partials <- c(final_partials, res$partial_constraints)
  }
  for (res in longitudinal_results) {
    final_partials <- c(final_partials, res$partial_constraints)
  }

  final_partials <- unique(final_partials)

  # ---------------------------------------------------------------------------
  # Final omnibus model across all G x T cells
  # ---------------------------------------------------------------------------
  message("Estimating final omnibus model across all group-by-time cells...")

  if (specific_noninv) {
    final_omnibus_model_syntax <- build_specific_noninv_syntax(model, all_cells, flagged_pairs)

    final_fit <- lavaan::cfa(
      model = final_omnibus_model_syntax,
      data = data,
      group = group_time_var,
      ...
    )
  } else {
    final_omnibus_model_syntax <- adjust_marker_variables(model, final_partials)

    final_fit <- lavaan::cfa(
      model         = final_omnibus_model_syntax,
      data          = data,
      group         = group_time_var,
      group.equal   = "loadings",
      group.partial = if (length(final_partials) > 0) final_partials else NULL,
      std.lv        = FALSE,
      ...
    )
  }

  phi_matrices <- lavaan::lavInspect(final_fit, "cov.lv")

  fit_measures_list <- list(
    step1_across_groups = lapply(cross_sectional_results, function(x) {
      rbind(
        configural = lavaan::fitMeasures(x$config_fit, fit_indices),
        metric = lavaan::fitMeasures(x$metric_fit, fit_indices)
      )
    }),
    step2_across_time = lapply(longitudinal_results, function(x) {
      rbind(
        configural = lavaan::fitMeasures(x$config_fit, fit_indices),
        metric = lavaan::fitMeasures(x$metric_fit, fit_indices)
      )
    }),
    final_model = lavaan::fitMeasures(final_fit, fit_indices)
  )

  return(list(
    step1_across_groups = cross_sectional_results,
    step2_across_time = longitudinal_results,
    group_levels = unique_groups,
    time_levels = unique_times,
    final_fit = final_fit,
    final_syntax = final_omnibus_model_syntax,
    phi_matrices = phi_matrices,
    final_partial_constraints = final_partials,
    flagged_global = flagged_global,
    flagged_pairs = flagged_pairs,
    fit_measures = fit_measures_list
  ))
}
