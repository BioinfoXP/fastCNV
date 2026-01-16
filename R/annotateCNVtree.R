#' Annotate a Phylogenetic Tree with CNV Events
#'
#' This function annotates a phylogenetic tree with copy number variation (CNV) events. It identifies
#' significant CNV events in the provided matrix, links them to clones and ancestral nodes, and
#' updates the tree with this information.
#'
#' @param tree A phylogenetic tree (of class \code{phylo}) that will be annotated.
#' @param cnv_mat A matrix of copy number variation (CNV) values, with samples as rows and regions as columns.
#' @param cnv_thresh A numeric threshold to filter significant CNV events. Default is 0.15.
#'
#' @return A data frame with the tree data, including annotations for CNV events.
#'
#' @examples
#' cnv_matrix <- structure(c(0.2, 0.4, 0, 0, 0.1, 0, 0.1, 0.2, 0.2), dim = c(3L,
#' 3L), dimnames = list(c("Clone 1", "Clone 2", "Clone 3"), c("Region 1",
#'                                                           "Region 2", "Region 3")))
#' tree <- buildCNVTree(cnv_matrix)
#' tree_data <- annotateCNVTree(tree, cnv_matrix)
#'
#' @importFrom dplyr left_join full_join pull
#' @importFrom ggplot2 fortify
#' @import stringr
#' @import dplyr
#' @importFrom utils combn
#' @importFrom purrr map2
#'
#' @export
#'
annotateCNVTree <- function(tree, cnv_mat, cnv_thresh = 0.15) {
  tree_data <- fortify(tree)
  # Get events
  major_events <- get_majorEvents(cnv_mat, cnv_thresh)
  clone_events <- get_cloneEvents(major_events)
  node_events <- get_ancestralEvents(tree, cnv_mat, major_events)
  # Annotate tree
  if(nrow(clone_events) < 1) {
    return(tree_data)
  }
  if(nrow(node_events) > 0) {
    event_data <- full_join(clone_events, node_events, by = c("label", "node", "all_events"))
  } else {
    event_data <- clone_events
  }
  tree_data <- left_join(tree_data, event_data, by = c("node", "label"))
  tree_data$events <- sapply(1:nrow(tree_data), .remove_parent_events, tree_data)
  tree_data$events <- sapply(tree_data$events, .aggregate_events)
  tree_data$events[tree_data$events == ""] <- NA_character_
  return(tree_data)
}


get_majorEvents <- function(cnv_mat, thresh) {
  gains <- .get_gains(cnv_mat, thresh)
  losses <- .get_losses(cnv_mat, thresh)
  major_events <- lapply(1:length(gains), function(i){
    events <- c(gains[[i]], losses[[i]])
    str_subset(events, "\\d|X") |> str_remove("\\.")
  })
  major_events <- .order_events(major_events)
  names(major_events) <- rownames(cnv_mat)
  return(major_events)
}

.get_gains <- function(cnv_mat, thresh) {
  gain_events <- apply(cnv_mat > thresh, 1, which)
  if(length(gain_events) == 0) {
    gain_events <- lapply(1:nrow(cnv_mat), function(i) { vector() })
  }
  clone_gains <- lapply(gain_events, function(events) {
    arms <- str_remove(names(events), "_CNV")
    paste0("+",arms)
  })
  return(clone_gains)
}
.get_losses <- function(cnv_mat, thresh) {
  loss_events <- apply(cnv_mat < -thresh, 1, which)
  if(length(loss_events) == 0) {
    loss_events <- lapply(1:nrow(cnv_mat), function(i) { vector() })
  }
  clone_losses <- lapply(loss_events, function(events) {
    arms <- str_remove(names(events), "_CNV")
    paste0("-",arms)
  })
  return(clone_losses)
}

.order_events <- function(all_events) {
  ordered_events <- lapply(all_events, function(events) {
    chrs <- str_extract(events, "[0-9]+|X") |> str_replace("X", "23") |> as.numeric()
    events[order(chrs)]
  })
  names(all_events) <- names(all_events)
  return(ordered_events)
}

get_cloneEvents <- function(major_events) {
  clone_events <- sapply(major_events, paste0, collapse = " ")
  clone_events <-
    tibble(label = names(clone_events),
           node = 1:length(clone_events),
           all_events = clone_events)

  return(clone_events)
}

.remove_parent_events <- function(i, tree_data) {
  child_events <- str_split(tree_data$all_events[i], " ")[[1]]
  parent <- tree_data$parent[i]
  parent_events <-  str_split(tree_data$all_events[parent], " ")[[1]]
  remaining_events <- setdiff(child_events, parent_events)
  remaining_events <- paste0(remaining_events, collapse = " ")
  return(remaining_events)
}
.aggregate_events <- function(events) {
  events <- str_split(events, " ")[[1]]
  chr_chrs <- str_remove(events, "p|q")
  chr_events <- chr_chrs[duplicated(chr_chrs)]
  arm_events <- events[ave(chr_chrs, chr_chrs, FUN = length) == 1]
  all_events <- c(chr_events, arm_events)
  all_events <- .order_events(list(all_events))[[1]]
  agg_events <- paste0(all_events, collapse = " ")
  return(agg_events)
}


get_ancestralEvents <- function(tree, cnv_mat, major_events) {
  cnv_combns <- combn(rownames(cnv_mat), 2) |>
    t() |>
    as.data.frame() |>
    tibble() |>
    mutate(
      shared_events = map2(.data$V1, .data$V2, function(c1,c2) {
        intersect(major_events[[c1]], major_events[[c2]])
      }))

  cnv_combns$V1_indices <- match(cnv_combns$V1, tree$tip.label)
  cnv_combns$V2_indices <- match(cnv_combns$V2, tree$tip.label)

  cnv_combns$node <- sapply(1:nrow(cnv_combns), function(i) {
    intersect(Ancestors(tree, cnv_combns$V1_indices[i]),
              Ancestors(tree, cnv_combns$V2_indices[i]))[1]
  })
  cnv_parent_events <- cnv_combns |>
    group_by(.data$node) |>
    reframe(all_events = if(length(.data$shared_events) > 1)
      Reduce(intersect, .data$shared_events)
      else unlist(.data$shared_events)) |>
    group_by(.data$node) |>
    mutate(all_events = paste0(.data$all_events, collapse = " "),
           label = NA) |>
    distinct() |>
    ungroup()

  cnv_parent_events$all_events <- sapply(
    str_split(cnv_parent_events$all_events, " "), function(events) {
      .order_events(events) |> paste(collapse = " ")})

  return(cnv_parent_events)
}
