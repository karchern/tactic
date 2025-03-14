load_yaml_to_global <- function(file_path) {
  # Read the YAML file
  yaml_data <- yaml::read_yaml(file_path)

  # Assign each key-value pair to the global environment
  for (name in names(yaml_data)) {
    assign(name, yaml_data[[name]], envir = .GlobalEnv)
  }
}


read96wMaps <- function(folder = "./input/maps/") {
  lapply(list.files(folder, pattern = "library_plate", full.names = T), fread) %>%
    rbindlist()
}


registerQuadrants <- function(
    form = "384w", plate_number = NULL,
    tech_rep = NULL) {
  if (!form %in% c("384w", "1536w")) stop("\nUnknown plate format. Should be '384w' or '1536'")
  if (is.null(plate_number)) stop("\nWas expecting plate numbers!")
  if (length(plate_number) %% 4 != 0) stop("\nWas expecting plates for all four quadrants!")

  if (is.null(tech_rep)) tech_rep <- rep(1, length(plate_number))

  if (length(plate_number) != length(tech_rep)) stop("\nNumber of techreps does not match the number of plates!")

  n <- length(plate_number) / 4
  mat <- matrix(plate_number, nrow = 4)
  mat <- cbind(1:4, mat)

  plates <- paste0("plt_", seq_len(n))
  colnames(mat) <- c("qrt", plates)

  out <- as.data.frame(mat) %>%
    pivot_longer(!contains("qrt"), names_to = "to", values_to = "from") %>%
    arrange(to, qrt) %>%
    # tidy up
    mutate(to = gsub("plt_", "", to)) %>%
    mutate(across(everything(), as.numeric))

  out <- cbind(out, tech_rep)

  if (form == "384w") {
    names(out) <- c("qrt384", "plt384", "plt96", "techrep96")
    if (length(unique(out$plt384)) == 1) out <- subset(out, select = -plt384)
  } else {
    names(out) <- c("qrt1536", "plt1536", "plt384", "techrep384")
    if (length(unique(out$plt1536)) == 1) out <- subset(out, select = -plt1536)
  }

  out
}



add384RowAndCol <- function(dat) {
  # add 384w coordinates
  rowcol384 <- rbind(
    data.table(
      qrt384 = 1,
      row96 = rep(LETTERS[1:8], each = 12),
      col96 = rep(1:12, 8),
      row384 = rep(LETTERS[seq(1, 16, 2)], each = 12),
      col384 = rep(seq(1, 24, 2), 8)
    ),
    data.table(
      qrt384 = 2,
      row96 = rep(LETTERS[1:8], each = 12),
      col96 = rep(1:12, 8),
      row384 = rep(LETTERS[seq(1, 16, 2)], each = 12),
      col384 = rep(seq(2, 24, 2), 8)
    ),
    data.table(
      qrt384 = 3,
      row96 = rep(LETTERS[1:8], each = 12),
      col96 = rep(1:12, 8),
      row384 = rep(LETTERS[seq(2, 16, 2)], each = 12),
      col384 = rep(seq(1, 24, 2), 8)
    ),
    data.table(
      qrt384 = 4,
      row96 = rep(LETTERS[1:8], each = 12),
      col96 = rep(1:12, 8),
      row384 = rep(LETTERS[seq(2, 16, 2)], each = 12),
      col384 = rep(seq(2, 24, 2), 8)
    )
  )
  left_join(dat, rowcol384, by = c("qrt384"), relationship = "many-to-many")
}

add1536RowAndCol <- function(dat) {
  # add 1536w coordinates
  rowcol1536 <- rbind(
    data.table(
      qrt1536  = 1,
      row384   = rep(LETTERS[1:16], each = 24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(1, 32, 2), each = 24),
      col1536  = rep(seq(1, 48, 2), 16)
    ),
    data.table(
      qrt1536  = 2,
      row384   = rep(LETTERS[1:16], each = 24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(1, 32, 2), each = 24),
      col1536  = rep(seq(2, 48, 2), 16)
    ),
    data.table(
      qrt1536 = 3,
      row384 = rep(LETTERS[1:16], each = 24),
      col384 = rep(seq(1, 24), 16),
      row1536 = rep(seq(2, 32, 2), each = 24),
      col1536 = rep(seq(1, 48, 2), 16)
    ),
    data.table(
      qrt1536 = 4,
      row384 = rep(LETTERS[1:16], each = 24),
      col384 = rep(seq(1, 24), 16),
      row1536 = rep(seq(2, 32, 2), each = 24),
      col1536 = rep(seq(2, 48, 2), 16)
    )
  )
  left_join(dat, rowcol1536, by = c("qrt1536", "row384", "col384"), relationship = "many-to-many")
}

loadIrisFiles <- function(folder) {
  name <- list.files(folder, pattern = ".iris$")
  path <- list.files(folder, pattern = ".iris$", full.names = T)
  out <- lapply(path, fread)
  names(out) <- name
  rbindlist(out, idcol = "file")
}

plotReplicateCorrelation <- function(dat) {
  # remove factor levels with all NA's ------
  percentage_of_nas <- dat %>%
    group_by(cond) %>%
    select(contains("rep"), cond) %>%
    summarize(across(everything(), ~ mean(is.na(.))))

  percentage_of_nas[percentage_of_nas == 1] <- NA
  cols_to_remove <- percentage_of_nas %>%
    select(where(~ any(is.na(.)))) %>%
    select(contains("rep")) %>%
    colnames()

  dat <- select(dat, !all_of(cols_to_remove))
  # ------------------------------------

  # plot replicate correlation
  n_plots <- ncol(dat) - 4

  has_zero <- apply(dat[, 4:ncol(dat)] == 0, 1, sum, na.rm = T) != 0

  dat_o <- dat[!has_zero] %>%
    # filter(`rep1` > 0 & `rep2` > 0)  %>%

    mutate(across(where(is.numeric), ~ log10(.x + min_colony_size_opacity)))

  # This is very heavy-handed, but I cannot think of a better way to suppress the warnings due to (some) NA values here
  # TODO: Clean this up
  plot_o <- suppressWarnings(dat_o %>%
    {
      ggpairs(.,
        columns = 4:ncol(.), aes(color = cond, alpha = 0.3),
        lower = list(size = 0.5 / n_plots),
        upper = list(continuous = wrap("cor", size = 6 / sqrt(n_plots))),
        progress = FALSE
      ) +
        xlab(expression(log[10] * (opacity))) +
        ylab(expression(log[10] * (opacity))) +
        scale_color_brewer(palette = "Dark2") +
        theme_bw()
    })
  return(list(dat_o, plot_o))
}


plotHC <- function(dat, meta_cols = NULL, pivlong = TRUE, value_col_name = NULL) {
  # plot hierachical clustering
  # remove factor levels with all NA's ------
  percentage_of_nas <- dat %>%
    select(contains("rep"), cond) %>%
    group_by(cond) %>%
    summarize(across(everything(), ~ mean(is.na(.))))

  percentage_of_nas[percentage_of_nas == 1] <- NA
  cols_to_remove <- percentage_of_nas %>%
    select(where(~ any(is.na(.)))) %>%
    select(contains("rep")) %>%
    colnames()

  # TODO: Fix deprecationWarning here
  dat <- select(dat, !all_of(cols_to_remove))
  # ------------------------------------

  # make data wider still
  dat_hc <- dat %>%
    # super ugly but whtaever
    {
      if (pivlong) {
        (.) %>% pivot_longer(cols = (4):ncol(.), names_to = "rep")
      } else {
        (.)
      }
    } %>%
    pivot_wider(
      id_cols = genename,
      names_from = all_of(meta_cols),
      # names_from = c(cond, numb, biorep96, techrep96, plate_replicate),
      values_from = all_of(value_col_name)
    ) %>%
    setDT()

  m <- dat_hc[, -1] %>% as.matrix()

  rownames(m) <- dat_hc[["genename"]]
  p <- colnames(m)

  meta <- data.frame(foo = p) %>%
    separate_wider_delim(
      foo, "_",
      names = meta_cols
    ) %>%
    data.frame(., row.names = p)

  # mat = m %>% cor(use='complete.obs')
  #
  # o = rownames(mat)
  # hc = hclust(as.dist(1 - mat))  # as.dist: use the correlation distance
  # mat = mat[hc$order, hc$order]
  # mat[lower.tri(mat)] = NA
  # mat = mat[o, o]

  getBrewerColors <- function(var_from_meta, palette_name) {
    varname <- unique(meta[[var_from_meta]])
    colorRampPalette(RColorBrewer::brewer.pal(n = 8, palette_name))(length(varname)) %>%
      setNames(varname)
  }

  all_c_map_names <- c(
    "Dark2",
    "Accent",
    "Paired",
    "Pastel1",
    "Pastel2",
    "Set1",
    "Set2",
    "Set3"
  )
  annCols <- map2(
    meta_cols,
    1:length(meta_cols),
    \(x, y) {
      getBrewerColors(x, all_c_map_names[y])
    }
    # ~
  )
  names(annCols) <- meta_cols
  mm <- m + min_colony_size_opacity
  mm[is.na(mm)] <- 0
  pheatmap(
    # TODO
    log10(m + min_colony_size_opacity), # assume multiplicative errors
    annotation_col = select(meta, all_of(meta_cols)),
    annotation_colors = annCols,
    show_colnames = T,
    fontsize_row = 2,
    fontsize_col = 2,
    cluster_rows = hclust(dist(t(mm + min_colony_size_opacity))),
    cluster_cols = hclust(dist(mm + min_colony_size_opacity)),
    # clustering_method='single',
    na_col = "#FFFFFF",
    heatmap_legend_param = list(title = expression(log[10] * (opacity)))
  )
}


getResultsFromLinearModel <- function(
    dat,
    folder,
    type = "biorep_all",
    plot_pca_qc = T,
    normalize_how = NULL) {
  fitness <- dat %>%
    pivot_wider(
      id_cols = genename,
      names_from = c(cond, numb, rep),
      values_from = all_of(normalize_how)
    ) %>%
    setDT()
  m <- fitness[, -1] %>%
    as.matrix()
  genes <- fitness[["genename"]]
  rownames(m) <- genes
  f <- data.frame(genes)
  rownames(f) <- genes
  p <- colnames(m)

  meta <- data.frame(foo = p) %>%
    separate_wider_delim(
      foo, "_",
      names = c("cond", "conc", "rep")
    ) %>%
    data.frame(., row.names = p) %>%
    mutate(rep = gsub("ep", "", rep))

  eset <- ExpressionSet(m, AnnotatedDataFrame(meta), AnnotatedDataFrame(f))

  if (plot_pca_qc) {
    pdf(paste0("./output/", folder, "/qc_", type, "_pca.pdf"), h = 5, w = 10)
    par(mfrow = c(1, 2), pty = "s")
    plotMDS(eset, labels = pData(eset)[, "rep"], cex = 0.75)
    plotMDS(eset, labels = pData(eset)[, "cond"], cex = 0.75)
    dev.off()
  }

  design <- model.matrix(~ 0 + cond, data = pData(eset))
  colnames(design) <- gsub("cond", "", colnames(design))
  cm <- makeContrasts(ara = AraIPTG - IPTG, levels = design)
  fit <- lmFit(eset, method = "robust", design, maxit = 10000)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, genelist = fit2$genes, number = Inf, sort.by = "none") %>%
    janitor::clean_names()

  res %>%
    {
      ggplot(., aes(log_fc, -log10(adj_p_val))) +
        geom_point(pch = 21, col = "grey30", fill = "grey60") +
        ggrepel::geom_text_repel(
          data = filter(., adj_p_val < 0.01),
          aes(label = genes), max.overlaps = Inf
        ) +
        labs(x = expression(log[2] * "FC"), y = expression(-log[10] * "(adj p-val)")) +
        theme_cowplot()
    }

  ggsave(paste0("./output/", folder, "/res_volcano_", type, "__normalization_method_", normalize_how, ".pdf"), width = 7, height = 7)

  fwrite(res, paste0("./output/", folder, "/res_", type, "__normalization_method_", normalize_how, ".csv"))
  return(res)
}

get_opacity_measurement_plots <- function(
    iris,
    control_gene_name = control_gene_name,
    value_to_plot = NULL,
    output_name_base = NULL,
    mi,
    ma,
    gtitle) {
  gfp_controls_plot <- ggplot(
    data = iris %>%
      filter(genename == control_gene_name), aes(x = numb, y = .data[[value_to_plot]], color = cond)
  ) +
    # geom_hline(yintercept = 0, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 25000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 50000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 75000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 100000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 125000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 150000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    geom_boxplot() +
    geom_text(
      data = iris %>%
        filter(genename == control_gene_name) %>%
        group_by(
          cond, numb, folder_plate_id, colony_on_edge
        ) %>%
        summarize(
          {{ value_to_plot }} := median(.data[[value_to_plot]]),
          .groups = "drop_last"
        ), aes(label = .data[[value_to_plot]]), y = max(iris[[value_to_plot]]), size = 3
    ) +
    theme_presentation() +
    facet_grid(folder_plate_id ~ colony_on_edge) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0)
    ) +
    scale_color_manual(values = cond_color_vector) +
    ylim(c(mi, ma * 1.2))
  ggsave(
    plot = gfp_controls_plot + ggtitle(str_c(gtitle, " (only gfp controls)")),
    filename = str_c("./output/", output_name_base, "_gfp_controls.pdf"),
    width = 4 * sqrt(length(unique(iris$numb))),
    height = 18
  )

  all_values_plot <- ggplot(
    data = iris, aes(x = numb, y = .data[[value_to_plot]], color = cond)
  ) +
    # geom_hline(yintercept = 0, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 25000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 50000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 75000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 100000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 125000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    # geom_hline(yintercept = 150000, linetype = "longdash", colour = "grey", alpha = 0.35) +
    geom_boxplot() +
    geom_text(
      data = iris %>%
        group_by(
          cond, numb, folder_plate_id, colony_on_edge
        ) %>%
        summarize(
          {{ value_to_plot }} := median(.data[[value_to_plot]]),
          .groups = "drop_last"
        ), aes(label = .data[[value_to_plot]]), y = max(iris[[value_to_plot]]), size = 3
    ) +
    theme_presentation() +
    facet_grid(folder_plate_id ~ colony_on_edge) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0)
    ) +
    scale_color_manual(values = cond_color_vector) +
    ylim(c(mi, ma * 1.2))

  ggsave(
    plot = all_values_plot + ggtitle(str_c(gtitle)),
    filename = str_c("./output/", output_name_base, ".pdf"),
    width = 4 * sqrt(length(unique(iris$numb))),
    height = 18
  )
}


assign_random_color <- function(n, s) {
  library(RColorBrewer)

  # Set seed for reproducibility
  set.seed(s)

  # Get the number of unique values in n
  unique_n <- unique(n)
  num_colors <- length(unique_n)

  # Generate a set of colors from RColorBrewer's Set3 palette
  colors <- suppressWarnings(brewer.pal(12, "Set3"))
  colors <- sample(colors, num_colors)

  # Return the vector of colors
  return(colors)
}

get_iptg_scatter <- function(
    iris_hits_input,
    effect_size_names = c("effect_size_iptghigh", "effect_size_iptglow"),
    plot_means = NULL,
    nudge_x = 0) {
  if (plot_means) {
    effect_size_names <- str_c(effect_size_names, "_mean")
  }

  p <- iptg_high_low_plot <- ggplot() +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed", alpha = 0.3
    ) +
    geom_point(
      data = iris_hits_input,
      aes(x = .data[[effect_size_names[1]]], y = .data[[effect_size_names[2]]], color = hit, alpha = hit_alpha)
    ) +
    geom_text(
      data = iris_wide %>%
        group_by(folder) %>%
        summarize(co = round(cor(effect_size_iptghigh, effect_size_iptglow, use = "complete.obs"), 3)),
      aes(x = min(iris_wide$effect_size_iptghigh, na.rm = TRUE), y = max(iris_wide$effect_size_iptglow, na.rm = TRUE), label = co), size = 3, nudge_x = nudge_x
    ) +
    geom_text(
      data = iris_hits %>%
        group_by(folder, hit) %>%
        tally() %>%
        group_by(folder) %>%
        mutate(
          # y_offset = seq(1, 4, length.out = n())
          y_offset = case_when(
            hit == "both" ~ 1,
            hit == "high" ~ 2,
            hit == "low" ~ 3,
            hit == "none" ~ 4
          )
        ) %>%
        arrange(folder, hit),
      aes(x = min(iris_wide$effect_size_iptghigh, na.rm = TRUE), y = max(iris_wide$effect_size_iptglow, na.rm = TRUE) - y_offset, label = n, color = hit), nudge_x = nudge_x, size = 3
    ) +
    theme_presentation() +
    facet_wrap(. ~ folder) +
    scale_color_manual(values = hit_color_mapping) +
    scale_alpha_continuous(guide = "none", range = c(0.15, 0.5)) +
    xlab("Effect size (IPTG: High)") +
    ylab("Effect size (IPTG: Low)") +
    ggtitle("Effect sizes > 0: relatively large colony upon retron/toxin induction\nEffect sizes < 0: relatively small colony upon retron/toxin induction") +
    # Make title font smaller
    theme(plot.title = element_text(size = 12)) +
    xlim(c(min(iris_wide$effect_size_iptghigh, iris_wide$effect_size_iptglow, na.rm = TRUE), max(iris_wide$effect_size_iptghigh, iris_wide$effect_size_iptglow, na.rm = TRUE))) +
    ylim(c(min(iris_wide$effect_size_iptghigh, iris_wide$effect_size_iptglow, na.rm = TRUE), max(iris_wide$effect_size_iptghigh, iris_wide$effect_size_iptglow, na.rm = TRUE))) +
    NULL

  return(p)
}

get_condition_wise_replicate_correlation_matrix <- function(m, folder = NULL) {
  condition_wise_replicate_correlations <- m %>%
    group_by(cond, numb) %>%
    nest() %>%
    mutate(pairwise_cor_matrix = map(data, \(x)  {
      x %>%
        # select(-genename) %>%
        as.data.frame() %>%
        column_to_rownames("genename") %>%
        cor(., use = "complete.obs") %>%
        as.data.frame() %>%
        as.matrix()
    }))

  hm_list <- list()
  for (i in 1:dim(condition_wise_replicate_correlations)[1]) {
    x <- condition_wise_replicate_correlations$pairwise_cor_matrix[[i]]
    cond <- condition_wise_replicate_correlations$cond[[i]]
    numb <- condition_wise_replicate_correlations$numb[[i]]

    fo <- folder # whatever
    tmp <- data.frame(biorep_all = rownames(x)) %>%
      left_join(iris %>% select(cond, numb, folder, biorep96, techrep96, plate_replicate, biorep_all) %>% distinct() %>% filter(folder == fo) %>% arrange(biorep_all), by = "biorep_all")
    tmp <- tmp[tmp$cond == cond, ]
    tmp <- tmp[tmp$numb == numb, ]
    tmp <- tmp[match(tmp$biorep_all, rownames(x)), ]
    # stopifnot(all(tmp$biorep_all == rownames(x)))
    if (!all(tmp$biorep_all == rownames(x))) {
      print("Some bioreps are missing, entering debug mode")
      browser()
    }

    a <- assign_random_color(unique(tmp$biorep96), s = 1)
    names(a) <- unique(tmp$biorep96)
    b <- assign_random_color(unique(tmp$techrep96), s = 23)
    names(b) <- unique(tmp$techrep96)
    c <- assign_random_color(unique(tmp$plate_replicate), s = 3123)
    names(c) <- unique(tmp$plate_replicate)
    ha <- HeatmapAnnotation(
      biorep96 = as.character(tmp$biorep96),
      techrep96 = as.character(tmp$techrep96),
      plate_replicate = as.character(tmp$plate_replicate),
      col = list(
        "biorep96" = a,
        "techrep96" = b,
        "plate_replicate" = c
      )
    )

    xx <- Heatmap(
      x,
      name = "correlation",
      # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
      # colors should range from 0 to 1, white to red
      col = colorRamp2(c(min(x), 1), c("white", "red")),
      show_row_names = TRUE,
      show_column_names = TRUE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      # add title
      column_title = str_c(cond, numb),
      top_annotation = ha
    ) %>%
      draw() %>%
      grid.grabExpr()
    hm_list[[length(hm_list) + 1]] <- xx
  }
  return(hm_list)
}

rep_cor_qc_and_limma <- function() {
  res <- map(unique(iris$folder), function(folder) {
    out_folder <- paste0("./output/", folder)
    if (!dir.exists(out_folder)) dir.create(out_folder, recursive = T)

    fold <- iris %>%
      filter(folder == {{ folder }}) %>%
      as.data.table()
    reps <- c(
      # biol_replicate_column_name,
      # tech_replicate_column_name,
      # plate_replicate_column_name,
      "biorep_all"
    )

    # # split to wide  ---------------------
    dat_wide <- lapply(reps, function(x) {
      form <- as.formula(paste("... ~", x))
      # fold[[x]] <- paste0("rep", fold[[x]]) # this is achieved further up there
      fold %>%
        # take median values over other replicates
        ####################################################################################
        # !!! when x == biorep96 this will collapse the gfp replicates to a single value !!!
        ####################################################################################
        .[, .(opacity = median(opacity)), by = c("genename", "cond", "numb", x)] %>%
        dcast(form, value.var = "opacity")
    }) %>% setNames(reps)


    # plot-n-store the rep cor plots
    # plotReplicateCorrelation might remove some replicates due to NAs
    cor_info <- lapply(reps, function(x) {
      d <- plotReplicateCorrelation(dat_wide[[x]])
      data <- d[[1]]
      ploto <- d[[2]]
      # This is very heavy-handed, but I cannot think of a better way to suppress the warnings due to (some) NA values here
      # TODO: Clean this up
      suppressWarnings(ggsave(paste0(out_folder, "/qc_", x, "_correlation_median_opacity_over_other_replicates.pdf"),
        ploto,
        h = 12, w = 12
      ))
      return(data)
    })
    names(cor_info) <- reps
    # represnt biorep_all pairwise correlation matrix (one per condition) to judge the quality of the replicates
    # color by biorep96, techrep96, plate_replicate
    print("a)")
    hm_list <- get_condition_wise_replicate_correlation_matrix(cor_info[["biorep_all"]], folder = folder)
    pdf(file = paste0(out_folder, "/pairwise_replicate_correlation_matrix.pdf"), h = 9, w = 16)
    plot(wrap_plots(hm_list) + plot_layout(nrow = 2, byrow = TRUE, guides = "collect"))
    dev.off()

    print("b)")
    # plot-n-store the hierarchical clustering plots
    lapply(reps, function(x) {
      pdf(file = paste0(out_folder, "/qc_", x, "_clustering_median_opacity_over_other_replicates.pdf"), h = 8, w = 6)
      plt <- plotHC(dat_wide[[x]], meta_cols = c("cond", "numb", "rep"), value_col_name = "value")
      print("c")
      browser()
      draw(plt)
      dev.off()
    })

    print("c)")
    # Also without averaging, over all replicates, with sensible color scale
    pdf(file = paste0(out_folder, "/qc_", "clustering_median_opacity_over_other_replicates.pdf"), h = 8, w = 6)
    plt <- plotHC(
      fold %>%
        filter(genename != control_gene_name) %>%
        select(genename, cond, numb, opacity, contains("rep")),
      meta_cols = c("cond", "numb", colnames(fold)[str_detect(colnames(fold), "rep")]),
      pivlong = FALSE,
      value_col_name = "opacity"
    )
    draw(plt)
    dev.off()

    print("d")
    # TODO: Clean this up and compare z-score analysis to Limma
    dat <- dat_wide[["biorep_all"]] %>%
      filter(numb %in% c("10010201", "100100201", "10010500201", "1001001", "100105001")) %>%
      mutate(cond = gsub("Spectet|Spectetamp", "", cond))

    dat_long <- dat %>%
      # make sure to retain values_drop_na or the plotMDS will not work
      pivot_longer(contains("rep"), names_to = "rep", values_to = "opacity", values_drop_na = T) %>%
      mutate(rep = gsub("ep", "", rep)) %>%
      # This normalizes by 'biorep_all' (which are nested within the plates)
      # This comes from Vallo, originally
      group_by(cond, numb, rep) %>%
      mutate(
        fitness_median_gfp = opacity / median(opacity[genename == "gfp"]),
        fitness_z_score_by_plate = scale(opacity)[, 1]
      )

    # Remove the ones you have information only from one condition
    # Why? lmFit can report back logFC NA, but method='robust' will trip 'rlm' up
    only_zeros_in_a_cond <- dat_long %>%
      group_by(cond, genename) %>%
      summarize(mu = mean(fitness_median_gfp, na.rm = T), .groups = "drop_last") %>%
      filter(
        mu == 0
      ) %>%
      ungroup() %>%
      select(genename) %>%
      distinct() %>%
      pull(genename)

    # print(str_c("Removing ", length(only_zeros_in_a_cond), " genes with only zeros in a condition (not sure what this does...)"))
    dat_long <- dat_long %>% filter(!genename %in% only_zeros_in_a_cond)

    fitness_long <- dat_long %>%
      filter(genename != control_gene_name) %>%
      select(-opacity) %>%
      mutate(
        fitness_median_gfp_log2 = log2(fitness_median_gfp)
      )

    res_median_gfp_log2 <- getResultsFromLinearModel(fitness_long, folder, type = "biorep_all", normalize_how = "fitness_median_gfp_log2")
    res_z_score <- getResultsFromLinearModel(fitness_long, folder, type = "biorep_all", normalize_how = "fitness_z_score_by_plate")
    return(list(
      res_median_gfp_log2 = res_median_gfp_log2,
      res_z_score = res_z_score
    ))
  })
  return(res)
}

remove_via_pertial_string_matching <- function(dat, vector_of_substrings) {
  # dat is a tibble, vector_of_substrings is a list of strings
  # filter dat on the 'comment' column, removing all rows that contain any of the substrings
  # in vector_of_substrings
  return(
    dat %>%
      filter(!str_detect(comment, str_c(vector_of_substrings, collapse = "|")))
  )
}
