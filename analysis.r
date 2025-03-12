# To install ggembl: https://git.embl.de/grp-zeller/ggembl
pacman::p_load(data.table, tidyverse, ggrepel, GGally, scales, ComplexHeatmap, limma, Biobase, cowplot, janitor, yaml, RColorBrewer, vegan, patchwork, ggembl, circlize, RColorBrewer)
source("./helper.R")
load_yaml_to_global("./input/config.yaml")
# This gets loaded from the yaml file
names(cond_color_vector) <- cond_color_vector_levels

# Specify replicate information
# Experimental design can be complicated, with biological replicates and nested technical replicates.
# Below, I implement sanity checks that specify the expected relationship/number of replicates in the experiment.
# This will be used for sanity checks down the line

# biological replicates
biol_replicate_column_name <- "biorep96"
expected_num_biol_replicates <- 2

# technical replicates
tech_replicate_column_name <- "techrep96"
# expected_num_tech_replicates <- 2
# If a vector, it means number of technical replicates can vary
expected_num_tech_replicates <- c(2, 3)

# plate replicate
plate_replicate_column_name <- "plate_replicate"
# Some systems have 2, some 3 plate replicates
# TODO: This would be better solved by a mapping of systems to plate replicates, but ok for now
expected_num_plate_replicates <- c(2, 3)


# Specify experimental design

# This code block does 2 things
# - It loads the initial 96 well plate pipeting scheme, relating genes to their positions in the plate. It then maps individual wells of the 96 well plate to wells in the 384 well plate, and then to the 1536 well plate.
# - Furthermore, and importantly, the registerQuadrants function specify which plate are biol. replicates, and how these are later on further technically replicated in 384 well plates

# It is important to note that the final replication step (what I call 'plate replicate') is not encoded here, but instead encoded in the iris file names (see 'tail' in the config.yaml)


# Plate maps -------------------

# 1) create/adjust the 96w maps at './input/maps/library_plate_n.csv'

# 2) go from 96w to 384w to 1536
# TODO: Document the logic here somehow - maybe best in the git repo.
# first 384w plate
# second 384w plate
# third 384w plate
# fourth 384w plate
plate_number_vector <- c(
  1, 2, 3, 2,
  4, 5, 6, 6,
  1, 2, 3, 1,
  4, 5, 6, 5
)
bio_rep_vector <- c(
  1, 1, 1, 1,
  2, 2, 2, 2,
  1, 1, 1, 1,
  2, 2, 2, 2
)
tech_rep_vector <- c(
  1, 1, 1, 2,
  1, 1, 1, 2,
  2, 3, 2, 3,
  2, 2, 3, 3
)

rep_information <- data.frame(
  plate_number = plate_number_vector,
  bio_rep = bio_rep_vector,
  tech_rep = tech_rep_vector
) %>%
  as_tibble()

print(str_c("We array a total of ", length(plate_number_vector), " 96 well plates with the following counts"))
tmp <- table(plate_number_vector)
names(tmp) <- str_c("Plate ", names(tmp))
print(tmp)

# Some diagnostics
stopifnot(length(bio_rep_vector) == length(tech_rep_vector))
if (length(unique(table(bio_rep_vector))) != 1) {
  stop("Number of biological replicates is not consistent across 96 well plates, check how you do define bio_rep_vector")
} else {
  print(str_c("Number of biological replicates (inferred from bio_rep_vector): ", length(unique(bio_rep_vector))))
  num_biol_repl <- length(unique(bio_rep_vector))
}

# TODO: I think what this table is printing is at the least misleading, if not completely wrong.
print("Overall, I'm seeing this:")
tmp <- rep_information %>%
  group_by(plate_number, bio_rep) %>%
  tally() %>%
  rename(tech_reps = n) %>%
  print()
w <- FALSE
if (length(unique(tmp$tech_reps)) != 1) {
  w <- TRUE
  warning("Number of technical replicates is not consistent over plates; this is possible but please make sure this is intended.")
}
if (w) {
  # seems idiotic but try putting it in the condition :)
  Sys.sleep(2)
}



map96to384quadrants <- registerQuadrants("384w",
  plate_number = plate_number_vector,
  bio_rep = bio_rep_vector,
  tech_rep = tech_rep_vector
)

map384to1536quadrants <- registerQuadrants("1536w",
  # TODO: TO be consistent with google scheme, make this work with letters
  plate_number = c(1, 2, 3, 4)
) # TODO: For safety reasons, maybe remove biorep384, techrep384???


folders <- list.files(iris_input_folder, full.names = T)

iris <- lapply(folders, loadIrisFiles) %>%
  setNames(gsub(paste0(iris_input_folder, "/"), "", folders)) %>%
  rbindlist(idcol = "folder") %>%
  group_by(folder) %>% # to annotate the same way in each folder
  mutate(plt1536 = match(file, unique(file))) %>%
  separate(file, column_names, sep = "-") %>%
  mutate(tail = gsub(".JPG.iris", "", tail)) %>%
  separate(tail, c("system_num", "plate_replicate")) %>% # This line generates the 'plate replicate'
  rename(row1536 = row, col1536 = column) %>%
  left_join(fread("./input/maps/systems.csv", colClasses = "character"), by = "system_num")

iris <- left_join(map96to384quadrants, map384to1536quadrants, relationship = "many-to-many", by = "plt384") %>%
  add384RowAndCol() %>%
  add1536RowAndCol() %>%
  right_join(iris, by = c("row1536", "col1536")) %>%
  as_tibble() %>%
  left_join(read96wMaps(), by = c("plt96", "biorep96", "row96", "col96")) %>%
  mutate(
    colony_id  = interaction(plt1536, row1536, col1536),
    biorep_all = interaction(across(all_of(c(biol_replicate_column_name, tech_replicate_column_name, plate_replicate_column_name)))) %>% as.numeric() %>% str_pad(2, pad = "0")
  ) %>%
  mutate(across(contains("rep"), \(x) paste0("rep", x))) %>%
  mutate(plate_id = interaction(folder, cond, plate_replicate, numb)) %>%
  setDT() %>%
  as_tibble()

# Clean up a bit
iris <- iris %>%
  select(
    cond, numb, genename, system_desc, opacity, colony_id, folder, row1536, col1536, plate_id,
    all_of(c(biol_replicate_column_name, tech_replicate_column_name, plate_replicate_column_name, "biorep_all"))
  ) %>%
  # mutate(across(contains("rep"), as.numeric)) %>%
  # TODO: Rename biol_repliate_column_name to bio_rep and tech_replicate_column_name to tech_rep
  as_tibble()


# Sanity check: iris object must not contain any NAs at this step
if (any(is.na(iris))) {
  stop("Collated results (iris object) contain NAs. This probably means something could not be joined in quadrant table joining. Please check the input files and ensure map files/configuration is correct.")
}

# Add edge information
mi_row_1536 <- min(iris$row1536)
ma_row_1536 <- max(iris$row1536)
mi_col_1536 <- min(iris$col1536)
ma_col_1536 <- max(iris$col1536)
iris <- iris %>%
  # TODO: It would probably be a good idea to leave gfp in here before z-scoring, but right now this isn't easy the way the data is formatted
  # filter(genename == control_gene_name) %>%
  # The edge definition should be in line with what Jacob describes in the tac/tic paper
  mutate(
    colony_on_edge =
      row1536 == mi_row_1536 |
        row1536 == (mi_row_1536 + 1) |
        row1536 == (mi_row_1536 + 2) |
        row1536 == (mi_row_1536 + 3) |
        row1536 == (mi_row_1536 + 4) |
        row1536 == ma_row_1536 |
        row1536 == (ma_row_1536 - 1) |
        row1536 == (ma_row_1536 - 2) |
        row1536 == (ma_row_1536 - 3) |
        row1536 == (ma_row_1536 - 4) |
        col1536 == mi_col_1536 |
        col1536 == (mi_col_1536 + 1) |
        col1536 == (mi_col_1536 + 2) |
        col1536 == (mi_col_1536 + 3) |
        col1536 == (mi_col_1536 + 4) |
        col1536 == ma_col_1536 |
        col1536 == (ma_col_1536 - 1) |
        col1536 == (ma_col_1536 - 2) |
        col1536 == (ma_col_1536 - 3) |
        col1536 == (ma_col_1536 - 4)
  ) %>%
  # mutate(plate_id = interaction(cond, plate_replicate, numb)) %>%
  mutate(folder_plate_id = interaction(folder, plate_replicate)) %>%
  mutate(folder_plate_id = factor(folder_plate_id, levels = sort(unique(as.character(folder_plate_id))))) %>%
  mutate(cond = factor(cond,
    levels =
      c(
        unique(cond)[str_detect(unique(cond), "Ara")],
        unique(cond)[!str_detect(unique(cond), "Ara")]
      )
  ))

# Save original iris object for later
iris_orig <- iris

# TODO: Move this into function
walk(unique(iris$folder), function(folder) {
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

  hm_list <- get_condition_wise_replicate_correlation_matrix(cor_info[["biorep_all"]], folder = folder)
  pdf(file = paste0(out_folder, "/pairwise_replicate_correlation_matrix.pdf"), h = 9, w = 16)
  plot(wrap_plots(hm_list) + plot_layout(nrow = 2, byrow = TRUE, guides = "collect"))
  dev.off()

  # plot-n-store the hierarchical clustering plots
  lapply(reps, function(x) {
    pdf(file = paste0(out_folder, "/qc_", x, "_clustering_median_opacity_over_other_replicates.pdf"), h = 8, w = 6)
    plt <- plotHC(dat_wide[[x]], meta_cols = c("cond", "numb", "rep"), value_col_name = "value")
    draw(plt)
    dev.off()
  })

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

  # TODO: Clean this up and compare z-score analysis to Limma
  return() # for now
  dat <- dat_wide[["biorep_all"]] %>%
    # TODO: This looks hard-coded and dangerous
    # TODO: We need to fix the fucking format of this, it is impossible for me to make sense of
    filter(numb %in% c("10010201", "10010500201", "1001001", "100105001")) %>%
    mutate(cond = gsub("Spectet|Spectetamp", "", cond))

  dat_long <- dat %>%
    # make sure to retain values_drop_na or the plotMDS will not work
    pivot_longer(contains("rep"), names_to = "rep", values_to = "opacity", values_drop_na = T) %>%
    mutate(rep = gsub("ep", "", rep)) %>%
    group_by(cond, numb, rep) %>%
    mutate(
      fitness_median_gfp = opacity / median(opacity[genename == "gfp"]),
      fitness_z_score_by_plate = scale(opacity)[, 1]
    )

  # Remove the ones you have information only from one condition
  # Why? lmFit can report back logFC NA, but method='robust' will trip 'rlm' up
  only_zeros_in_a_cond <- dat_long %>%
    group_by(cond, genename) %>%
    summarize(mu = mean(fitness_median_gfp, na.rm = T)) %>%
    filter(
      mu == 0
    ) %>%
    ungroup() %>%
    select(genename) %>%
    distinct() %>%
    pull(genename)

  print(str_c("Removing ", length(only_zeros_in_a_cond), " genes with only zeros in a condition (not sure what this does...)"))
  dat_long <- dat_long %>% filter(!genename %in% only_zeros_in_a_cond)

  fitness_long <- dat_long %>%
    filter(genename != control_gene_name) %>%
    select(-opacity) %>%
    mutate(
      fitness_median_gfp_log2 = log2(fitness_median_gfp)
    )

  res_z_score <- getResultsFromLinearModel(fitness_long, folder, type = "biorep_all", normalize_how = "z_score_by_plate")
  res_median_gfp <- getResultsFromLinearModel(fitness_long, folder, type = "biorep_all", normalize_how = "median_gfp_log2")
})

##################################
# Z-score-based analysis
##################################

# Reload original iris object
iris <- iris_orig %>%
  # Remove the doubleplasmid experiment for now
  # TODO: We need to change the encoding for the doubleplasmid experiments and then we can include this
  filter(!str_detect(folder, "doubleplasmid")) %>%
  {
    if (visualize_only_control_and_treatment_conditions) {
      (.) %>% filter(grepl("10010201", numb) | grepl("1001021", numb) | grepl("1001001", numb) | grepl("100101", numb)) # SpectetAraIPTG (library + system induced, with low and high concentrations) and SpectetIPTG (only system induced))
    } else {
      (.)
    }
  }

# Plot raw opacities by system/plate and edge
value_tp <- "opacity"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])

get_opacity_measurement_plots(
  iris,
  "gfp",
  value_tp,
  "opacity_raw",
  mi,
  ma,
  "Raw opacity values"
)

# Normalize opacity values edge effects wrt the center of the plates (ALL values, not just gfp controls)
edge_correction_factor_per_plate <- iris %>%
  group_by(
    cond, numb, folder, plate_replicate, colony_on_edge
  ) %>%
  summarize(opacity = median(opacity), .groups = "drop_last") %>%
  pivot_wider(
    id_cols = c(cond, numb, folder, plate_replicate),
    names_from = colony_on_edge,
    values_from = opacity
  ) %>%
  mutate(edge_correction_factor = `TRUE` / `FALSE`) %>%
  select(folder, cond, numb, plate_replicate, edge_correction_factor)
iris <- iris %>%
  left_join(edge_correction_factor_per_plate, by = c("cond", "numb", "folder", "plate_replicate")) %>%
  mutate(
    opacity_edge_corrected = ifelse(
      colony_on_edge,
      as.integer(opacity / edge_correction_factor),
      opacity
    )
  )
# Visualize edge-corrected opacity values
value_tp <- "opacity_edge_corrected"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])
get_opacity_measurement_plots(iris, "gfp", value_tp, "opacity_edge_corrected", mi, ma, "Edge-corrected opacity values")

iris <- iris %>%
  filter(grepl("10010201", numb) | grepl("1001021", numb) | grepl("1001001", numb) | grepl("100101", numb)) # SpectetAraIPTG (library + system induced, with low and high concentrations) and SpectetIPTG (only system induced)

iris <- iris %>%
  select(-opacity) %>% # We have edge-controlled opacity now, so we dont need this anymore
  group_by(
    folder, cond, numb, folder_plate_id
  ) %>%
  mutate(
    # Get z-scores of edge-corrected opacity values by *plate* (this should account for plate effects)
    ec_opacity_z_scored_by_plate_including_gfp = round(scale(opacity_edge_corrected)[, 1], 3),
    # Also get ec_opacity corrected by a plates (median) gfp opacity - I think this should give us very similar results in the end
    ec_opacity_corrected_by_plate_using_gfp = round(opacity_edge_corrected / median(opacity_edge_corrected[genename == control_gene_name]), 3)
  ) %>%
  group_by(
    folder, cond, numb, folder_plate_id, biorep96
  ) %>%
  mutate(
    # Get z-scores of edge-corrected opacity values by *plate* (this should account for plate effects)
    ec_opacity_z_scored_by_plate_and_biorep_including_gfp = round(scale(opacity_edge_corrected)[, 1], 3),
    # Also get ec_opacity corrected by a plates (median) gfp opacity - I think this should give us very similar results in the end
    ec_opacity_corrected_by_plate_and_biorep_using_gfp = round(opacity_edge_corrected / median(opacity_edge_corrected[genename == control_gene_name]), 3)
  )

# Visualize the edge-corrected and normalized values
value_tp <- "ec_opacity_z_scored_by_plate_including_gfp"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])
get_opacity_measurement_plots(iris, "gfp", value_tp, "ec_opacity_z_scored_by_plate_including_gfp", mi, ma, "Z-scored, edge-corrected opacity values (by plate)")
value_tp <- "ec_opacity_z_scored_by_plate_and_biorep_including_gfp"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])
get_opacity_measurement_plots(iris, "gfp", value_tp, "ec_opacity_z_scored_by_plate_and_biorep_including_gfp", mi, ma, "Z-scored, edge-corrected opacity values (by plate and biorep)")
value_tp <- "ec_opacity_corrected_by_plate_using_gfp"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])
get_opacity_measurement_plots(iris, "gfp", value_tp, "ec_opacity_corrected_by_plate_using_gfp", mi, ma, "edge-corrected opacity values\nnormalized by plate-wise median gfp")
value_tp <- "ec_opacity_corrected_by_plate_and_biorep_using_gfp"
mi <- min(iris[[value_tp]])
ma <- max(iris[[value_tp]])
get_opacity_measurement_plots(iris, "gfp", value_tp, "ec_opacity_corrected_by_plate_using_gfp", mi, ma, "edge-corrected opacity values\nnormalized by plate- and biorep-wise median gfp")

p <- iris %>%
  ungroup() %>%
  select(
    folder,
    ec_opacity_z_scored_by_plate_including_gfp,
    # ec_opacity_corrected_by_plate_using_gfp,
    ec_opacity_z_scored_by_plate_and_biorep_including_gfp
    # ec_opacity_corrected_by_plate_and_biorep_using_gfp
  ) %>%
  rename(
    `ec_opacity_z_scored_by\nplate` = ec_opacity_z_scored_by_plate_including_gfp,
    # `ec_opacity_corrected_by\nplate_using_gfp` = ec_opacity_corrected_by_plate_using_gfp,
    `ec_opacity_z_scored_by\nplate_and_biorep` = ec_opacity_z_scored_by_plate_and_biorep_including_gfp,
    # `ec_opacity_corrected_by\nplate_and_biorep_using_gfp` = ec_opacity_corrected_by_plate_and_biorep_using_gfp
  ) %>%
  ggplot() +
  geom_point(aes(
    x = `ec_opacity_z_scored_by\nplate`,
    y = `ec_opacity_z_scored_by\nplate_and_biorep`
  ), alpha = 0.05) +
  facet_wrap(
    . ~ folder,
    ncol = 2
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    alpha = 0.3
  ) +
  theme_presentation() +
  xlab("Normalized z-score differences (i.e. effect size)\nby PLATE") +
  ylab("Normalized z-score differences (i.e. effect size)\nby PLATE and BIOREP")

ggsave(
  plot = p,
  filename = str_c(out_folder, "/", "normalized_effect_sizes_comparison.pdf"),
  width = 5.5,
  height = 6.5,
  dpi = 300
)

iris <- iris %>%
  # ! remove the control gene from the analysis
  filter(genename != control_gene_name) %>%
  mutate(iptgconc = case_when(
    numb == "1001001" ~ "iptg_low",
    numb == "100101" ~ "iptg_high",
    numb == "10010201" ~ "iptg_low",
    numb == "1001021" ~ "iptg_high",
    TRUE ~ NA
  ))

if (any(is.na(iris$iptgconc))) {
  stop("iris contains NAs in iptgconc. This should never happen. This probably means something in your encoding changed. Please talk to Nic.")
}

iris <- iris %>%
  group_by(
    folder, iptgconc, genename, biorep_all
  ) %>%
  # Filter 1: Remove strains where opacity is less than min_colony_size_opacity in IPTG plates (both low- and high concentration)
  filter(opacity_edge_corrected[cond == "SpectetIPTG"] > min_colony_size_opacity) %>%
  # TODO: Filter 2: Remove mucoid colonies - I'm actually skipping this for now since I don't quite understand what the measure is
  # (from Bobonis et al: "[...] (2) mucoid in the control plates44 (colony densities of both replicates > 51 [..]").
  # TODO: Filter 3: noisy strains in control/treatment plates standard deviation > 23000 (this is taken 1-1 from Bobonis et al.)
  # I'm also skipping this for now
  ungroup()

iris <- iris %>%
  mutate(
    cond = case_when(
      numb == "10010201" ~ "SpectetAraIPTGLow",
      numb == "1001021" ~ "SpectetAraIPTGHigh",
      numb == "1001001" ~ "SpectetIPTGLow",
      numb == "100101" ~ "SpectetIPTGHigh"
    )
  )

# iris_wide <-
iris_wide <- iris %>%
  # TODO: This is hacky for now and needs to be cleaned up, but probably on the raw data sid
  pivot_wider(
    id_cols = c(folder, genename, biorep_all),
    names_from = cond,
    values_from = all_of(effect_size_name_for_non_limma_analyis) # any of: 'ec_opacity_z_scored_by_plate_including_gfp', 'ec_opacity_corrected_by_plate_using_gfp', ""
  ) %>%
  mutate(
    effect_size_iptghigh = SpectetAraIPTGHigh - SpectetIPTGHigh,
    effect_size_iptglow = SpectetAraIPTGLow - SpectetIPTGLow
  ) %>%
  arrange(desc(effect_size_iptghigh))

iris_hits <- iris_wide %>%
  group_by(folder, genename) %>%
  summarize(
    effect_size_iptghigh_mean = mean(effect_size_iptghigh, na.rm = T),
    effect_size_iptghigh_sd = sd(effect_size_iptghigh, na.rm = T),
    lower_iptghigh = effect_size_iptghigh_mean - effect_size_iptghigh_sd,
    upper_iptghigh = effect_size_iptghigh_mean + effect_size_iptghigh_sd,
    # Same but for iptghigh
    effect_size_iptglow_mean = mean(effect_size_iptglow, na.rm = T),
    effect_size_iptglow_sd = sd(effect_size_iptglow, na.rm = T),
    lower_iptglow = effect_size_iptglow_mean - effect_size_iptglow_sd,
    upper_iptglow = effect_size_iptglow_mean + effect_size_iptglow_sd,
    .groups = "drop_last"
  )

if (any(is.na(iris_hits))) {
  stop("iris_hits contains NAs. I was afraid this might be the case, but I haven't implemented a solution yet. Talk to Nic please.")
}

iris_hits <- iris_hits %>%
  # call hits and add color helpers for plotting
  mutate(
    is_hit_iptg_high = ifelse(!str_detect(folder, "Retron"), effect_size_iptghigh_mean > z_score_cutoff, effect_size_iptghigh_mean < -z_score_cutoff),
    is_hit_iptg_low = ifelse(!str_detect(folder, "Retron"), effect_size_iptglow_mean > z_score_cutoff, effect_size_iptglow_mean < -z_score_cutoff),
    ## for blockers we use blue "#377eb8ff", for triggers green
    # color = ifelse(is_hit_iptg_high, "green4", "grey65")
    hit = case_when(
      is_hit_iptg_high & is_hit_iptg_low ~ "both",
      is_hit_iptg_high ~ "high",
      is_hit_iptg_low ~ "low",
      TRUE ~ "none"
    ),
    hit_alpha = case_when(
      is_hit_iptg_high & is_hit_iptg_low ~ 1,
      is_hit_iptg_high ~ 1,
      is_hit_iptg_low ~ 1,
      TRUE ~ 0.15
    ),
    color = case_when(
      is_hit_iptg_high & is_hit_iptg_low ~ "green4",
      is_hit_iptg_high ~ "#9c5703",
      is_hit_iptg_low ~ "#e6b967",
      TRUE ~ "grey65"
    )
  ) %>%
  mutate(hit = factor(hit, levels = c("both", "high", "low", "none")))

# get the color mapping like so, bit ugly but what can you do
tmp <- iris_hits %>%
  ungroup() %>%
  select(hit, color) %>%
  distinct()
hit_color_mapping <- tmp$color
names(hit_color_mapping) <- tmp$hit
hit_alpha_mapping <- map2_dbl(hit_color_mapping, names(hit_color_mapping), \(entry, na)
if (na == "none") {
  return(0.25)
} else {
  return(0.75)
})


# Scatterplots of effect sizes between high- and low induction.
# Make sure to include the effect size name in the output file
iptg_high_low_plot_one_point_per_gene <- get_iptg_scatter(
  iris_hits_input = iris_hits %>%
    arrange(desc(hit)),
  plot_means = TRUE
)

ggsave(
  plot = iptg_high_low_plot_one_point_per_gene,
  filename = str_c(out_folder, "/iptg_high_low_scatter_one_point_per_gene_", effect_size_name_for_non_limma_analyis, ".pdf"),
  width = 7,
  height = 4.875
)

## I'm disabling this for now since it's throwing warnings about missing values
## These come from the low colony size filter above
## This plot is anyway nto that useful, I think.
# iptg_high_low_plot_all_measurements <- get_iptg_scatter(
#   iris_hits_input = iris_wide %>%
#     left_join(iris_hits, by = c("folder", "genename")) %>%
#     arrange(desc(hit)),
#   plot_means = FALSE
# )

# ggsave(
#   plot = iptg_high_low_plot_all_measurements,
#   filename = str_c(out_folder, "/iptg_high_low_scatter_all_measurements_", effect_size_name_for_non_limma_analyis, ".pdf"),
#   width = 8,
#   height = 4.5
# )

z_score_plot_high <- ggplot(data = iris_hits) +
  geom_pointrange(
    aes(x = genename, y = effect_size_iptghigh_mean, ymin = effect_size_iptghigh_mean - effect_size_iptghigh_sd, ymax = effect_size_iptghigh_mean + effect_size_iptghigh_sd),
    pch = 21, size = 0.5
  ) +
  geom_point(aes(x = genename, y = effect_size_iptghigh_mean), shape = 21, size = 0.5, alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_identity() +
  labs(
    x = "Genes",
    y = "Triggering/Blocking score"
  ) +
  geom_text_repel(data = iris_hits %>% filter(is_hit_iptg_high), aes(x = genename, y = effect_size_iptghigh_mean, label = genename), max.overlaps = Inf) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.2, "cm")
  ) +
  facet_wrap(. ~ folder, ncol = 1) +
  geom_hline(data = data.frame(folder = unique(iris_hits$folder)) %>%
    mutate(z_score_cutoff = ifelse(!str_detect(folder, "Retron"), z_score_cutoff, -z_score_cutoff)), aes(yintercept = z_score_cutoff), linetype = "longdash", colour = "grey", linewidth = 0.5)

z_score_plot_low <- ggplot(data = iris_hits) +
  geom_pointrange(
    aes(x = genename, y = effect_size_iptglow_mean, ymin = effect_size_iptglow_mean - effect_size_iptglow_sd, ymax = effect_size_iptglow_mean + effect_size_iptglow_sd),
    pch = 21, size = 0.5
  ) +
  geom_point(aes(x = genename, y = effect_size_iptglow_mean), shape = 21, size = 0.5, alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_identity() +
  labs(
    x = "Genes",
    y = "Triggering/Blocking score"
  ) +
  geom_text_repel(data = iris_hits %>% filter(is_hit_iptg_low), aes(x = genename, y = effect_size_iptglow_mean, label = genename), max.overlaps = Inf) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.2, "cm")
  ) +
  facet_wrap(. ~ folder, ncol = 1) +
  geom_hline(data = data.frame(folder = unique(iris_hits$folder)) %>%
    mutate(z_score_cutoff = ifelse(!str_detect(folder, "Retron"), z_score_cutoff, -z_score_cutoff)), aes(yintercept = z_score_cutoff), linetype = "longdash", colour = "grey", linewidth = 0.5)


ggsave(
  plot = (z_score_plot_high + ggtitle("IPTG: High") + theme(plot.title = element_text(size = 20, face = "bold"))) + (z_score_plot_low + ggtitle("IPTG: Low") + theme(plot.title = element_text(size = 20, face = "bold"))) + plot_layout(ncol = 2),
  filename = str_c(out_folder, "/z_score_plot.pdf"),
  width = 16,
  height = 20
)

#########################
# Produce gene-wise plots
########################
mi <- min(iris[[effect_size_name_for_non_limma_analyis]])
ma <- max(iris[[effect_size_name_for_non_limma_analyis]])
data <- iris %>%
  as_tibble() %>%
  group_by(folder, genename) %>%
  mutate(cond = factor(cond, levels = c(
    unique(cond)[str_detect(unique(cond), "Ara")],
    unique(cond)[!str_detect(unique(cond), "Ara")]
  ))) %>%
  mutate(cond_ara = ifelse(str_detect(cond, "Ara"), "ara_yes", "ara_no")) %>%
  mutate(biorep96 = factor(biorep96)) %>%
  mutate(techrep96 = factor(techrep96)) %>%
  mutate(plate_replicate = factor(plate_replicate)) %>%
  mutate(iptgconc_plate_id = interaction(iptgconc, plate_replicate, sep = ".plate")) %>%
  mutate(iptgconc_plate_id = factor(iptgconc_plate_id, levels = sort(unique(as.character(iptgconc_plate_id))))) %>%
  arrange(plate_id)

data_plot <- data %>%
  group_by(folder, genename) %>%
  mutate(plate_replicate = str_c("Plate: ", plate_replicate)) %>%
  mutate(biorep96 = str_c("biorep: ", biorep96)) %>%
  nest() %>%
  mutate(
    plottt = map(data, \(x) {
      return(
        ggplot(x, aes(x = plate_replicate, y = .data[[effect_size_name_for_non_limma_analyis]], color = cond_ara, group = cond_ara, shape = techrep96)) +
          # geom_boxplot(outlier.color = NA) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.125, jitter.height = 0), alpha = 0.5) +
          theme_presentation() +
          facet_wrap(iptgconc ~ biorep96, ncol = 2, scales = "free_x") +
          labs(title = genename) +
          ylim(c(-5, ma * 1.01)) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
          ylab("Effect size")
      )
    })
  )

pwalk(list(data_plot$plottt, data_plot$folder, data_plot$genename), \(plo, fol, gn) {
  # if str_c(out_folder, "/", fol, '/opacities_per_gene') doesnt exit, generate it
  if (!dir.exists(str_c(out_folder, "/", fol, "/opacities_per_gene"))) {
    dir.create(str_c(out_folder, "/", fol, "/opacities_per_gene"), recursive = TRUE)
  }
  ggsave(
    plot = plo,
    filename = str_c(out_folder, "/", fol, "/opacities_per_gene", "/", gn, "_", effect_size_name_for_non_limma_analyis, ".pdf"),
    width = 6,
    height = 5
  )
})
