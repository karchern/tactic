# Load z score results
zscoreresults <- read_tsv("output/z_score_effect_sizes.tsv") %>%
    rename(
        system = folder,
        gene = genename
    ) %>%
    select(system, gene, effect_size_iptghigh_mean, effect_size_iptglow_mean)
# load limma results
limmaresults <- list.files("output", full.names = TRUE, recursive = TRUE) %>%
    grep("res_", ., value = TRUE) %>%
    grep("normalization_method", ., value = TRUE) %>%
    grep(".csv", ., value = TRUE)
limmaresults <- map(limmaresults, \(x) {
    return(
        read_csv(x) %>%
            mutate(tmp = x) %>%
            mutate(system = str_split_fixed(tmp, "/", 3)[, 2]) %>%
            mutate(fitted_value = str_split_fixed(tmp, "/", 3)[, 3]) %>%
            mutate(fitted_value = str_replace(fitted_value, "res_biorep_all__normalization_method_fitness_", "")) %>%
            mutate(fitted_value = str_replace(fitted_value, ".csv", "")) %>%
            select(-tmp)
    )
}) %>%
    bind_rows() %>%
    rename(
        gene = genes,
        system = system
    ) %>%
    select(system, gene, log_fc)
# Compare the two
p <- inner_join(
    zscoreresults,
    limmaresults
) %>%
    ggplot() +
    geom_point(
        aes(
            x = effect_size_iptghigh_mean,
            y = log_fc
        ),
        alpha = 0.2
    ) +
    xlab("Z score effect size (Z score)") +
    ylab("Limma effect size (log fc)") +
    geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        alpha = 0.35
    ) +
    facet_wrap(~system, nrow = 2) +
    theme_presentation()

ggsave(
    plot = p,
    filename = "tmp_compare_zscore_to_limma.pdf",
    width = 8,
    height = 6
)
