library(tidyverse)
library(ggembl)

prep <- function(x) {
    return(
        x %>%
            select(folder, genename, effect_size_iptghigh_mean, effect_size_iptghigh_sd, effect_size_iptglow_mean, effect_size_iptglow_sd) %>%
            mutate(
                genename = str_replace(genename, "T5_GG_", "T5_")
            )
    )
}

screen_1_with_outliers <- read_tsv("/Users/karcher/tactic/tagged_archives_unpacked/tactic-screen_1_WITH_wrong_clones/output/screen_1/z_score_effect_sizes.tsv")
screen_1_without_outliers <- read_tsv("/Users/karcher/tactic/tagged_archives_unpacked/tactic-screen_1_WITHOUT_wrong_clones/output/screen_1/z_score_effect_sizes.tsv")
screen_2 <- read_tsv("/Users/karcher/tactic/tagged_archives_unpacked/tactic-screen_2/output/screen_2/z_score_effect_sizes.tsv")

screen_1_with_outliers <- screen_1_with_outliers %>%
    prep()
screen_1_without_outliers <- screen_1_without_outliers %>%
    prep()
screen_2 <- screen_2 %>%
    prep()

data <- inner_join(
    screen_1_with_outliers,
    screen_1_without_outliers,
    by = c("folder", "genename"),
    suffix = c("_screen_1_with_outliers", "_screen_1_without_outliers")
) %>%
    inner_join(
        screen_2 %>%
            rename(
                effect_size_iptghigh_mean_screen_2 = effect_size_iptghigh_mean,
                effect_size_iptghigh_sd_screen_2 = effect_size_iptghigh_sd,
                effect_size_iptglow_mean_screen_2 = effect_size_iptglow_mean,
                effect_size_iptglow_sd_screen_2 = effect_size_iptglow_sd
            ),
        by = c("folder", "genename"),
    )

genes_to_label <- screen_2 %>%
    mutate(
        show_labels = abs(effect_size_iptghigh_mean) > 3
    ) %>%
    select(folder, genename, show_labels)

data <- left_join(data, genes_to_label)



# screen_1: with outliers vs. without outliers
p1 <- data %>%
    ggplot(aes(x = effect_size_iptghigh_mean_screen_1_with_outliers, y = effect_size_iptghigh_mean_screen_1_without_outliers)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    geom_text_repel(
        data = data %>% filter(show_labels), aes(label = genename), max.overlaps = Inf
    ) +
    labs(x = "mean effect size (screen 1 with outliers)", y = "mean effect size (screen 1 without outliers)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

p1.sd <- data %>%
    ggplot(aes(x = effect_size_iptghigh_sd_screen_1_with_outliers, y = effect_size_iptghigh_sd_screen_1_without_outliers)) +
    geom_text_repel(
        data = data %>% filter(show_labels), aes(label = genename), max.overlaps = Inf
    ) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "SD effect size (screen 1 with outliers)", y = "SD effect size (screen 1 without outliers)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

# p1
# p1.sd

# screen_1: with outliers vs. screen_2
p2 <- data %>%
    ggplot(aes(x = effect_size_iptghigh_mean_screen_1_with_outliers, y = effect_size_iptghigh_mean_screen_2)) +
    geom_point(alpha = 0.3) +
    geom_text_repel(
        data = data %>% filter(show_labels), aes(label = genename), max.overlaps = Inf
    ) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "mean effect size (screen 1 with outliers)", y = "mean effect size (screen 2)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

p2.sd <- data %>%
    ggplot(aes(x = effect_size_iptghigh_sd_screen_1_with_outliers, y = effect_size_iptghigh_sd_screen_2)) +
    geom_point(alpha = 0.3) +
    geom_text_repel(
        data = data %>% filter(show_labels), aes(label = genename), max.overlaps = Inf
    ) +
    labs(x = "mean effect size (screen 1 with outliers)", y = "mean effect size (screen 1 without outliers)") +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "SD effect size (screen 1 with outliers)", y = "SD effect size (screen 2)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

# p2
# p2.sd

# screen_1: without outliers vs. screen_2
p3 <- data %>%
    ggplot(aes(x = effect_size_iptghigh_mean_screen_1_without_outliers, y = effect_size_iptghigh_mean_screen_2)) +
    geom_point(alpha = 0.3) +
    geom_text_repel(
        data = data %>% filter(show_labels), aes(label = genename), max.overlaps = Inf
    ) +
    labs(x = "mean effect size (screen 1 with outliers)", y = "mean effect size (screen 1 without outliers)") +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "mean effect size (screen 1 without outliers)", y = "mean effect size (screen 2)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

p3.sd <- data %>%
    ggplot(aes(x = effect_size_iptghigh_sd_screen_1_without_outliers, y = effect_size_iptghigh_sd_screen_2)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "SD effect size (screen 1 without outliers)", y = "SD effect size (screen 2)") +
    facet_wrap(. ~ folder, ncol = 2) +
    theme_presentation()

# p3
# p3.sd

ggsave(
    (p2 + ggtitle("screen 2 vs\nscreen 1 with outliers [Effect sizes]")) +
        (p3 + ggtitle("screen 2 vs\nscreen 1 without outliers [Effect sizes]")) +
        (p1 + ggtitle("screen 1 with outliers vs\nscreen 1 without outliers [Effect sizes]")) +
        (p1.sd + ggtitle("screen 1 with outliers vs\nscreen 1 without outliers [SD effect sizes]")),
    filename = "/Users/karcher/tactic/output/comparing_screen_1_to_screen_2/mean_effect_size_screen_1_vs_screen_2.pdf",
    width = 12,
    height = 10
)
