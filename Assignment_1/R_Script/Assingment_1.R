## 1 - Packages used ----
#Uncomment and run if you need to install the package or load the library 
#installed packages for bold data acquisition
#install.packages("BiocManager")
#BiocManager::install("sangerseqR")
#install.packages("remotes")
#remotes::install_github("ropensci/bold")
#library(bold)
#installed package for world map
#install.packages("rworldmap")
library(tidyverse)
library(rworldmap)
library(rworldxtra)
library(viridis)
conflicted::conflict_prefer("filter", "dplyr")

## 2 - BOLD data acquiring and saving ----
# df_bold = bold_specimens(taxon = c("Cricetidae", "Chromadorea"))
# write_tsv(df_bold, "Bold.tsv")

## 3 - Reading Data ----

# Remember to set your working directory

# Reading in data and Selecting relevant columns for the study

df_bold <- read_tsv("../data/Bold.tsv") %>%
  select(class_name, species_name, lat, lon) %>%
  drop_na() %>%
  filter(!if_any(where(is.character), ~ . == "" | . == " ")) # rather than simply checking, remove rows with empty strings or space

# It is a bit confusing that Cricetidae are called "Mammalia" in the tsv file and in the rest of the code. For clarity of code comprehension I am changing these entries to "Cricetidae", which carries over to the rest of the code
df_bold$class_name[df_bold$class_name == "Mammalia"] <- "Cricetidae"

## 4 - Exploring Data and Checking Formatting ----
head(df_bold, 3) # more concise way of showing colnames, classes and variable examples

# Summary of of overall data structure
summary(df_bold)

## 5 - Map representing worldwide distribution of species of both Chromadorea and Cricetidae species----
# Tibble for Chromadorea species count grouped based on latitude and longitude
species_latlon <- function(df = df_bold, class) {
  df %>%
    filter(class_name == class) %>%
    mutate(lat = round(lat), lon = round(lon)) %>%
    group_by(lat, lon) %>%
    summarise(
      species_count = n_distinct(species_name),
      .groups = "drop"
    )
}
df_chro_ll <- species_latlon(class = "Chromadorea")
df_cric_ll <- species_latlon(class = "Cricetidae")

# note that most of the count data is confined to a small range
summary(df_chro_ll$species_count)
quantile(df_chro_ll$species_count, 0.995) # 99.5% of the Chromadorea data is in the range 1-50, Cricetidae is confined to an even smaller range (code not shown).

# Get world map data and converting it to data frame for ggplot2
df_world_map <- getMap(resolution = "high") %>%
  fortify()
# Plot Chromadorea species geographical distribution on world map
ggplot() +
  geom_polygon(
    data = df_world_map, aes(x = long, y = lat, group = group),
    fill = "lightgray", color = "gray"
  ) +
  geom_point(data = df_chro_ll, aes(x = lon, y = lat, color = species_count), shape = 10, size = 0.9, alpha = 0.6) +
  geom_point(data = df_cric_ll, aes(x= lon, y = lat, color = species_count), shape = 25, size = 0.9, alpha = 0.6) +
  scale_color_viridis_c(option = "turbo", name = "Number of Species", transform = "log2", limits = c(1, 50)) + # including 99.5% of the data, on log scale for visibility
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  labs(
    title = "Distribution of Chromadorea and Cricetidae Worldwide",
    subtitle = "Number of Species transformed to log2 scale",
    x = "Longitude (deg)",
    y = "Latitude (deg)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.3, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

# Exporting the Plot with size writable for A4 paper
#ggsave("Chromadorea Cricetidae Species Geographical Distribution.PNG", dpi = 600)

## 6 - Correlation between distance from equator (degree) and species richness of Chromadorea and Cricetidae ----
# Tibble with the species count of each targeted taxa grouped based on distance from equator (degree)
df_dist_eq <- df_bold %>%
  mutate(dist_eq_deg = round(abs(lat))) %>% # Converting latitude into distance from equator (degrees)
  group_by(dist_eq_deg) %>%
  summarise(
    Cricetidae = n_distinct(species_name[class_name == "Cricetidae"]),
    Chromadorea = n_distinct(species_name[class_name == "Chromadorea"])
  )

# Create function for Shapiro-Wilk test
shapiroWilk <- function(name, df = df_dist_eq, log = FALSE) {
  if (log == TRUE) {
    temp <- df %>%
      filter(name != 0) %>%
      pull(name)
    col <- ifelse(temp == 0, 0, log(temp))
    return (col %>%
              shapiro.test %>%
              print())}
  else df %>%
    filter(name != 0 ) %>%
    pull(name) %>%
    shapiro.test %>%
    print()
}

# Shapiro-Wilk test of normality to determine if the Cricetidae data follows normal distribution
shapiroWilk("Cricetidae") # W = 0.93 (very good fit) and p < 0.001
shapiroWilk("Chromadorea") # W = 0.72 (somewhat good fit) and p < 4e-11

# Combining Cricetidae and Chromadorea species count into a single column
df_dist_eq_long <- df_dist_eq %>%
  pivot_longer(
    cols = c(Cricetidae, Chromadorea),
    names_to = "group",
    values_to = "species_richness"
  )

# Scatter plot represents correlation between species richness (Chromadorea and Cricetidae) and distance from the equator
ggplot(df_dist_eq_long, aes(x = dist_eq_deg, y = species_richness)) +
  geom_point() +
  labs(
    title = "Correlation between species richness vs distance from equator",
    x = "Distance from Equator (degrees)",
    y = "Species Richness"
  ) +
  facet_wrap(~group, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray100", size = 0.25),
    strip.text = element_text(size = 12, face = "bold")
  )

# Exporting the Plot with size suitable for A4 paper
# ggsave("Correlation between species richness vs distance from equator.PNG", dpi = 600)

# with log-transformed data the test statistics are both good fits
shapiroWilk("Cricetidae", log = TRUE) # W = 0.91 (very good fit) and p < 2e-05
shapiroWilk("Chromadorea", log = TRUE) # W = 0.96 (near-perfect fit) and p < 0.025

# extract r^2 from quadratic model for both classes, to be used on plot downstream
r2_cricetidae <- lm(data = df_dist_eq,
                    ifelse(Cricetidae == 0, 0, log(Cricetidae)) ~ poly(dist_eq_deg, 2)) %>% 
  summary %>%
  .$r.squared
r2_chromadorea <- lm(data = df_dist_eq,
                     ifelse(Chromadorea == 0, 0, log(Chromadorea)) ~ poly(dist_eq_deg, 2)) %>% 
  summary %>%
  .$r.squared

# Visualize transformed data
ggplot(data = df_dist_eq) +
  # points for Cricetidae
  geom_point(
    aes(
      x=dist_eq_deg, 
      y = ifelse(Cricetidae == 0,
                 Cricetidae,
                 log(Cricetidae)),
      shape = "Cricetidae", color = "Cricetidae")) +
  # polynomial regression fitted to log-transformed Cricetidae count data
  geom_smooth(
    data = df_dist_eq,
    size = 0.5,
    aes(
      x = dist_eq_deg,
      y = ifelse(Cricetidae == 0, Cricetidae, log(Cricetidae)),
      color = "Cricetidae"
    ),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = FALSE # don't need to show standard error bars
  ) +
  # points for Chromadorea
  geom_point(
    aes(
      x=dist_eq_deg, 
      y = ifelse(Chromadorea == 0,
                 Chromadorea,
                 log(Chromadorea)),
      shape = "Chromadorea", color = "Chromadorea")) +
  scale_color_manual(
    values = c("Cricetidae" = "blueviolet", "Chromadorea" = "darkgoldenrod3"),
    name = "Species"
  ) +
  # polynomial regression fitted to log-transformed Chromadorea count data
  geom_smooth(
    data = df_dist_eq,
    size = 0.5,
    aes(
      x = dist_eq_deg,
      y = ifelse(Chromadorea == 0, Chromadorea, log(Chromadorea)),
      color = "Chromadorea"
    ), 
    method = "lm",
    formula = y ~ poly(x, 2),
    se = FALSE # don't need standard error bars shown
  ) +
  scale_shape_manual(
    values = c("Cricetidae" = 0, "Chromadorea" = 1),
    name = "Species"
  ) +
  # Add R^2 for Cricetidae quadratic regression
  annotate(
    "text", x = 63, y = 4.35, label = paste0("Cricetidae R² = ", round(r2_cricetidae, 3)),
    color = "black", size = 2.5, hjust = 0
  ) +
  # Add R^2 for Chromadorea quadratic regression
  annotate(
    "text", x = 63, y = 4.20, label = paste0("Chromadorea R² = ", round(r2_chromadorea, 3)),
    color = "black", size = 2.5, hjust = 0
  ) +
  labs(
    title = "Species count of Cricetidae and Chromadorea vs. distance from the equator",
    subtitle = "Quadratic regression lines fitted overtop data",
    x = "Distance from Equator (deg)",
    y = "ln(species count)") +
  theme(
    plot.title = element_text(hjust = 0.3, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  ) +
  theme_bw()

# save plot
#ggsave("Species count vs distance.PNG", dpi = 600)

# Although the distance from the equator appears to follow a normal distribution for both Cricetidae and Chromadorea, there does not appear to be a linear relationship between species count and distance from the equator. As such, pearson's product moment correlation is not an appropriate correlation test. Either spearman or kendall's rank correlation should be used- but which of the two?

## 7 - Rank Correlation ---- 
# Spearman's Rank Correlation to determine if there is a correlation between latitude and Cricetidae species richness.


spearman_cric <- cor.test(df_dist_eq$dist_eq_deg[df_dist_eq$Cricetidae != 0], df_dist_eq$Cricetidae[df_dist_eq$Cricetidae != 0], method = "spearman") %>%
  print()

# ==> Result Spearman's Rank Correlation
# Rho (ρ) value = -0.46 ==> moderate negative correlation (with the bounds of [-1,1])
# p-value = 2.567e-05 ==> correlation is statistically significant

# Running Spearman's rank correlation test returns the warning that there are ties in the data. With only a few tied data points, this would not affect the interpretability of the results. Checking how many ties there are in the data:
paste0("There are ", sum(duplicated(df_dist_eq$Cricetidae)), " ties in Cricetidae counts.") %>% cat
paste0("There are ", sum(duplicated(df_dist_eq$Chromadorea)), " ties in Chromadorea counts.") %>% cat

# With just 80 rows, 47.5% to 61.25% of the data has ties. As such, Kendall's Tau will be more interpretable than Spearman's Rho for so much tied data:
Kendall_cor_Cric <- cor.test(df_dist_eq$dist_eq_deg, df_dist_eq$Cricetidae, method = "kendall") %>%
  print() # tau = -0.355 (slight-moderate negative correlation), p < 4e-06. 
Kendall_cor_Chro <- cor.test(df_dist_eq$dist_eq_deg, df_dist_eq$Chromadorea, method = "kendall") %>%
  print() # tau = -0.19 (very slight negative correlation), p < 0.025 

# There is a statistically significant moderate negative correlation between the distance from the equator (degrees) and species richness
# however, in data with discrete random variables (e.g. species counts) the bounds of Rho and Tau change -> Denuit, Mesfioui & Trufin wrote a statistical proof showing just this, with the real-world implication being that in some cases a small Rho or Tau actually represents a strong correlation (Denuit, M., Mesfioui, M., & Trufin, J. (2019). Concordance-based predictive measures in regression models for discrete responses. Scandinavian Actuarial Journal, 1–13. doi:10.1080/03461238.2019.1624274)
# Kendall's tau is represented by 𝛕[Y,X] = P(Y1 - Y2)(X1-X2) > 0] - P(Y1-Y2)(X1-X2) < 0] i.e. the probability that values are not tied minus the probability that they are tied. cor.test() assumes bounds of [-1,1] but Denuit, Mesfioui & Trufin show that this is not necessarily the case, especially in data with many observations.
# The upper bound can be derived as 𝛕[Y,X] ≤ 2*min{E[F_Y(Y^-)], E[F_X(X^-)]} where F_Y(Y^-) and F_X(X^-) are the CDFs of Y and X evaluated at the left limit.

## 8 - Calculating adjusted Tau ----
#Applying this to the dist_eq_deg and Cricetidae columns (Y and X, respectively):
# 1. calculate E[F_Y(Y^-)]
pmf_dist <- as.data.frame(table(df_dist_eq$dist_eq_deg)/80) # probability mass function, the probability of getting any one row should be the same as they are indices
colnames(pmf_dist) <- c("Obs", "Prob")
pmf_dist_lower <- c(0, cumsum(pmf_dist$Prob)[-nrow(pmf_dist)]) # cumulative probabilities up to i-1
cdfDist <- sum(pmf_dist$Prob * pmf_dist_lower) # E[F_Y(Y^-)] = 0.49375

# 2. calculate E[F_X(X^-)]
pmf_Cric <- as.data.frame(table(df_dist_eq$Cricetidae)/80, stringsAsFactors = F) # creating a table with each value that appears in the data, and a probability based on how frequently that value appears
colnames(pmf_Cric) <- c("Obs", "Prob")
pmf_Cric$Obs <- as.numeric(pmf_Cric$Obs)
pmf_Cric_lower <- c(0, cumsum(pmf_Cric$Prob)[-nrow(pmf_Cric)]) # ""
cdfCric <- sum(pmf_Cric$Prob * pmf_Cric_lower) # 0.4745

# and for Chromadorea
pmf_Chro <- as.data.frame(table(df_dist_eq$Chromadorea)/80, stringsAsFactors = F) # creating a table with each value that appears in the data, and a probability based on how frequently that value appears
colnames(pmf_Chro) <- c("Obs", "Prob")
pmf_Chro$Obs <- as.numeric(pmf_Chro$Obs)
pmf_Chro_lower <- c(0, cumsum(pmf_Chro$Prob)[-nrow(pmf_Chro)]) # ""
cdfChro <- sum(pmf_Chro$Prob * pmf_Chro_lower) # 0.4745

# 3. calculate true upper bound for Kendall's tau
Cric_bound <- 2*min(cdfDist, cdfCric)
Chro_bound <- 2*min(cdfDist, cdfChro)

# adjust Tau for true upper bound
unadjusted_tau_Cric <- Kendall_cor_Cric$estimate[[1]] 
adjusted_tau_Cric <- unadjusted_tau_Cric/Cric_bound # now representative of Tau in [-1,1] bounds
unadjusted_tau_Chro <- Kendall_cor_Chro$estimate[[1]] 
adjusted_tau_Chro <- unadjusted_tau_Chro/Chro_bound # now representative of Tau in [-1,1] bounds

paste0("The adjusted tau for Cricetidae is ", round(adjusted_tau_Cric, 2), " from -0.35") %>% cat # adjusted is about 2% greater than unadjusted
paste0("The adjusted tau for Chromadorea is ", round(adjusted_tau_Chro, 2), " from -0.19") %>% cat # adjusted is about 1% greater than unadjusted

## 9 - Correlation between Cricetidae and Chromadorea species distributions ----
# Kendall's Rank Correlation -- To determine if there is a correlation between Chromadorea and Circetidae species richness
correlation_result <- cor.test(df_dist_eq$Cricetidae, df_dist_eq$Chromadorea, method = "kendall") %>%
  print()

# ==> Result Spearman's Rank Correlation
# Tau (𝛕) value = 0.381 ==> moderate positive correlation
# p-value = 1.2e-06 ==> correlation is  statistically significant
# There is statistically significant moderate positive correlation between the geographical distributions of Cricetidae and Chromadorea
