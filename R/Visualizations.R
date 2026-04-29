Sys.setlocale("LC_TIME", "English")
library(pacman)
p_load(here,
       readxl,
       dplyr,
       lubridate,
       eurostat,
       tidyverse)

# ── 1. Downloading annual rates of change ────────────────────────────────────────
rates <- get_eurostat(
  id          = "prc_hicp_manr",
  time_format = "date",
  cache       = FALSE,
  filters     = list(
    geo    = "EA20",
    unit   = "RCH_A",
    coicop = c("CP00", "NRG", "FOOD", "IGD_NNRG", "SERV")
  )
)

# ── 2. Downloading item weights ───────────────────────────────────────────────────
weights <- get_eurostat(
  id          = "prc_hicp_inw",
  time_format = "date",
  cache       = FALSE,
  filters     = list(
    geo    = "EA20",
    coicop = c("CP00", "NRG", "FOOD", "IGD_NNRG", "SERV")
  )
)

# ── 3. Merging and computing contributions ────────────────────────────────────────
weights_clean <- weights %>%
  mutate(year = year(time)) %>%
  group_by(coicop, year) %>%
  summarise(weight = mean(values, na.rm = TRUE), .groups = "drop")

contributions <- rates %>%
  filter(coicop != "CP00",
         time >= as.Date("2002-01-01")) %>%
  mutate(year = year(time)) %>%
  left_join(weights_clean, by = c("coicop", "year")) %>%
  mutate(
    contribution = (weight / 1000) * values,
    label = recode(coicop,
                   "NRG"      = "Energy",
                   "FOOD"     = "Food",
                   "IGD_NNRG" = "Non-energy goods",
                   "SERV"     = "Services"
    )
  ) %>%
  group_by(time, label) %>%
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = factor(label, levels = c("Non-energy goods", "Services", "Energy", "Food")))

headline <- rates %>%
  filter(coicop == "CP00",
         time >= as.Date("2002-01-01")) %>%
  select(time, headline = values)

# ── 4. Plot ────────────────────────────────────────────────────────────────────
ecb_colours <- c(
  "Food"             = "#4472C4",
  "Energy"           = "#CC3399",
  "Non-energy goods" = "#70AD47",
  "Services"         = "#FFA500"
)

ggplot() +
  geom_bar(
    data     = contributions,
    aes(x = time, y = contribution, fill = label),
    stat     = "identity",
    position = "stack"
  ) +
  geom_line(
    data      = headline,
    aes(x = time, y = headline, linetype = "Headline inflation"),
    colour    = "black",
    linewidth = 0.7
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  scale_fill_manual(values = ecb_colours) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Headline inflation" = "solid")
  ) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%y",
    expand      = c(0.01, 0)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    title   = "Euro Area Headline Inflation and Component Contributions",
    x       = NULL,
    y       = "Percentage points",
    fill    = NULL,
    caption = paste0(
      "Source: Eurostat (prc_hicp_manr, prc_hicp_inw). ",
      "Note: Bars show weighted percentage point contributions to headline HICP inflation (year-on-year, %)."
    )
  ) +
  guides(
    fill     = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(hjust = 0.5, size = 12),
    plot.caption       = element_text(size = 7, hjust = 0, colour = "grey40"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA)
  )

ggsave(
  filename = "Output/Figures/hicp_contributions.png",
  width    = 10,
  height   = 5.5,
  dpi      = 300
)
