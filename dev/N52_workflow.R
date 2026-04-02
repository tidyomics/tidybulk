library(tidyverse)
library(tidybulk)

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		aspect.ratio=1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

counts =
	readRDS("dev/N52.rds") %>%
	filter(symbol %>% is.na %>% `!`) %>%
	group_by(cell_type_formatted) %>%
	do(	keep_abundant((.), sample, symbol, count, minimum_proportion = 1/2)	) %>%
	ungroup() %>%
	group_by(symbol) %>%
	mutate(n = n()) %>%
	ungroup() %>%
	filter(n == 52) %>%
	select(-n)

# Create tidybulk object
tt = counts %>% tidybulk(sample, symbol, count)

# Normalise
tt_scaled = tt %>% scale_abundance()

# Plot densities
tt_scaled %>%
	pivot_longer(values_to = "count", names_to = "Normalisation", cols = c(count, `count_scaled`)) %>%
	ggplot(aes(count + 1, group=sample, color=cell_type_formatted)) +
	facet_grid(~Normalisation) +
	geom_density() + scale_x_log10() +
	my_theme

# Visualise MDS for cell types
tt_scaled %>%
	reduce_dimensions(assay = "counts_scaled", method="MDS", .dims = 3) %>%
	select(contains("Dim"), sample, CAPRA_TOTAL,  cell_type_formatted) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=cell_type_formatted))

# Vidualise MDS CAPRA
tt_scaled %>%
	group_by(cell_type_formatted) %>%
	do(reduce_dimensions((.), assay = "counts_scaled", method="MDS")) %>%
	select(contains("Dim"), sample, CAPRA_TOTAL,  cell_type_formatted) %>%
	distinct() %>%
	ggplot(aes(x=Dim1, y = Dim2, color=CAPRA_TOTAL)) +
	geom_point() +
	facet_wrap(~cell_type_formatted)
