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

library("pasilla")
pasCts <- system.file("extdata",
											"pasilla_gene_counts.tsv",
											package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
											 "pasilla_sample_annotation.csv",
											 package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]

# Create tidybulk object
tt =
	cts %>%
	as_tibble(rownames = "transcript") %>%
	pivot_longer(names_to = "sample", values_to = "count", cols=-transcript) %>%
	left_join(
		coldata %>%
			as_tibble(rownames = "sample") %>%
			mutate(sample = gsub("fb", "", sample))
	) %>%

	# attach symbol
	left_join(
		tibble(
			transcript = (.) %>% distinct(transcript) %>% pull(1),
			symbol = flybaseR::id.converter((.) %>% distinct(transcript) %>% pull(1), symbols = TRUE, convert.into = "genes")
		)
	) %>%

	# tidybulk
	mutate_if(is.character, as.factor) %>%
	tidybulk(sample, transcript, count)

# Normalise
tt_scaled = tt %>% scale_abundance()

# Plot densities
p1 =
	tt_scaled %>%
	rename(`Raw counts` = count, `Scaled counts` = count_scaled) %>%
	pivot_longer(values_to = "Counts", names_to = "Normalisation", cols = c(`Raw counts`, `Scaled counts`)) %>%
	ggplot(aes(Counts + 1, group=sample, color=type)) +
	facet_grid(~Normalisation) +
	geom_density() +
	scale_x_log10() +
	scale_color_brewer(palette = "Set1") +
	ylab("Density") +
	my_theme

tt_mds = tt_scaled %>% 	reduce_dimensions(assay = "counts_scaled", method="MDS", .dims = 3)

# Visualise MDS for cell types
p2 =
	tt_mds %>%
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	rename(`Dimension 1` = Dim1, `Dimension 2` = Dim2, `Dimension 3` = Dim3) %>%
	GGally::ggpairs(
		columns = 1:3,
		ggplot2::aes(colour=condition),
		upper = list(continuous = GGally::wrap("cor", size = 3))
	) +
	#scale_color_brewer(palette = "Set1") +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		legend.position="bottom",
		strip.background = element_blank()
	)

# Visualise MDS for cell types
tt_mds %>%
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=type), upper = list(continuous = GGally::wrap("cor", size = 3))) +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		legend.position="bottom",
		strip.background = element_blank()
	)

# Adjust abundance for type, and check MDS densities
tt_adj = tt_mds %>%	adjust_abundance(~ condition + type)

p3 =
	tt_adj %>%
	filter( `count_scaled_adjusted` %>% is.na %>% `!`) %>%
	reduce_dimensions(assay = "counts_scaled_adjusted", method="MDS", .dims = 3) %>%

	# Plot
	select(contains("Dim"), sample, type,  condition ) %>%
	distinct() %>%
	pivot_longer(names_to = "Dim", values_to = ".value", cols = contains("Dim")) %>%
	separate(Dim, c("Dim", "Adj"), sep="\\.") %>%
	mutate(Dim = gsub("Dim", "Dimension ", Dim)) %>%
	mutate(Adj = ifelse(Adj == "y", "Non-adjusted", "adjusted") %>% factor(labels = c("Non-adjusted", "adjusted"))) %>%
	rename(Method = type, `Treatment status` = condition) %>%
	pivot_longer(names_to = "variation", values_to = "Type of variation", cols = c(Method, `Treatment status`)) %>%
	mutate(`Type of variation` = stringr::str_to_sentence(`Type of variation`)) %>%
	ggplot(aes(y = .value, x = variation, fill = `Type of variation`)) +
	geom_boxplot() +
	scale_fill_brewer(palette="Set1") +
	facet_grid(Adj ~ Dim ) +
	xlab("Annotation") +
	ylab("Dimension value") +
	my_theme

# DE analyses
tt_test = tt_adj %>% test_differential_abundance(~ condition + type)

# MA plot
p4 =
	tt_test %>%
	filter(!`lowly_abundant`) %>%
	mutate(`Differentially abundant` = FDR<0.05 & abs(logFC) >=2) %>%
	distinct(transcript, symbol, logCPM, logFC, `Differentially abundant`) %>%
	mutate(symbol = ifelse(`Differentially abundant` & abs(logFC)>2, as.character(symbol), NA)) %>%
	ggplot(aes(x = logCPM, y = logFC, label=symbol)) +
	geom_point(aes(color = `Differentially abundant`, size = `Differentially abundant`, alpha=`Differentially abundant`)) +
	ggrepel::geom_text_repel() +
	scale_color_manual(values=c("black", "#e11f28")) +
	scale_size_discrete(range = c(0, 2)) +
	xlab("Logarithm of count per million") +
	ylab("Logarithm of fold change") +
	my_theme

# Plot top genes for raw, scaled and ajusted counts
p5 =
	tt_test %>%
	inner_join( (.) %>% distinct(symbol, PValue) %>% arrange(PValue) %>% head(6)) %>%
	pivot_longer(names_to = "Stage", values_to = "Counts", cols = starts_with("count")) %>%
	mutate(Stage = ifelse(Stage == "count", "Raw", Stage)) %>%
	mutate(Stage = ifelse(Stage == "count_scaled", "Scaled", Stage)) %>%
	mutate(Stage = ifelse(Stage == "count_scaled_adjusted", "Adjusted", Stage)) %>%
	mutate(Stage = factor(Stage, levels = c("Raw", "Scaled", "Adjusted"))) %>%
	rename(`Treatment status` = condition) %>%
	mutate(`Treatment status` = stringr::str_to_sentence(`Treatment status`)) %>%
	mutate(symbol = toupper(symbol)) %>%
	ggplot(aes(x = Stage, y = Counts + 1, fill = `Treatment status`)) +
	geom_boxplot() +
	#ggpol::geom_boxjitter(jitter.width = 0.1,  jitter.size = 0.2	) +
	facet_wrap(~symbol) +
	scale_y_log10() +
	scale_fill_brewer(palette = "Set1") +
	my_theme +
	theme(axis.text.x =element_text(angle = 60, hjust = 1, vjust = 1))

# Heatmap
pdf("dev/pasilla_p6.pdf", width = "91.5mm", height = "91.5mm",)
tt_test %>%
filter(FDR < 0.05 & abs(logFC) > 2) %>%
rename(`Counts\nscaled\nadjusted\n` = count_scaled_adjusted) %>%
tidyHeatmap::heatmap(		.horizontal = sample,		.vertical = symbol,		.abundance = `Counts\nscaled\nadjusted\n` ,
	annotation = c(condition, type),
	log_transform = TRUE
)
dev.off()

# Integrate all plots
ggsave(	"dev/pasilla_p1.pdf",	plot = p1,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)
ggsave(	"dev/pasilla_p2.pdf",	plot = p2,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 , heigh = 183/2,	limitsize = FALSE)
ggsave(	"dev/pasilla_p3.pdf",	plot = p3,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)
ggsave(	"dev/pasilla_p4.pdf",	plot = p4,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)
ggsave(	"dev/pasilla_p5.pdf",	plot = p5,	useDingbats=FALSE,	units = c("mm"),	width = 183/2 ,	limitsize = FALSE)
