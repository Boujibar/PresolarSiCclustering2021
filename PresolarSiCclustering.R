# Author: Asmaa Boujibar
# email: aboujibar@carnegiescience.edu  asmaa.boujibar@gmail.com


# load library ------------------------------------------------------------
# If a library is not available yet click on "Tools" (top menu) then "install packages" and look for it to install
library(readxl)
library(janitor)
library(tidyr)
library(factoextra)
library(FactoMineR)
library(mclust)
library(plotly)
library(xlsx)
library(tidyverse)

# Clearing variables and plots
rm(list = ls(all.names = TRUE))
dev.off()


# Click on "Session" at the top menu, then "Set Working Directory" and "Choose directory".
# Choose the folder where the the database excel file is located. Make sure to include a folder named "Figures" in it.


######### Load data ---------------------------------------------------------------

data_path <- "PGD_SiC_2020-08-18.xlsx"

df <- read_excel(data_path, sheet = "PGD-SIC")

######### Select data ---------------------------------------------------------

# clean names
df_name <- df %>% 
  clean_names()

# Check names
names(df_name)


# Selection attributes and errors

df_select <- df_name %>% 
  dplyr::select(pgd_id,
                pgd_type,
                x12c_13c,
                err_12c_13c,
                x14n_15n,
                err_14n_15n,
                d_29si_28si,
                err_d_29si_28si,
                d_30si_28si,
                err_d_30si_28si,
  ) %>% 
  dplyr::rename(id = pgd_id,
                ntype = pgd_type,
                "12C/13C" = x12c_13c,
                "err12C/13C" = err_12c_13c,
                "14N/15N" = x14n_15n,
                "err14N/15N" = err_14n_15n,
                "29Si/28Si" = d_29si_28si,
                "err29Si/28Si" = err_d_29si_28si,
                "30Si/28Si" = d_30si_28si,
                "err30Si/28Si" = err_d_30si_28si)

# Making sure isotopic compositions are numeric values

df_select[, 3:10] <- sapply(df_select[, 3:10], as.numeric)
sapply(df_select, class)

# Print first lines

head(df_select)

# Removing M grains with large Si errors

df_select_large_Si_error <- subset(df_select, `err29Si/28Si` > 10 & `err30Si/28Si` > 10 & ntype == "M")

idslargeerr <- df_select_large_Si_error$id

for (ro in 1:dim(df_select_large_Si_error)[1]) {
  df_select <- subset(df_select, id != idslargeerr[ro])
}

# Select data for clustering

df_scenario <- df_select %>% 
  dplyr::filter(complete.cases(.[c("12C/13C", "30Si/28Si", "29Si/28Si", "14N/15N")])) %>%
  dplyr::select("id", "12C/13C", "30Si/28Si", "29Si/28Si", "14N/15N", "ntype")

# Excluding C and U grains 

df_scenario <- df_scenario %>%
  dplyr::filter(grepl("X|N|AB|M|Y|Z", df_scenario$ntype))
table(df_scenario$ntype)

# Excluding contaminated grains

df_scenario_contam <- subset(df_scenario, `12C/13C` < 93.56 & `12C/13C` > 88.87)
df_scenario_contam <- subset(df_scenario_contam, `14N/15N` < 339.94 & `14N/15N` > 248)
df_scenario_contam <- subset(df_scenario_contam, `30Si/28Si` < 50 & `30Si/28Si` > -50)
df_scenario_contam <- subset(df_scenario_contam, `29Si/28Si` < 50 & `29Si/28Si` > -50)

dim_df_contam <- dim(df_scenario_contam)
idscontam <- df_scenario_contam$id

for (ro in 1:dim_df_contam[1]) {
  df_scenario <- subset(df_scenario, id != idscontam[ro])
}

str(df_scenario)

########### Conversions -------------------------------------------------------------

C12_C13_S <- 89 
N14_15_S <- 440
N14_15_E <- 272 
Si29_28_0 <- 0.0506331
Si30_28_0 <- 0.0334744

# Transform delta Si into ratios

df_scenario <- df_scenario %>% 
  mutate(`29Si/28Si` = ((`29Si/28Si`/1000)+1) * Si29_28_0,
         `30Si/28Si` = ((`30Si/28Si`/1000)+1) * Si30_28_0
  )

# Remove negative 30Si

df_scenario <- df_scenario %>%
  dplyr::filter(
    `30Si/28Si` > 0)

# Calculate log of ratios

df_log <- df_scenario %>% 
  mutate(`12C/13C` = log10(`12C/13C`),
         `30Si/28Si` = log10(`30Si/28Si`),
         `29Si/28Si` = log10(`29Si/28Si`),
         `14N/15N` = log10(`14N/15N`)
  )

# scaling 

df_log_scale <- as.data.frame(scale(df_log[2:5], center = TRUE, scale = TRUE))

df_log_scale$Type <- df_log[[6]]

df_log_scale$Id <- df_log[[1]]

head(df_log_scale)


######## Cluster analysis with Mclust -------------------------------------------------

#Setting seed for reproducibility

set.seed(42)

# Choose the optimum mixed models and number of clusters

BIC <- mclustBIC(df_log_scale[1:4], G = 1:20)

#Actually run the model based clustering

mc <- Mclust(df_log_scale[1:4], x=BIC)

fviz_mclust(mc, 
            "BIC", 
            palette = "jco",
            legend = "bottom")

# Saving cluster results in dataframe

df_log_scale$cluster <- as.factor(mc$classification)

# Re-convert to original values then plots

df_original <- df_log %>% 
    mutate(`12C/13C` = 10^(`12C/13C`),
           `14N/15N` = 10^(`14N/15N`),
           `30Si/28Si` = 1000*((10^(`30Si/28Si`)/Si30_28_0)-1),
           `29Si/28Si` = 1000*((10^(`29Si/28Si`)/Si29_28_0)-1))

df_original <- df_original %>% 
    dplyr::rename(d30Si = "30Si/28Si",
                  d29Si = "29Si/28Si")

# Assign new cluster label

df_original$cluster <- NA
df_original$cluster[which(df_log_scale$cluster == 1)] <- 5
df_original$cluster[which(df_log_scale$cluster == 2)] <- 4
df_original$cluster[which(df_log_scale$cluster == 3)] <- 9
df_original$cluster[which(df_log_scale$cluster == 4)] <- 1
df_original$cluster[which(df_log_scale$cluster == 5)] <- 2
df_original$cluster[which(df_log_scale$cluster == 6)] <- 8
df_original$cluster[which(df_log_scale$cluster == 7)] <- 7
df_original$cluster[which(df_log_scale$cluster == 8)] <- 3
df_original$cluster[which(df_log_scale$cluster == 9)] <- 6
df_original$cluster <- as.factor(df_original$cluster)

# Plots

cbPalette <- c("#10BFBF", "#00FF00", "#FFFF00", "#FF9900", "#1200FF", "#999999", "#FF0000", "#FF00CC", "#136330")
symbols <- c(22, 25, 17, 12, 9, 8, 10, 11, 18)

ggplot(df_original) +
  geom_point(aes(x = `12C/13C`, y = `14N/15N`, shape = cluster, color = cluster), size = 2.5) +
  stat_ellipse(aes(x = `12C/13C`, y= `14N/15N`, color = cluster),type = "norm") +
  scale_shape_manual(values = symbols) +
  scale_colour_manual(values = cbPalette) +
  scale_x_log10(limits = c(1,5000),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(limits = c(5,19950),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_vline(xintercept= `C12_C13_S`)+
  geom_hline(yintercept= `N14_15_S`)+
  geom_hline(yintercept= `N14_15_E`)+
  labs(x = bquote(""^12*"C/"^13*"C"), y = bquote(""^14*"N/"^15*"N")) +
  theme_bw() +
  theme(aspect.ratio = 0.9,
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        plot.margin = ggplot2::margin(l=0.5, t=0.5, b=0.3, r=0.5, unit = 'cm'),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()
  )
dev.copy2pdf(file = file.path("Figures", paste0("C-N-clusters", ".pdf")), width = 6, height = 6, useDingbats=FALSE)

ggplot(df_original) +
  geom_point(aes(x = `d30Si`, y = `d29Si`, shape = cluster, color = cluster), size = 2.5) +
  stat_ellipse(aes(x = `d30Si`, y= `d29Si`, color = cluster),type = "norm") +
  scale_colour_manual(values = cbPalette) +
  scale_x_continuous(sec.axis = sec_axis(~ . * 1), limits = c(-750,600)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 1), limits = c(-600,250)) +
  geom_vline(xintercept= `Si30_28_0`)+
  geom_hline(yintercept= `Si29_28_0`)+
  scale_shape_manual(values = symbols) +
  labs(x = expression(paste(delta, "("^{30}, "Si/"^{28}, "Si) (???)")), y = expression(paste(delta, "("^{29}, "Si/"^{28}, "Si) (???)"))) +
  theme_bw() +
  theme(aspect.ratio = 0.9,
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        plot.margin = ggplot2::margin(l=0.5, t=0.5, b=0.3, r=0.5, unit = 'cm'),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.text.x = element_text(size = 14, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),
        axis.text.y = element_text(size = 14, margin = unit(c(t = 2.5, r = 3, b = 3, l = 3), "mm")),
        axis.text.x.top = element_blank(), axis.text.y.right = element_blank()
  )
dev.copy2pdf(file = file.path("Figures", paste0("Si-Si_clusters", ".pdf")), width = 6, height = 6, useDingbats=FALSE)


ggplot(df_original) +
  geom_point(aes(x = `d30Si`, y = `d29Si`, shape = cluster, color = cluster), size = 2.5) +
  stat_ellipse(aes(x = `d30Si`, y= `d29Si`, color = cluster),type = "norm") +
  scale_colour_manual(values = cbPalette) +
  scale_x_continuous(sec.axis = sec_axis(~ . * 1), limits = c(-100,200)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 1), limits = c(-100,240)) +
  geom_vline(xintercept= `Si30_28_0`)+
  geom_hline(yintercept= `Si29_28_0`)+
  scale_shape_manual(values = symbols) +
  labs(x = expression(paste(delta, "("^{30}, "Si/"^{28}, "Si) (???)")), y = expression(paste(delta, "("^{29}, "Si/"^{28}, "Si) (???)"))) +
  theme_bw() +
  theme(aspect.ratio = 0.9,
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        plot.margin = ggplot2::margin(l=0.5, t=0.5, b=0.3, r=0.5, unit = 'cm'),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.text.x = element_text(size = 14, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),
        axis.text.y = element_text(size = 14, margin = unit(c(t = 2.5, r = 3, b = 3, l = 3), "mm")),
        axis.text.x.top = element_blank(), axis.text.y.right = element_blank()
  )
dev.copy2pdf(file = file.path("Figures", paste0("Si-Si-zoom_clusters", ".pdf")), width = 6, height = 6, useDingbats=FALSE)


# Probability in each cluster

prob <- data.frame(mc$z)

prob <- prob %>% 
  dplyr::select(X1,
                X2,
                X3,
                X4,
                X5,
                X6,
                X7,
                X8,
                X9
  ) %>% 
  dplyr::rename(X1 = X5,
                X2 = X4,
                X3 = X9,
                X4 = X1,
                X5 = X2,
                X6 = X8,
                X7 = X7,
                X8 = X3,
                X9 = X6,
                
  )

prob <- prob[, c(5, 4, 9, 1, 2, 8, 7, 3, 6)]

# Probability in assigned cluster

probUnique <- vector(length=length(mc$classification))
for (i in 1:length(mc$classification)) { probUnique[i] <- max(prob[i, ])}

df_original$prob <- probUnique

# Save data in excel file

write.xlsx(df_original, "Figures\\ResultsClusters.xlsx", sheetName = "data and cluster", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(prob, "Figures\\ResultsClusters.xlsx", sheetName = "prob in each cluster", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


######## Principal component analysis -------------------------------------------------

pca_out<- PCA(df_log_scale[1:4], scale.unit = FALSE, ncp = 5, graph = FALSE)

pcaPalette <- c("#1200FF", "#FF9900", "#136330", "#10BFBF", "#00FF00", "#FF00CC", "#FF0000", "#FFFF00", "#999999")
pcasymbols <- c(9, 12, 18, 22, 25, 11, 10, 17, 8)

fviz_pca_biplot(pca_out, 
                label="var", 
                habillage=as.factor(mc$classification),
                pointsize = 3,
                palette=pcaPalette,
                addEllipses=TRUE, 
                col.var = "black",
                ellipse.level=0.95) +
  theme_bw() + scale_shape_manual(values=pcasymbols)+
  theme(legend.position = "none",
        aspect.ratio = 0.9,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(file = file.path("Figures", paste0("pca_cluster", ".pdf")), width = 8, height = 8, useDingbats=FALSE)

# Variable contribution

var <- get_pca_var(pca_out)

fviz_contrib(pca_out, choice = "var", axes = 1)
dev.copy2pdf(file = file.path("Figures", paste0("Contrib_var1_pca_DB4", ".pdf")), width = 8, height = 8, useDingbats=FALSE)

fviz_contrib(pca_out, choice = "var", axes = 2)
dev.copy2pdf(file = file.path("Figures", paste0("Contrib_var2_pca_DB4", ".pdf")), width = 8, height = 8, useDingbats=FALSE)


