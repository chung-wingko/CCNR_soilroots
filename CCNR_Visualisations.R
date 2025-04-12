### Set up ###

#clear workspace
rm (list =ls())
setwd("/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Data Analysis/Manuscript")

#load libraries
library(readr) #read_csv
library(dplyr) #data cleaning
library(tidyr)
library(plyr) #ddply function
library(tibble) #data cleaning
library(ggplot2) #plots
library(plotrix) #standard error
library(corrplot) #correlation matrix
library(MetBrewer) #colour palette
library("GGally") #correlogram
library(FactoMineR) #PCA
library(factoextra) #ggplot2 based visualisation
library(corrplot) #visualisation of correlation matrix
library(ggpubr) #graph aesthetics
library(cowplot) #add density plots to axes
library(lme4) #mixed models
library(lmerTest) #p values
library(LMERConvenienceFunctions) #diagnostic plots
library(emmeans) #posthoc tukey test
library(ggeffects) #predicted EMMs for continuous variable
library(ggbreak) #insert axis break
library(performance) #diff diagnostic plots
library(patchwork) #annotate multiple plots
library(stringr) #RDA
library(ggforce) #abline in lmer
library(Hmisc) #p values for correlation matrices
library(vegan) #dataset & NMDS & RDA
library(ggvegan) #ggplot for vegan
library(vegan3d) #visualise NMDS 3D
library(goeveg) #NMDS screeplot
library(ggrepel) #labels
library(ggord) #creating ordination plots with ggplot2 - RDA
library(CCA) #CCA
library(VennDiagram) #for variance partitioning
library(gridExtra)
library(agricolae) #tukey test
library(ggcorrplot)


#### Import data 

soildata <- read_csv("/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Cleaned Data/Soil_final.csv")
soildata$FT <- factor(soildata$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

rootdata <- read_csv("/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Cleaned Data/Root_final.csv")
rootdata$FT <- factor(rootdata$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

CCNRalldata <- soildata %>% full_join(rootdata)
CCNRalldata$FT <- factor(CCNRalldata$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

ft <- soildata %>% dplyr::select(Plot, FT) %>% unique()
ft$FT <- as.factor(ft$FT)
ft$FT <- factor(ft$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

ft_colors <- c("#dd5129", "#0f7ba2", "#43b284")


#### Plot measurements

summary(soildata$pH)
boxplot(soildata$pH)
hist(soildata$pH)

summary(rootdata$rN_P)
hist(rootdata$rN_P)

summary(soildata$sP)
summary(soildata$PO4_NF)

tree_counts <- read_csv("/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Cleaned Data/Tree_Counts.csv") 
tree_counts <- tree_counts %>% column_to_rownames("Plot")

# tree species richness
sprich <- tree_counts
sprich$sprich <- rowSums(sprich!= 0)
sprich <- sprich %>% rownames_to_column("Plot")
sprich <- sprich %>% dplyr::select(Plot, sprich) %>% full_join(ft)
sprich$FT <- as.factor(sprich$FT)
richanova <- aov(sprich ~ FT, data = sprich)
tukey <- HSD.test(richanova, "FT", group = TRUE)
print(tukey)
summary(richanova)
ggplot(sprich, aes(x = FT, y = sprich)) + geom_boxplot() + theme_bw()
aggregate(sprich ~ FT, data = sprich, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

# Shannon's diversity
shannon <- tree_counts
shannon$shannon <- diversity(shannon, index = "shannon")
shannon <- shannon %>% rownames_to_column("Plot") %>% dplyr::select(Plot, shannon) %>% full_join(ft)
shannon$FT <- as.factor(shannon$FT)
shannonanova <- aov(shannon ~ FT, data = shannon)
tukey <- HSD.test(shannonanova, "FT", group = TRUE)
print(tukey)
summary(shannonanova)
ggplot(shannon, aes(x = FT, y = shannon)) + geom_boxplot() + theme_bw()
aggregate(shannon ~ FT, data = shannon, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

# stem density
stemdens <- tree_counts
stemdens$stemdens <- rowSums(stemdens)
#convert to stem density/ha
stemdens$stemdens <- stemdens$stemdens *10000/400
stemdens <- stemdens %>% rownames_to_column("Plot") %>% dplyr::select(Plot, stemdens) %>% full_join(ft)
stemdens$FT <- as.factor(stemdens$FT)
stemdensanova <- aov(stemdens ~ FT, data = stemdens)
tukey <- HSD.test(stemdensanova, "FT", group = TRUE)
print(tukey)
summary(stemdensanova)
ggplot(stemdens, aes(x = FT, y = stemdens)) + geom_boxplot() + theme_bw()
aggregate(stemdens ~ FT, data = stemdens, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#basal area
BAanova <- aov(BA ~ FT, data = soildata)
tukey <- HSD.test(BAanova, "FT", group = TRUE)
print(tukey)
summary(BAanova)
ggplot(soildata, aes(x = FT, y = BA)) + geom_boxplot() + theme_bw()
aggregate(BA ~ FT, data = soildata, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

plot_measurements <- sprich %>% full_join(shannon) %>% full_join(stemdens) %>% full_join(soildata) %>% dplyr::select(Plot, FT, Direction, sprich, shannon, stemdens, BA)

CCNRalldata <- CCNRalldata %>% full_join(plot_measurements)


#### All means

alllength <- as.data.frame(colSums(!is.na(CCNRalldata))) %>% tibble::rownames_to_column("variable")
alllength <- alllength[c(4:55),]

all_means <- CCNRalldata %>% 
  pivot_longer(-c(Plot, FT, Direction), names_to = "variable") %>% 
  group_by(FT) %>% full_join(alllength) %>% dplyr::rename(n = 'colSums(!is.na(CCNRalldata))') 

all_means_FT <- all_means %>%
  group_by(FT, variable) %>%
  dplyr::reframe(mean = mean(value, na.rm = TRUE),
                sd = sd(value, na.rm = TRUE),
                se = sd / sqrt(n))
all_means_FT <- unique(all_means_FT)

total_means <- all_means %>%
  group_by(variable) %>%
  dplyr::reframe(mean = mean(value, na.rm = TRUE),
                sd = sd(value, na.rm = TRUE),
                se = sd / sqrt(n))
total_means <- unique(total_means)
total_means$FT <- "Overall"
total_means <- total_means %>% dplyr::select(FT, variable, mean, sd, se)
options(scipen = 999)
allsummary <- rbind(total_means, all_means_FT)
allsummary <- allsummary %>% dplyr::select(-sd)
allsummary <- allsummary %>% pivot_wider(names_from = FT, values_from = c(mean, se))
#allsummary <- allsummary %>% arrange(variable)

write.csv(allsummary, "mean_soil_values.csv")

soildata <- soildata %>% dplyr::select(-sC_percent, -sN_percent)
rootdata <- rootdata %>% dplyr::select(-rC_percent, -rN_percent)
CCNRalldata <- CCNRalldata %>% dplyr::select(-sC_percent, -sN_percent, -rC_percent, -rN_percent)


#### Correlograms

#root vs soil nutrients
CCNRalldata_corr <- CCNRalldata %>% dplyr::select(sN, sP, NO3_c, NH4_c, PO4_NF, sAl, sCa, sFe, sK, sMg, rC, rN, rP, rAl, rCa, rFe, rK, rMg)
corr <- round(cor(CCNRalldata_corr, use = "pairwise.complete.obs"),3) 
p.mat <- cor_pmat(CCNRalldata_corr)
corr <- corr[c(1:10),c(11:18)]
p.mat <- p.mat[c(1:10),c(11:18)]
ggcorrplot(corr, p.mat = p.mat, insig = "blank", lab = TRUE,
           ggtheme = ggplot2::theme_classic, colors = c("#E46726", "white", "#6D9EC1"))

#soil/root nutrient correlogram 
ggpairs(CCNRalldata[,c("sN", "sP", "NO3_c", "NH4_c", "PO4_NF", "sAl", "sCa", "sFe", "sK", "sMg", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg")], columnLabels = c("sN", "sP", "NO3", "NH4", "PO4", "sAl", "sCa", "sFe", "sK", "sMg", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg"), aes(color=soildata$FT), upper = list(continuous = wrap('cor', size = 3)), lower = list(continuous = wrap("smooth")), diag = list(continuous = wrap("densityDiag", alpha = 0.3)), title = "Correlogram of Soil/Root Nutrients") + scale_fill_manual(values=met.brewer("Egypt", 3)) + scale_color_manual(values=met.brewer("Egypt", 3)) + theme_bw()


#### NMDS

#import data
tree_counts <- read_csv("/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Cleaned Data/Tree_Counts.csv") 
tree_counts <- tree_counts %>% tibble::column_to_rownames("Plot")

#match plot level averages
soilparameters <- soildata %>% dplyr::select("Plot", "Direction", "pH", "sC", "sN", "sP", "TOC_NF", "NO3_c", "NH4_c", "PO4_NF", "sAl", "sCa", "sFe", "sK", "sMg", "CEC", "TEB", "AlSat")
soilparameters <- soilparameters %>% dplyr::rename(C = sC, N = sN, P = sP, TOC = TOC_NF, NO3 = NO3_c, NH4 = NH4_c, PO4 = PO4_NF, Al = sAl, Ca = sCa, Fe = sFe, K = sK, Mg = sMg)
soilparameters_avg <- soilparameters %>% dplyr::select(-Direction)
soilparameters_avg <- soilparameters_avg %>% group_by(Plot) %>% summarise_all(mean, na.rm = TRUE)
soilparameters_avg <- soilparameters_avg %>% tibble::column_to_rownames("Plot")

#scree plot
dimcheckMDS(tree_counts, distance = "bray", k = 6, trymax = 100, autotransform = FALSE)
#NMDS
NMDS_WYK <- metaMDS(tree_counts, k = 3, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
stressplot(NMDS_WYK)
NMDS_WYK

#envfit
en = envfit(NMDS_WYK, soilparameters_avg, permutations = 999, na.rm = TRUE)
head(en)

#extract NMDS scores (x and y coordinates) for sites from vegan package
data.scores = as.data.frame(scores(NMDS_WYK)$sites)
data.scores <- data.scores %>% rownames_to_column("Plot") %>% full_join(ft)
data.scores$FT <- factor(data.scores$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))

#extract NMDS vectors for envfit variables
env.scores <- as.data.frame(scores(en, display = "vectors"))
env.scores <- env.scores %>% rownames_to_column("Variable")
env.scores <- cbind(env.scores, pval = en$vectors$pvals)

## Figure S1A
NMDS_plot = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
     geom_point(data = data.scores, aes(fill = FT), colour = "black", pch = 21, size = 5, alpha = 1) +
  stat_ellipse(aes(fill = FT), geom = "polygon", alpha = 0.2) +
     scale_fill_manual(values=met.brewer("Egypt", 3))  + 
  scale_color_manual(values=met.brewer("Egypt", 3))  +
     theme(axis.title = element_text(size = 20), 
       panel.background = element_blank(), panel.border = element_rect(fill = NA), 
       axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
       legend.title = element_text(size = 20), 
       legend.text = element_text(size =15)) + 
     labs(fill = "Forest Type")

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S1A", 
    width = 700, height = 600, units = "px")
NMDS_plot
dev.off()


#### PERMANOVA

adonis2(tree_counts ~ FT, data = ft)
anosim <- anosim(tree_counts, ft$FT, distance = "bray", permutations = 999)


#### PERMDISP

#create dissimilarity matrix
dis <- vegdist(tree_counts)
ft$FT <- factor(ft$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

#check for non-equal dispersions using betadisper() and permutest()
dispersion_model <- betadisper(dis, ft$FT)

permutest(dispersion_model, pairwise = TRUE) #include pairwise comparisons

#visualise differences in dispersions between groups
## Figure S1B
plot(dispersion_model, col = c("#dd5129", "#0f7ba2", "#43b284"), label.cex = 0.7)
disp <- recordPlot()

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S1B", 
    width = 450, height = 600, units = "px")
disp
dev.off()


## Hypothesis I: Soil nutrients

Despite being on the same parent material, secondary forests will have lower concentrations of plant-available nutrients than old-growth forests due to loss of nutrient pools in topsoil and plant biomass during past agricultural land use. 

#### Soil PCA

#select only relevant columns 
soilPCA <- soildata[c("sN", "sP", "NO3_c", "NH4_c", "PO4_NF", "sCa", "sFe", "sK", "sMg", "CEC", "TEB", "AlSat")] 

soilPCA <- soilPCA %>% dplyr::rename(N = sN, P = sP, NO3 = NO3_c, NH4 = NH4_c, PO4 = PO4_NF, Ca = sCa, Fe = sFe, K = sK, Mg = sMg)

#run PCA
res.pca.soil <- PCA(soilPCA, scale.unit = TRUE, ncp = 5, graph = FALSE)
print(res.pca.soil) # view PCA output

# extract eigenvalues
eig.val.soil <- get_eigenvalue(res.pca.soil)
# extract variables
var.soil <- get_pca_var(res.pca.soil)
# extract PCA individuals/observations
ind.soil <- get_pca_ind(res.pca.soil)

#scree plot
fviz_eig(res.pca.soil, addlabels = TRUE, ylim = c(0, 50), barfill = "#a3c585", barcolor = "black") + theme_classic() + ggtitle("Scree Plot (Soil PCA)") 
get_eigenvalue(res.pca.soil)

# visualise contribution of variables to PCs
var.soil$contrib
corrplot(var.soil$contrib, is.corr=FALSE, method = "circle", col = met.brewer("VanGogh3"), tl.col = "black", cl.align.text = "l" )

#PCA -- variables
fviz_pca_var(res.pca.soil, col.var = "contrib", 
             gradient.cols = met.brewer("Tam")) + 
  labs(colour = "Contribution (%)", x = "PC1 (44.7%)", y = "PC2 (14.1%)", 
       title = "Soil PCA (Variables Only)")

#PCA -- individuals
fviz_pca_ind(res.pca.soil, col.ind = soildata$FT, geom = "point") +
  labs(colour = "Forest Type", x = "PC1 (44.7%)", y = "PC2 (14.1%)", 
       title = "Soil PCA (Individual Plots)") +
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) + 
  scale_color_manual(values=met.brewer("Egypt", 3))

#PCA -- biplot
soilbiplot <- fviz_pca_biplot(res.pca.soil,
                col.ind = soildata$FT,
                col.var = "black", 
                geom = "point",
                addEllipses = TRUE, 
                #ellipse.level = 0.95,
                label = "var",
                labelsize = 6,
                repel = TRUE,
                mean.point = FALSE,
                legend.title = "Forest Type") + 
  labs(#colour = "Contribution (%)", 
       x = "PC1 (44.7%)", y = "PC2 (14.1%)", title = NULL) +
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  scale_fill_manual(values=met.brewer("Egypt", 3)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.position = "none") + theme_classic2() +
  theme(legend.title = element_text(size=18), legend.text = element_text(size=14))

## Figure 5A
soilbiplot

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_5A", 
    width = 750, height = 500, units = "px")
soilbiplot
dev.off()


#### LMMs of soil nutrients across FT

Linear mixed models will be run on all soil variables. Linear mixed models assume independence, linearity, normality, and equal variance. Posthoc tests using Tukey’s HSD method will be run to test for significant differences between forest types. 

#Total C
sCmodel <- lmer(log(sC) ~ FT + (1|Plot), data = soildata)
check_model(sCmodel)
anova(sCmodel, type = 2)

#Total N
sNmodel <- lmer(log(sN) ~ FT + (1|Plot), data = soildata)
check_model(sNmodel)
anova(sNmodel, type = 2)

#Total P
sPmodel <- lmer(log(sP) ~ FT + (1|Plot), data = soildata)
check_model(sPmodel)
anova(sPmodel, type = 2)

#NO3
NO3model <- lmer(log(NO3_c) ~ FT + (1|Plot), data = soildata)
check_model(NO3model)
anova(NO3model, type = 2)

#NH4
NH4model <- lmer(log(NH4_c) ~ FT + (1|Plot), data = soildata)
check_model(NH4model)
anova(NH4model, type = 2)

#PO4
PO4model <- lmer(sqrt(PO4_NF) ~ FT + (1|Plot), data = soildata)
check_model(PO4model)
anova(PO4model, type = 2)

#Al
sAlmodel <- lmer((sAl) ~ FT + (1|Plot), data = soildata)
check_model(sAlmodel)
anova(sAlmodel, type = 2)

#Ca
sCamodel <- lmer(log(sCa) ~ FT + (1|Plot), data = soildata)
check_model(sCamodel)
anova(sCamodel, type = 2)

#Fe
sFemodel <- lmer(log(sFe) ~ FT + (1|Plot), data = soildata)
check_model(sFemodel)
anova(sFemodel, type = 2)

#K
sKmodel <- lmer(log(sK) ~ FT + (1|Plot), data = soildata)
check_model(sKmodel)
anova(sKmodel, type = 2)

#Mg
sMgmodel <- lmer(log(sMg) ~ FT + (1|Plot), data = soildata)
check_model(sMgmodel)
anova(sMgmodel, type = 2)

#pH
pHmodel <- lmer((pH) ~ FT + (1|Plot), data = soildata)
check_model(pHmodel)
anova(pHmodel, type = 2)

#CEC
CECmodel <- lmer(log(CEC) ~ FT + (1|Plot), data = soildata)
check_model(CECmodel)
anova(CECmodel, type = 2)

#TEB
TEBmodel <- lmer((TEB) ~ FT + (1|Plot), data = soildata)
check_model(TEBmodel)
anova(TEBmodel, type = 2)

#AlSat
AlSatmodel <- lmer((AlSat) ~ FT + (1|Plot), data = soildata)
check_model(AlSatmodel)
anova(AlSatmodel, type = 2)


#### Nutrient ratios

#soil ratios (totals, molar ratio)
histogram(soildata$sC_N)
histogram(soildata$sN_P)
histogram(soildata$sC_P)
# soil NO3:NH4 (SEAL)
soildata$NO3_NH4 <- soildata$NO3_c / soildata$NH4_c
histogram(soildata$NO3_NH4)

#root ratios (totals, molar ratio)
histogram(rootdata$rC_N)
histogram(rootdata$rN_P)
histogram(rootdata$rC_P)

#### LMMs ####

# soil C:N
sCNmodel <- lmer((sC_N) ~ FT + (1|Plot), data = soildata)
check_model(sCNmodel)
anova(sCNmodel, type = 2)
summary(sCNmodel)
soildata$FT <- relevel(soildata$FT, ref = "Old-growth")
sCNmodel <- lmer((sC_N) ~ FT + (1|Plot), data = soildata)
anova(sCNmodel, type = 2)
summary(sCNmodel)
library(emmeans)
emmeans(sCNmodel, list(pairwise~FT), adjust = "tukey")
3.53/(10.05 + 3.53) #explains 26% of variation left over

# soil C:P **SIGNIFICANT b/w FT 2/4**
sCPmodel <- lmer((sC_P) ~ FT + (1|Plot), data = soildata)
check_model(sCPmodel)
anova(sCPmodel, type = 2)
summary(sCPmodel)
soildata$FT <- relevel(soildata$FT, ref = "Old-growth")
sCPmodel <- lmer((sC_P) ~ FT + (1|Plot), data = soildata)
anova(sCPmodel, type = 2)
summary(sCPmodel)
346734/(346734 + 469595) #explains 42.5% of variation left over

# soil N:P **SIGNIFICANT b/w FT 2/4**
sNPmodel <- lmer((sN_P) ~ FT + (1|Plot), data = soildata)
check_model(sNPmodel)
anova(sNPmodel, type = 2)
summary(sNPmodel)
soildata$FT <- relevel(soildata$FT, ref = "Early secondary")
sNPmodel <- lmer((sN_P) ~ FT + (1|Plot), data = soildata)
anova(sNPmodel, type = 2)
summary(sNPmodel)
438.8/(438.8 + 738.1) #explains 37.3% of variation left over

# soil NO3:NH4
sNratiosmodel <- lmer(sqrt(NO3_NH4) ~ FT + (1|Plot), data = soildata)
check_model(sNratiosmodel)
anova(sNratiosmodel, type = 2)

# root N:P
rNPmodel <- lmer(log(rN_P) ~ FT + (1|Plot), data = rootdata)
check_model(rNPmodel)
anova(rNPmodel, type = 2)

# root C:P
rCPmodel <- lmer(log(rC_P) ~ FT + (1|Plot), data = rootdata)
check_model(rCPmodel)
anova(rCPmodel, type = 2)

# root C:N
rCNmodel <- lmer(log(rC_N) ~ FT + (1|Plot), data = rootdata)
check_model(rCNmodel)
anova(rCNmodel, type = 2)

#### EMMs ####

#soil C:N
sCN_emmeans <- emmeans(sCNmodel, ~FT)
sCN_emmeans <- as.data.frame(sCN_emmeans)
#order forest types in order
sCN_emmeans$FT <- factor(sCN_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
#plot
sCN_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = sCN_emmeans) +
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = sCN_emmeans) + 
  geom_jitter(aes(x = FT, y = sC_N, color = FT), width = 0.25, alpha = 0.5, 
             #position = position_nudge(x = -0.2, y = 0), 
             data = soildata) + 
  labs(y = "Soil C:N") +
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#soil C:P
sCP_emmeans <- emmeans(sCPmodel, ~FT)
sCP_emmeans <- as.data.frame(sCP_emmeans)
#order forest types in order
sCP_emmeans$FT <- factor(sCP_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(sCPmodel, list(pairwise~FT), adjust = "tukey")
#plot
sCP_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = sCP_emmeans) +
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = sCP_emmeans) + 
  geom_jitter(aes(x = FT, y = sC_P, color = FT), width = 0.25, alpha = 0.5, 
             #position = position_nudge(x = -0.2, y = 0), 
             data = soildata) + 
  labs(y = "Soil C:P", x = "Forest Type") +
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.title.x = element_text(vjust=-0.5))

#soil N:P

sNP_emmeans <- emmeans(sNPmodel, ~FT)
sNP_emmeans <- as.data.frame(sNP_emmeans)
#order forest types in order
sNP_emmeans$FT <- factor(sNP_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(sNPmodel, list(pairwise~FT), adjust = "tukey")
#plot
sNP_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = sNP_emmeans) +
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = sNP_emmeans) + 
  geom_jitter(aes(x = FT, y = sN_P, color = FT), width = 0.25, alpha = 0.5, 
             #position = position_nudge(x = -0.2, y = 0), 
             data = soildata) + 
  labs(y = "Soil N:P") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#plot

nutrientratioCIs <- (sCN_plot / sNP_plot / sCP_plot)

## Figure 2

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_2", 
    width = 450, height = 600, units = "px")
nutrientratioCIs + plot_annotation(tag_levels = "A")
dev.off()


## Hypothesis II: Root traits

Differences in tree species composition due to land-use change corresponds to differences in root morphological traits between forest types. Early-successional secondary forests are dominated by pioneer species and will thus have more acquisitive root morphological traits, relying on rapid turnover and soil exploration as a nutrient acquisition strategy. On the other hand, trees in old-growth forest plots, some of which are ECM-associated, will have more conservative root traits and increased “collaboration” with mycorrhizal partners, relying on established mycorrhizal networks for nutrient mining advantages.

#### Root PCA

#select only relevant columns 
rootPCA <- rootdata[c("rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg", "PME", "TotalRootBiomass", "SRL", "SRA", "D", "RTD", "Tips")]

rootPCA <- rootPCA %>% dplyr::rename(C = rC, N = rN, P = rP, Al = rAl, Ca = rCa, Fe = rFe, K = rK, Mg = rMg)

#run PCA
res.pca.root <- PCA(rootPCA, scale.unit = TRUE, ncp = 5, graph = FALSE)
print(res.pca.root) # view PCA output

#invert PCA axes
#res.pca.root$var$coord <- res.pca.root$var$coord*-1
#res.pca.root$ind$coord <- res.pca.root$ind$coord*-1

# extract eigenvalues
eig.val.root <- get_eigenvalue(res.pca.root)
# extract variables
var.root <- get_pca_var(res.pca.root)
# extract PCA individuals/observations
ind.root <- get_pca_ind(res.pca.root)

#scree plot
fviz_eig(res.pca.root, addlabels = TRUE, ylim = c(0, 50), barfill = "#a3c585", barcolor = "black") + theme_classic() + ggtitle("Scree Plot (Root PCA)") 
get_eigenvalue(res.pca.root)

# visualise contribution of variables to PCs
var.root$contrib
corrplot(var.root$contrib, is.corr=FALSE, method = "circle", col = met.brewer("VanGogh3"), tl.col = "black", cl.align.text = "l" )

#PCA -- variables (PC1 and PC2)
fviz_pca_var(res.pca.root, col.var = "contrib", 
             gradient.cols = met.brewer("Tam")) + 
  labs(colour = "Contribution (%)", x = "PC1 (18.7%)", y = "PC2 (16.1%)", 
       title = "Root PCA (Variables Only)") + scale_y_reverse()

##PC 1 vs PC 3
fviz_pca_var(res.pca.root, axes = c(1, 3), col.var = "contrib", 
             gradient.cols = met.brewer("Tam")) +
  labs(colour = "Contribution (%)", x = "PC1 (18.7%)", y = "PC3 (12.8%)", 
       title = "Root PCA (Variables Only)")+ scale_y_reverse()

##PC 2 vs PC 3
fviz_pca_var(res.pca.root, axes = c(2, 3), col.var = "contrib", 
             gradient.cols = met.brewer("Tam")) +
  labs(colour = "Contribution (%)", x = "PC2 (16.1%)", y = "PC3 (12.8%)", 
       title = "Root PCA (Variables Only)")

#PCA -- individuals
fviz_pca_ind(res.pca.root, col.ind = rootdata$FT, geom = "point",
             pointsize = rootdata$BA_ECM) +
  labs(colour = "Forest Type", x = "PC1 (18.7%)", y = "PC2 (16.1%)", 
       title = "Root PCA (Individual Plots)") +
   labs(size = "ECM Basal Area (%)") +
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) + 
  scale_color_manual(values=met.brewer("Egypt", 3))+ scale_y_reverse()

rootdata$FT <- factor(rootdata$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))

#PCA -- biplot PC1 vs PC2
rootbiplot12 <- fviz_pca_biplot(res.pca.root,
                col.ind = rootdata$FT,
                col.var = "black",
                geom = "point",
                pointsize = rootdata$BA_ECM,
                addEllipses = TRUE, 
                label = "var",
                labelsize = 6,
                repel = TRUE,
                mean.point = FALSE,
                legend.title = "Forest Type") + 
  labs(x = "PC1 (18.7%)", y = "PC2 (16.1%)", title = NULL) +
  labs(size = "ECM Basal Area (%)") +
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  scale_fill_manual(values=met.brewer("Egypt", 3)) +
  #theme(legend.position = "none") + 
  theme_classic2() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=15))+ scale_y_reverse()

## Figure 3A
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_3A", 
    width = 700, height = 600, units = "px")
rootbiplot12
dev.off()


#### LMMs for root morphological traits

rootdata$FT <- relevel(rootdata$FT, ref = "Old-growth")

#SRL -- significant
SRLmodel <- lmer(log(SRL) ~ FT + (1|Plot), data = rootdata)
check_model(SRLmodel) #assumptions met
anova(SRLmodel, type = 2)
summary(SRLmodel)
rootdata$FT <- relevel(rootdata$FT, ref = "Early secondary")
rootdata$FT <- relevel(rootdata$FT, ref = "Old-growth")

## EMMs for SRL
SRL_emmeans <- emmeans(SRLmodel, ~FT)
SRL_emmeans <- as.data.frame(SRL_emmeans)
#back transform
SRL_emmeans$emmean <- exp(SRL_emmeans$emmean)
SRL_emmeans$lower.CL <- exp(SRL_emmeans$lower.CL)
SRL_emmeans$upper.CL <- exp(SRL_emmeans$upper.CL)
SRL_emmeans

#order forest types in order
SRL_emmeans$FT <- factor(SRL_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(SRLmodel, list(pairwise~FT), adjust = "tukey")
#plot
SRL_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = SRL_emmeans) + scale_y_continuous(trans='log10') + geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = SRL_emmeans) + 
  geom_jitter(aes(x = FT, y = (SRL), color = FT), alpha = 0.5, width = 0.2,
             #position = position_nudge(x = -0.2, y = 0), 
             data = rootdata) + 
  labs(x = "Forest Type", y = "Specific Root Length (m/g)") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=10)) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#SRA -- significant
SRAmodel <- lmer(log(SRA) ~ FT + (1|Plot), data = rootdata)
check_model(SRAmodel) #assumptions met
anova(SRAmodel, type = 2)
summary(SRAmodel)
rootdata$FT <- relevel(rootdata$FT, ref = "Early secondary")
rootdata$FT <- relevel(rootdata$FT, ref = "Old-growth")

## EMMs for SRA
SRA_emmeans <- emmeans(SRAmodel, ~FT)
SRA_emmeans <- as.data.frame(SRA_emmeans)
#back transform
SRA_emmeans$emmean <- exp(SRA_emmeans$emmean)
SRA_emmeans$lower.CL <- exp(SRA_emmeans$lower.CL)
SRA_emmeans$upper.CL <- exp(SRA_emmeans$upper.CL)
SRA_emmeans

#order forest types in order
SRA_emmeans$FT <- factor(SRA_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(SRAmodel, list(pairwise~FT), adjust = "tukey")
#plot
SRA_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = SRA_emmeans) + scale_y_continuous(trans='log10') + 
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = SRA_emmeans) + 
  geom_jitter(aes(x = FT, y = (SRA), color = FT), alpha = 0.5, width = 0.2,
             #position = position_nudge(x = -0.2, y = 0), 
             data = rootdata) + 
  labs(x = "Forest Type", y = "Specific Root Area (cm^2)") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=10)) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#D -- significant
Dmodel <- lmer(log(D) ~ FT + (1|Plot), data = rootdata)
check_model(Dmodel) #assumptions met
anova(Dmodel, type = 2)
summary(Dmodel)
rootdata$FT <- relevel(rootdata$FT, ref = "Early secondary")
rootdata$FT <- relevel(rootdata$FT, ref = "Old-growth")

## EMMs for D
D_emmeans <- emmeans(Dmodel, ~FT)
D_emmeans <- as.data.frame(D_emmeans)
#back transform
D_emmeans$emmean <- exp(D_emmeans$emmean)
D_emmeans$lower.CL <- exp(D_emmeans$lower.CL)
D_emmeans$upper.CL <- exp(D_emmeans$upper.CL)
D_emmeans

#order forest types in order
D_emmeans$FT <- factor(D_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(Dmodel, list(pairwise~FT), adjust = "tukey")
#plot
D_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = D_emmeans) + scale_y_continuous(trans='log10') + 
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = D_emmeans) + 
  geom_jitter(aes(x = FT, y = (D), color = FT), alpha = 0.5, width = 0.2,
             #position = position_nudge(x = -0.2, y = 0), 
             data = rootdata) + 
  labs(x = "Forest Type", y = "Root Diameter (mm)") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=10))

#RTD -- significant
RTDmodel <- lmer(log(RTD) ~ FT + (1|Plot), data = rootdata)
check_model(RTDmodel) #assumptions met
anova(RTDmodel, type = 2)
summary(RTDmodel)
rootdata$FT <- relevel(rootdata$FT, ref = "Early secondary")
rootdata$FT <- relevel(rootdata$FT, ref = "Old-growth")

## EMMs for RTD
RTD_emmeans <- emmeans(RTDmodel, ~FT)
RTD_emmeans <- as.data.frame(RTD_emmeans)
#back transform
RTD_emmeans$emmean <- exp(RTD_emmeans$emmean)
RTD_emmeans$lower.CL <- exp(RTD_emmeans$lower.CL)
RTD_emmeans$upper.CL <- exp(RTD_emmeans$upper.CL)
RTD_emmeans

#order forest types in order
RTD_emmeans$FT <- factor(RTD_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
emmeans(RTDmodel, list(pairwise~FT), adjust = "tukey")
#plot
RTD_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = RTD_emmeans) + scale_y_continuous(trans='log10') + 
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = RTD_emmeans) + 
  geom_jitter(aes(x = FT, y = (RTD), color = FT), alpha = 0.5, width = 0.2,
             #position = position_nudge(x = -0.2, y = 0), 
             data = rootdata) + 
  labs(x = "Forest Type", y = "Root Tissue Density (g/cm^3)") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=10))

#Tips
Tipsmodel <- lmer(log(Tips) ~ FT + (1|Plot), data = rootdata)
check_model(Tipsmodel) #assumptions met
anova(Tipsmodel, type = 2)

## Figure 4
roottraitsCIs <- (SRL_plot | SRA_plot ) / (D_plot | RTD_plot)

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_4", 
    width = 850, height = 500, units = "px")
roottraitsCIs + plot_annotation(tag_levels = "A")
dev.off()


#### LMMs for other root traits

#root enzyme activity -- not significantly different
enzymemodel <- lmer(log(PME) ~ FT + (1|Plot), data = rootdata)
check_model(enzymemodel)
anova(enzymemodel, type = 2)

#root biomass
rootbiomassmodel <- lmer((TotalRootBiomass) ~ FT + (1|Plot), data = CCNRalldata)
check_model(rootbiomassmodel)
anova(rootbiomassmodel, type = 2)


## Hypothesis III: Root variability

#### Extract all PCs

#extract coordinates of soil PCA
soilcoord <- data.frame(res.pca.soil$ind$coord)
soilcoord$row_num <- seq.int(nrow(soilcoord)) #add column to join
#join by row number
soildatadf <- soildata
soildatadf$row_num <- seq.int(nrow(soildatadf)) 

soilcoord <- soildatadf %>% right_join(soilcoord) %>%
  dplyr::select(row_num, Plot, Direction, FT, Dim.1, Dim.2, BA_ECM) %>% dplyr::rename(soilDim.1 = Dim.1) %>% dplyr::rename(soilDim.2 = Dim.2)

#extract coordinates of root PCA
rootcoord <- data.frame(res.pca.root$ind$coord)
rootcoord$Dim.2 <- rootcoord$Dim.2*-1
rootcoord$row_num <- seq.int(nrow(rootcoord)) #add column to join
#join to soil PCA by row number
soilrootcoord <- soilcoord %>% right_join(rootcoord) 
soilrootcoord <- soilrootcoord %>%
  dplyr::select(Plot, Direction, FT, soilDim.1, soilDim.2, Dim.1, Dim.2, Dim.3, BA_ECM) %>% dplyr::rename(rootDim.1 = Dim.1, rootDim.2 = Dim.2, rootDim.3 = Dim.3)

#change FT reference group to FT4
soilrootcoord$FT <- relevel(soilrootcoord$FT, ref = "Old-growth")


#### LMMs for soil PCs vs FT and root PCs against soil and FT

##### Soil PCs

#Soil PC1
soil1 <- lmer((soilDim.1) ~ FT + (1|Plot), data = soilrootcoord)
check_model(soil1)
anova(soil1, type = 2)
summary(soil1)

#Soil PC2
soil2 <- lmer((soilDim.2) ~ FT + (1|Plot), data = soilrootcoord)
check_model(soil2)
anova(soil2, type = 2)
summary(soil2)


##### Root PC1 (Collaboration gradient)

#against Soil PC1
root1soil1 <- lmer((rootDim.1) ~ soilDim.1+FT + (1|Plot), data = soilrootcoord)
check_model(root1soil1)
anova(root1soil1, type = 2)
summary(root1soil1)
emmeans(root1soil1, list(pairwise~FT), adjust = "tukey")

root1soil1ECM <- lmer((rootDim.1) ~ soilDim.1*BA_ECM + (1|Plot), data = soilrootcoord)
check_model(root1soil1ECM)
anova(root1soil1ECM, type = 2)

#against Soil PC2
root1soil2 <- lmer((rootDim.1) ~ soilDim.2*FT + (1|Plot), data = soilrootcoord)
check_model(root1soil2)
anova(root1soil2, type = 2)

root1soil2ECM <- lmer((rootDim.1) ~ soilDim.2*BA_ECM + (1|Plot), data = soilrootcoord)
check_model(root1soil2ECM)
anova(root1soil2ECM, type = 2)


#soil PC1 vs root PC1
soilrootcoord$FT <- factor(soilrootcoord$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

#plot EMMs
rootcollab_FT_emmeans <- emmeans(root1soil1, ~FT)
rootcollab_FT_emmeans <- as.data.frame(rootcollab_FT_emmeans)
#order forest types in order
rootcollab_FT_emmeans$FT <- factor(rootcollab_FT_emmeans$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))

## Figure 3B
rootcollabFT_plot <- ggplot() + geom_point(aes(x = FT, y = emmean, color = FT), 
  position = position_dodge(width = 0.5), size = 5, data = rootcollab_FT_emmeans) +
  geom_errorbar(aes(x = FT, ymin = lower.CL, ymax = upper.CL, color = FT),
                width = 0.2, size = 1, 
                position = position_dodge(width = 0.5), 
                data = rootcollab_FT_emmeans) + 
  geom_jitter(aes(x = FT, y = rootDim.1, color = FT), alpha = 0.5, width = 0.2,
            #position = position_nudge(x = -0.2, y = 0), 
            data = soilrootcoord) + 
  labs(y = "Predicted Root Collaboration Gradient", x = "Forest Type") + 
  scale_color_manual(values=met.brewer("Egypt", 3)) + 
  theme_classic2() + theme(legend.position = "none") + 
  coord_flip() +
  theme(axis.text.y = element_text(angle=90, hjust = 0.5)) + 
  theme(axis.title.x =element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x =element_blank())

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_3B", 
    width = 450, height = 600, units = "px")
rootcollabFT_plot
dev.off()


##### Root PC2 (Conservation gradient)

soilrootcoord$FT <- relevel(soilrootcoord$FT, ref = "Old-growth")

#against Soil PC1
root2soil1 <- lmer((rootDim.2) ~ soilDim.1+FT + (1|Plot), data = soilrootcoord)
check_model(root2soil1)
anova(root2soil1, type = 2)
summary(root2soil1)

root2soil1ECM <- lmer((rootDim.2) ~ soilDim.1*BA_ECM + (1|Plot), data = soilrootcoord)
check_model(root2soil1ECM)
anova(root2soil1ECM, type = 2)

#against Soil PC2
root2soil2 <- lmer((rootDim.2) ~ soilDim.2*FT + (1|Plot), data = soilrootcoord)
check_model(root2soil2)
anova(root2soil2, type = 2)

root2soil2ECM <- lmer((rootDim.2) ~ soilDim.2*BA_ECM + (1|Plot), data = soilrootcoord)
check_model(root2soil2ECM)
anova(root2soil2ECM, type = 2)


soilrootcoord$FT <- factor(soilrootcoord$FT, level = c("Early secondary", "Mature secondary", "Old-growth"))

#raw results
soil1root2pc <- soilrootcoord %>% ggplot(aes(x = soilDim.1, y = rootDim.2)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm", alpha = .2) + ylim(-5, 5) +
  labs(x = "Soil Fertility Gradient", y = "Root Conservation Gradient") + theme(axis.title=element_text(size=14)) + 
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) +
     guides(size = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0)))) + 
  theme(axis.text = element_text(size=15), 
        axis.title=element_text(size=18)) + theme(legend.position = "none")
soil1root2pc

#EMMs
predictEMM <- predict_response(root2soil1, terms = "soilDim.1")

## Figure 5B
soil1root2pc <- ggplot(predictEMM, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Predicted Root Conservation Gradient", x = "Soil Fertility Gradient") + theme_classic2() +
  geom_point(aes(x = soilDim.1, y = rootDim.2), alpha = 0.5, 
             #position = position_nudge(x = -0.2, y = 0), 
             data = soilrootcoord) + 
    theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_5B", 
    width = 600, height = 350, units = "px")
soil1root2pc
dev.off()


#### Soil gradient vs root traits

#combine soil PC1 with root data
soilgradient <- soilcoord %>% dplyr::select(Plot, Direction, soilDim.1, soilDim.2)
rootsoilgradientdf <- rootdata %>% right_join(soilgradient) %>% dplyr::select(-PercentNecro)
mydata <- rootsoilgradientdf[,c(4:18, 24, 25)] #check columns

#correlation matrix 
mydata.cor <- cor(mydata, use = "pairwise.complete.obs") #ignoring NAs
round(mydata.cor, 2) #round to two decimals

#correlation matrix with p values
mydata.rcorr <- rcorr(as.matrix(mydata))
mydata.coeff <- mydata.rcorr$r #extract the correlation coefficients
mydata.p <- mydata.rcorr$P #extract p-values
mydata.p[is.na(mydata.p)] <- 1 #replace NA with 1 in p values

#visualise with corrplot
corrplot(mydata.cor, type = "upper",
         tl.col = "black", tl.srt = 45)

# Insignificant correlation are crossed
par(bg = "#FFFFFF")
corrplot(mydata.cor, type="lower", method = "color", addgrid.col = 'black', tl.col = "black",
         p.mat = mydata.p, sig.level = 0.05, insig = "blank")

#correlogram like above
ggpairs(rootsoilgradientdf[,c("PME", "TotalRootBiomass", "SRL", "SRA", "D", "RTD", "Tips", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg", "soilDim.1")],
        columnLabels = c("PME", "Biomass", "SRL", "SRA", "D", "RTD", "# Tips", "C", "N", "P", "Al", "Ca", "Fe", "K", "Mg", "Soil.Gradient"),
        upper = list(continuous = wrap('cor', size = 3)),
        lower = list(continuous = wrap("smooth", size = 0.1, alpha = 0.3)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.3)),
        title = "Correlogram of Root Traits with Soil Fertility") + theme_bw()


#### LMMs: Root nutrients vs. each soil nutrient

#rC vs. sC
rCsoilmodel <- lmer((rC) ~ sC*FT + (1|Plot), data = CCNRalldata)
check_model(rCsoilmodel)
anova(rCsoilmodel, type = 2)

#rC vs. soil TOC -- significant
rCsoilTOCmodel <- lmer((TOC_NF) ~ rC*FT + (1|Plot), data = CCNRalldata)
check_model(rCsoilTOCmodel)
anova(rCsoilTOCmodel, type = 2)
#plot
predictrCEMM <- predict_response(rCsoilTOCmodel, terms = "rC")
rC_plot <- ggplot(predictrCEMM, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Predicted Total Organic Carbon \n (mg/kg soil)", x = "Total Root Carbon (ug/g)") + theme_classic2() +
  geom_point(aes(x = rC, y = (TOC_NF)), alpha = 0.5, 
             data = CCNRalldata)

#total root biomass vs. soil TOC -- significant
biomasssoilTOCmodel <- lmer((TOC_NF) ~ TotalRootBiomass*FT + (1|Plot), data = CCNRalldata)
check_model(biomasssoilTOCmodel)
anova(biomasssoilTOCmodel, type = 2)
#plot
predictbiomassEMM <- predict_response(biomasssoilTOCmodel, terms = "TotalRootBiomass")
biomassTOC_plot <- ggplot(predictbiomassEMM, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Predicted Total Organic Carbon \n (mg/kg soil)", x = "Total Root Biomass (g)") + theme_classic2() +
  geom_point(aes(x = TotalRootBiomass, y = (TOC_NF)), alpha = 0.5, 
             data = CCNRalldata)

#rN vs sN
rNsoilmodel <- lmer((rN) ~ sN*FT + (1|Plot), data = CCNRalldata)
check_model(rNsoilmodel)
anova(rNsoilmodel, type = 2)

#rN vs NO3
rNsoilNO3model <- lmer((rN) ~ NO3_c*FT + (1|Plot), data = CCNRalldata)
check_model(rNsoilNO3model)
anova(rNsoilNO3model, type = 2)

#rN vs NH4
rNsoilNH4model <- lmer((rN) ~ NH4_c*FT + (1|Plot), data = CCNRalldata)
check_model(rNsoilNH4model)
anova(rNsoilNH4model, type = 2)

#rP vs sP
rPsoilmodel <- lmer((rP) ~ sP*FT + (1|Plot), data = CCNRalldata)
check_model(rPsoilmodel)
anova(rPsoilmodel, type = 2)
emmeans(rPsoilmodel, list(pairwise~FT), adjust = "tukey")
#not significant when FT alone
rPsoilmodel <- lmer((rP) ~ FT + (1|Plot), data = CCNRalldata)
check_model(rPsoilmodel)
anova(rPsoilmodel, type = 2)

#rP vs soil PO4
rPsoilPO4model <- lmer((rP) ~ PO4_NF*FT + (1|Plot), data = CCNRalldata)
check_model(rPsoilPO4model)
anova(rPsoilPO4model, type = 2)

#rAl vs sAl -- significant
rAlsoilmodel <- lmer((rAl) ~ AlSat*FT + (1|Plot), data = CCNRalldata)
check_model(rAlsoilmodel)
anova(rAlsoilmodel, type = 2)
#plot
predictAlEMM <- predict_response(rAlsoilmodel, terms = "AlSat")
rAlsoil_plot <- ggplot(predictAlEMM, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Predicted Root Al (ug/g)", x = "Soil Al Saturation (% CEC)") + theme_bw() +
  theme(axis.text = element_text(size=15), 
        axis.title=element_text(size=18))

#rCa vs sCa
rCasoilmodel <- lmer((rCa) ~ sCa*FT + (1|Plot), data = CCNRalldata)
check_model(rCasoilmodel)
anova(rCasoilmodel, type = 2)

#rFe vs sFe
rFesoilmodel <- lmer((rFe) ~ sFe*FT + (1|Plot), data = CCNRalldata)
check_model(rFesoilmodel)
anova(rFesoilmodel, type = 2)

#rK vs sK
rKsoilmodel <- lmer((rK) ~ sK*FT + (1|Plot), data = CCNRalldata)
check_model(rKsoilmodel)
anova(rKsoilmodel, type = 2)

#rMg vs sMg
rMgsoilmodel <- lmer(log(rMg) ~ sMg*FT + (1|Plot), data = CCNRalldata)
check_model(rMgsoilmodel)
anova(rMgsoilmodel, type = 2)

#rCa vs. CEC
rCaCECsoilmodel <- lmer((rCa) ~ CEC*FT + (1|Plot), data = CCNRalldata)
check_model(rCaCECsoilmodel)
anova(rCaCECsoilmodel, type = 2)

#rFe vs. CEC
rFeCECsoilmodel <- lmer((rFe) ~ CEC*FT + (1|Plot), data = CCNRalldata)
check_model(rFeCECsoilmodel)
anova(rFeCECsoilmodel, type = 2)
root_FeCEC_plot <- CCNRalldata %>% ggplot(aes(CEC, rFe)) + geom_point() + geom_smooth(method = lm) +
  xlab("CEC") + ylab("Root Fe") + 
  scale_shape_manual(name = "Forest Type", values = c(19, 19, 19)) +
  scale_color_manual(values=met.brewer("Egypt", 1)) + 
  scale_fill_manual(values=met.brewer("Egypt", 1)) + theme_classic() + theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))

#rK vs. CEC
rKCECsoilmodel <- lmer((rK) ~ CEC*FT + (1|Plot), data = CCNRalldata)
check_model(rKCECsoilmodel)
anova(rKCECsoilmodel, type = 2)

#rMg vs. CEC
rMgCECsoilmodel <- lmer((rMg) ~ CEC*FT + (1|Plot), data = CCNRalldata)
check_model(rMgCECsoilmodel)
anova(rMgCECsoilmodel, type = 2)

rootnutrientplots <- (rC_plot | biomassTOC_plot)
rootnutrientplots + plot_annotation(tag_levels = "A")


#### Plotting Soil fertility x FT vs each root trait

#root PME
PMEmodel <- lmer((PME) ~ soilDim.1 * FT + (1|Plot), data = rootsoilgradientdf)
check_model(PMEmodel)
anova(PMEmodel, type = 2)

#root biomass
biomassmodel <- lmer((TotalRootBiomass) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(biomassmodel)
anova(biomassmodel, type = 2)
#**significant for soil fertility

#SRL
SRLmodel <- lmer((SRL) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(SRLmodel)
anova(SRLmodel, type = 2)
#**significant for FT

#SRA
SRAmodel <- lmer((SRA) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(SRAmodel)
anova(SRAmodel, type = 2)
#**significant for FT

#D
Dmodel <- lmer((D) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(Dmodel)
anova(Dmodel, type = 2)
#**significant for FT

#RTD
RTDmodel <- lmer((RTD) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(RTDmodel)
anova(RTDmodel, type = 2)
#**significant for FT

#Tips
tipsmodel <- lmer((Tips) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(tipsmodel)
anova(tipsmodel, type = 2)

#root N
rNmodel <- lmer((rN) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rNmodel)
anova(rNmodel, type = 2)

#root C
rCmodel <- lmer((rC) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rCmodel)
anova(rCmodel, type = 2)

#root P
rPmodel <- lmer((rP) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rPmodel)
anova(rPmodel, type = 2)
#**significant for soil fertility

#root Al
rAlmodel <- lmer((rAl) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rAlmodel)
anova(rAlmodel, type = 2)
#**significant for soil fertility

#root Ca
rCamodel <- lmer((rCa) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rCamodel)
anova(rCamodel, type = 2)
#**significant for soil fertility

#root Fe
rFemodel <- lmer((rFe) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rFemodel)
anova(rFemodel, type = 2)
#**significant for soil fertility

#root K
rKmodel <- lmer((rK) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rKmodel)
anova(rKmodel, type = 2)
#**significant for FT
emmeans(rKmodel, list(pairwise~FT), adjust = "tukey")

#root Mg
rMgmodel <- lmer((rMg) ~ soilDim.1 + FT + (1|Plot), data = rootsoilgradientdf)
check_model(rMgmodel)
anova(rMgmodel, type = 2)

#predict biomass
predictbiomass <- predict_response(biomassmodel, terms = "soilDim.1")
predictbiomass$group <- gsub('1','TotalRootBiomass',predictbiomass$group)
predictbiomass$group <- as.factor(predictbiomass$group)

biomass_plot <- ggplot(predictbiomass, aes(x = x, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  labs(y = "Predicted Total Root Biomass (g)", x = "Soil Fertility Gradient") + theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#plot
predictrP <- predict_response(rPmodel, terms = "soilDim.1")
predictrP$group <- gsub('1','P',predictrP$group)
predictrCa <- predict_response(rCamodel, terms = "soilDim.1")
predictrCa$group <- gsub('1','Ca',predictrCa$group)
predictrAl <- predict_response(rAlmodel, terms = "soilDim.1")
predictrAl$group <- gsub('1','Al',predictrAl$group)
predictrFe <- predict_response(rFemodel, terms = "soilDim.1")
predictrFe$group <- gsub('1','Fe',predictrFe$group)

rootcation_predict <- rbind(predictrCa, predictrAl, predictrFe)
rootcation_predict$group <- as.factor(rootcation_predict$group)

ggplot(rootcation_predict, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  labs(y = "Predicted Root Cation (ug/g)", x = "Soil Fertility Gradient") + theme_bw() +
  theme(axis.text = element_text(size=15), 
        axis.title=element_text(size=18))

rootcation_predict <- rbind(predictrP, predictrCa, predictrAl, predictrFe)
rootcation_predict$group <- factor(rootcation_predict$group, levels=c("Al", "Ca", "Fe", "P"))

rootcation_predicted <- ggplot(rootcation_predict, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) + ylim(0.25, 14) +
  scale_y_break(c(0.5, 1.2), scales = 4) + 
  labs(y = "Predicted Root Content (ug/g)", x = "Soil Fertility Gradient", color = "Chemical Trait", fill = "Chemical Trait") + theme_bw() +
  scale_color_manual(values=met.brewer("Austria", 5)) +
  scale_fill_manual(values=met.brewer("Austria", 5)) +
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
#+ theme(legend.position = "none")

## Figure 5C
rootgradients <- (biomass_plot / rootcation_predicted)

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_5C", 
    width = 500, height = 700, units = "px")
rootgradients
dev.off()

#root vs soil nutrients
rootsoilgradient_corr <- rootsoilgradientdf %>% dplyr::select(soilDim.1, TotalRootBiomass, rP, rAl, rCa, rFe)
corr <- round(cor(rootsoilgradient_corr, use = "pairwise.complete.obs"),3) 
p.mat <- cor_pmat(rootsoilgradient_corr)
ggcorrplot(corr, p.mat = p.mat, lab = TRUE,
           ggtheme = ggplot2::theme_classic, colors = c("#E46726", "white", "#6D9EC1"))


#### RDA of soil var vs root traits

##### Sort data

# sort soil data
soilparameters_RDA <- soilgradient
soilparameters_RDA$Sample = paste(soilparameters_RDA$Plot, soilparameters_RDA$Direction, sep="-")
soilparameters_RDA <- soilparameters_RDA %>% dplyr::select(-Plot, -Direction) %>% dplyr::select(Sample, everything())
soilparameters_RDA <- na.omit(soilparameters_RDA)

# sort root data
roottraits_RDA <- rootdata %>% dplyr::select("Plot", "Direction", "PME", "TotalRootBiomass", "SRL", "SRA", "D", "RTD", "Tips", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg")
roottraits_RDA$Sample = paste(roottraits_RDA$Plot, roottraits_RDA$Direction, sep="-")
roottraits_RDA <- roottraits_RDA %>% dplyr::select(-Plot, -Direction) %>% dplyr::select(Sample, everything())
roottraits_RDA <- na.omit(roottraits_RDA)

# join to match complete rows
RDA_fulldata <- soilparameters_RDA %>% inner_join(roottraits_RDA)
soilparameters_RDA <- RDA_fulldata %>% dplyr::select("Sample", "soilDim.1", "soilDim.2") %>% column_to_rownames("Sample")
roottraits_RDA <- RDA_fulldata %>% dplyr::select("Sample", "PME", "TotalRootBiomass", "SRL", "SRA", "D", "RTD", "Tips", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg") %>% column_to_rownames("Sample")

#scale and center soil parameters due to diff units
soilparameters_RDA_trans <- decostand(soilparameters_RDA, method = "standardize")
#variables are now centered around a mean of 0
round(apply(soilparameters_RDA_trans, 2, mean), 1)
#and scaled to have a standard deviation of 1
apply(soilparameters_RDA_trans, 2, sd)
soilparameters_RDA_trans <- soilparameters_RDA_trans %>% rownames_to_column("Sample")
soilparameters_RDA_trans <- soilparameters_RDA_trans %>% tidyr::separate(Sample, into = c("Plot", "Direction"), sep = "-")
soilparameters_RDA_trans <- soilparameters_RDA_trans %>% full_join(ft)
soilparameters_RDA_trans$FT <- factor(soilparameters_RDA_trans$FT)
soilparameters_RDA_trans$Sample = paste(soilparameters_RDA_trans$Plot, soilparameters_RDA_trans$Direction, sep="-")
soilparameters_RDA_trans <- soilparameters_RDA_trans %>% dplyr::select(-Plot, -Direction) %>% dplyr::select(Sample, everything())

#and for root traits
#scale and center soil parameters due to diff units
roottraits_RDA_trans <- decostand(roottraits_RDA, method = "standardize")
#variables are now centered around a mean of 0
round(apply(roottraits_RDA_trans, 2, mean), 1)
#and scaled to have a standard deviation of 1
apply(roottraits_RDA_trans, 2, sd)


##### Run RDA

RDA <- rda(roottraits_RDA_trans~FT + soilDim.1 + soilDim.2, data = soilparameters_RDA_trans)
screeplot(RDA)

R.sum <- summary(RDA)
#explanatory variables explain 10.42% of the variation in root traits 
R.sum$cont   # Prints the "Importance of components" table
R.sum$cont$importance[2, "RDA1"]
# 0.04939833
R.sum$cont$importance[2, "RDA2"]
# 0.03849879

constrained_eig <- RDA$CCA$eig/RDA$tot.chi*100
unconstrained_eig <- RDA$CA$eig/RDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')


#looking for model fit for constrained ordination
RsquareAdj(RDA)$r.squared #unadjusted R^2
RsquareAdj(RDA)$adj.r.squared #adjusted R^2 measures the unbiased amt of explained variation

#so total is 8.25% of root variation explained by soils/FT


##### RDA Plot

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)

#two types of sites scores
#‘lc’: orthogonal linear combinations of the explanatory variable
#‘wa’: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.

#two types of scaling
#Scaling 1- distance biplot (object focused): distance between objects are eudlidean distances (objects closer together have similar variable values), angles between vectors of response variables are meaningless. Angles between vectors of response variables and explanatory variables reflect linear correlation.
#Scaling 2- correlation biplot (response variable focused): distances between objects are not approximate Euclidean distances. Angles between all vectors reflect linear correlation.

#we care about relationships and not values so scaling 2!
plot(RDA, scaling=2, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(RDA, choices=1:2, scaling=2, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')
text(RDA, "species", col="black", cex=0.3, scaling = "sites")


## more aesthetic
roottraits_RDA <- roottraits_RDA %>% rownames_to_column("Sample")

full_ft <- soildata %>% dplyr::select(Plot, Direction, FT)
full_ft$Sample = paste(full_ft$Plot, full_ft$Direction, sep="-")
full_ft <- full_ft %>% inner_join(roottraits_RDA) %>% dplyr::select(Sample, FT)

#colors for FT
full_ft$FT <- as.factor(full_ft$FT)
full_ft$color <- full_ft$FT %>% as.character()
full_ft$color[full_ft$color == "Early secondary"] <- "#dd5129"
full_ft$color[full_ft$color == "Mature secondary"] <- "#0f7ba2"
full_ft$color[full_ft$color == "Old-growth"] <- "#43b284"
summary(full_ft)

## Figure 6A
#generate biplot for root traits and soil var
plot(c(-1.5, 1), c(-3, 1), xlab="RDA1 (47.42%)", ylab="RDA2 (36.96%)", type="n")
abline(h = 0, v = 0, lty = 2, lwd =.5)
text(RDA, dis="cn", cex = 0.8, scaling = 2, col = "blue", pos = 4)
points(RDA, "species", pch=21, col="black", cex=2, scaling = 2)
text(RDA, "species", col="black", cex=0.8, scaling = 2, pos = 2)
RDA_plot <- recordPlot()

png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_6A", 
    width = 750, height = 650, units = "px")
RDA_plot
dev.off()


##### Testing significance

# Test of RDA result
set.seed(123)  # For reproducibility
perm_test <- anova.cca(RDA, permutations = 999, strata = soilparameters_RDA_trans$Plot)
print(perm_test)
#Holm's correction
test_axis.adj <- perm_test
test_axis.adj$`Pr(>F)` <- p.adjust (perm_test$`Pr(>F)`, method = 'holm')
test_axis.adj

# Test of all canonical axes
set.seed(123)  # For reproducibility
perm_test <- anova.cca(RDA, by='axis', permutations = 999, strata = soilparameters_RDA_trans$Plot)
print(perm_test)
#Holm's correction
test_axis.adj <- perm_test
test_axis.adj$`Pr(>F)` <- p.adjust (perm_test$`Pr(>F)`, method = 'holm')
test_axis.adj

#each explanatory variable
set.seed(123)  # For reproducibility
perm_test <- anova.cca(RDA, by='margin', permutations = 999, strata = soilparameters_RDA_trans$Plot)
print(perm_test)
#Holm's correction
test_axis.adj <- perm_test
test_axis.adj$`Pr(>F)` <- p.adjust (perm_test$`Pr(>F)`, method = 'holm')
test_axis.adj


#### Variance partitioning

# subset into FT vs. soil
env_FT <- subset(soilparameters_RDA_trans, select = c(FT))
env_soil <- subset(soilparameters_RDA_trans, select = c(soilDim.1, soilDim.2))

# rda(Y ~ X + Condition(Z))
# perform an RDA of a response matrix "Y" and an explanatory matrix "X" after partialling out the effects of a conditioning matrix "Z"
FT_RDA <- rda(roottraits_RDA_trans~FT + Condition(soilDim.1 + soilDim.2), soilparameters_RDA_trans)
soilPC1_RDA <- rda(roottraits_RDA_trans~soilDim.1 + Condition(FT + soilDim.2), soilparameters_RDA_trans)
soilPC2_RDA <- rda(roottraits_RDA_trans~soilDim.2 + Condition(FT + soilDim.1), soilparameters_RDA_trans)

ordiplot(FT_RDA, scaling = 2, main = "Soil RDA")
ordiplot(soilPC1_RDA, scaling = 2, main = "FT RDA")
ordiplot(soilPC2_RDA, scaling = 2, main = "FT RDA")


#### Variance Partitioning

# Partition the variation in root traits
spe.part.all <- varpart(roottraits_RDA_trans, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_6B", 
    width = 500, height = 500, units = "px")
## Figure 6B
draw.pairwise.venn(area1 = 4.6, area2 = 3.7, cross.area = 0.3, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-0.5,-1.5), c(0.5, -0.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)

dev.off()

#significance testing
anova.cca(rda(roottraits_RDA_trans, env_FT, strata = soilparameters_RDA_trans$Plot)) #FT without controlling for soil
anova.cca(rda(roottraits_RDA_trans, env_FT, env_soil, strata = soilparameters_RDA_trans$Plot)) #FT while controlling for soil
anova.cca(rda(roottraits_RDA_trans, env_soil)) #soil without controlling for FT
anova.cca(rda(roottraits_RDA_trans, env_soil, env_FT)) #soil while controlling for FT
#all are significant

#specific root traits

#SRL
spe.part.all <- varpart(roottraits_RDA_trans$SRL, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
SRL_var <- recordPlot()
#D
spe.part.all <- varpart(roottraits_RDA_trans$D, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
D_var <- recordPlot()
#RTD
spe.part.all <- varpart(roottraits_RDA_trans$RTD, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
RTD_var <- recordPlot()
#biomass
spe.part.all <- varpart(roottraits_RDA_trans$TotalRootBiomass, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
biomass_var <- recordPlot()
#C
spe.part.all <- varpart(roottraits_RDA_trans$rC, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
rC_var <- recordPlot()
#N
spe.part.all <- varpart(roottraits_RDA_trans$rN, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
rN_var <- recordPlot()
#P
spe.part.all <- varpart(roottraits_RDA_trans$rP, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
rP_var <- recordPlot()
#Ca
spe.part.all <- varpart(roottraits_RDA_trans$rCa, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5, textsize = 20)
rCa_var <- recordPlot()
#Tips
spe.part.all <- varpart(roottraits_RDA_trans$Tips, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
tips_var <- recordPlot()
#PME
spe.part.all <- varpart(roottraits_RDA_trans$PME, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
PME_var <- recordPlot()
#Al
spe.part.all <- varpart(roottraits_RDA_trans$rAl, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
Al_var <- recordPlot()
#Fe
spe.part.all <- varpart(roottraits_RDA_trans$rFe, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
Fe_var <- recordPlot()
#K
spe.part.all <- varpart(roottraits_RDA_trans$rK, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
K_var <- recordPlot()
#Mg
spe.part.all <- varpart(roottraits_RDA_trans$rMg, env_FT, env_soil)
spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("FT", "Soil"), # name the partitions
     bg = c("#43b284", "#0f7ba2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
Mg_var <- recordPlot()

## Figure S2 (all below)

#SRL
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/SRL", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 20, area2 = 4.1, cross.area = 3.3, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-2), c(1, -2)), print.mode = "raw", fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#D
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/D", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 8.3, area2 = 3.2, cross.area = 2.3, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-0.5,-1.5), c(1.5, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#RTD
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/RTD", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 6.8, area2 = 3.5, cross.area = 0.8, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-2), c(1.5, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#biomass
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/biomass", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 0.9, area2 = 10.2, cross.area = 0.9, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(0,0), c(1, -1.5)), print.mode = "raw", 
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#rC
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rC", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 3.4, area2 = 1.9, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-0.5,-1.5), c(1.5, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#rN
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rN", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 3.8, area2 = 0.8, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-2), c(1.5, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#rP
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rP", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 1.9, area2 = 2.7, cross.area = 0.3, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-0.5,-1.5), c(1.5, -2)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4, 
                 rotation.degree = 180)
dev.off()

#rCa
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rCa", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 4.4, area2 = 5.1, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-2), c(1.5, -2)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4,
                   rotation.degree = 180)
dev.off()

#tips
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/tips", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 0.2, area2 = 1.7, cross.area = 0.2, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-1), c(1, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4, 
                   rotation.degree = 180)
dev.off()

#PME
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/PME", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 0.5, area2 = 0.4, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-2), c(1, -1.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4)
dev.off()

#Al
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rAl", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 0.3, area2 = 21.9, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(0,-0.5), c(1.5, -2)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4, 
                   rotation.degree = 180)
dev.off()

#Fe
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rFe", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 0, area2 = 5.8, cross.area = 0, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-1,-1.5), c(1.5, -2)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4, 
                   rotation.degree = 180)
dev.off()

#K
grid.newpage()
png(filename = "/Users/chungwingko/Library/Mobile Documents/com~apple~CloudDocs/CCProject/Figures_Final/Fig_S2/rK", 
    width = 500, height = 500, units = "px")
draw.pairwise.venn(area1 = 7.9, area2 = 0.4, cross.area = 0.4, category = c("FT", "Soil"), cex = 1.5, cat.cex = 1.5, euler.d = TRUE, scaled = TRUE, cat.just = list(c(-0.5,-1.5), c(0.5, -0.5)), print.mode = "raw",
                   fill = c("#175f5d","#ae8548"), alpha = 0.4, 
                   rotation.degree = 180)
dev.off()


#### Plot-level partial RDA + tree community

##### Prep data

tree_counts <- data.scores %>% dplyr::select(Plot, NMDS1, NMDS2, NMDS3)

# avg soil data across plot
soilparameters_RDA_avg <- soilcoord %>% dplyr::select(Plot, FT, soilDim.1, soilDim.2)
soilparameters_RDA_avg$FT <- factor(soilparameters_RDA_avg$FT, levels=c("Early secondary", "Mature secondary", "Old-growth"))
soilparameters_RDA_avg$FT <- as.numeric(soilparameters_RDA_avg$FT)
soilparameters_RDA_avg <- soilparameters_RDA_avg %>% group_by(Plot) %>% summarise_all(mean, na.rm = TRUE)
soilparameters_RDA_avg <- soilparameters_RDA_avg %>% tibble::column_to_rownames("Plot")
#soilparameters_RDA_avg <- soilparameters_RDA_avg[!(row.names(soilparameters_RDA_avg) %in% c("W04","WC1")),]
soilparameters_RDA <- na.omit(soilparameters_RDA)

# avg root data across plot
roottraits_RDA_avg <- rootdata %>% dplyr::select("Plot", "PME", "TotalRootBiomass", "SRL", "SRA", "D", "RTD", "Tips", "rC", "rN", "rP", "rAl", "rCa", "rFe", "rK", "rMg")
roottraits_RDA_avg <- roottraits_RDA_avg %>% group_by(Plot) %>% summarise_all(mean, na.rm = TRUE)
roottraits_RDA_avg <- roottraits_RDA_avg %>% tibble::column_to_rownames("Plot")
#roottraits_RDA_avg <- roottraits_RDA_avg[!(row.names(roottraits_RDA_avg) %in% c("W04","WC1")),]
roottraits_RDA_avg <- na.omit(roottraits_RDA_avg)

#scale and center soil parameters due to diff units
soilparameters_avg_trans <- decostand(soilparameters_RDA_avg, method = "standardize")
#variables are now centered around a mean of 0
round(apply(soilparameters_avg_trans, 2, mean), 1)
#and scaled to have a standard deviation of 1
apply(soilparameters_avg_trans, 2, sd)
soilparameters_avg_trans <- soilparameters_avg_trans %>% rownames_to_column("Plot")

#combine full explanatory matrix
combined_treesoil <- tree_counts %>% full_join(soilparameters_avg_trans)

#and for root traits
#scale and center soil parameters due to diff units
roottraits_RDA_trans <- decostand(roottraits_RDA_avg, method = "standardize")
#variables are now centered around a mean of 0
round(apply(roottraits_RDA_trans, 2, mean), 1)
#and scaled to have a standard deviation of 1
apply(roottraits_RDA_trans, 2, sd)


##### Partial RDA time

# subset into FT vs. soil
env_tree <- subset(combined_treesoil, select = c(NMDS1, NMDS2, NMDS3))
env_soil <- subset(combined_treesoil, select = c(soilDim.1, soilDim.2))
env_FT <- subset(combined_treesoil, select = c(FT))

# rda(Y ~ X + Condition(Z))
# perform an RDA of a response matrix "Y" and an explanatory matrix "X" after partialling out the effects of a conditioning matrix "Z"
soil_RDA <- rda(roottraits_RDA_trans~soilDim.1 + soilDim.2 + Condition(FT + NMDS1 + NMDS2 + NMDS3), combined_treesoil)
tree_RDA <- rda(roottraits_RDA_trans~NMDS1 + NMDS2 + NMDS3 + Condition(FT + soilDim.1 + soilDim.2), combined_treesoil)
FT_RDA <- rda(roottraits_RDA_trans~FT + Condition(soilDim.1 + soilDim.2 + NMDS1 + NMDS2 + NMDS3), combined_treesoil)

# Partition the variation in root traits
spe.part.all <- varpart(roottraits_RDA_trans, ~NMDS1 + NMDS2 + NMDS3, ~FT, ~soilDim.1 + soilDim.2, data = combined_treesoil)

spe.part.all$part  # access results!
# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("Tree", "FT", "Soil"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "brown"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

#significance testing
anova.cca(rda(roottraits_RDA_trans, env_FT)) #FT without controlling for soil
anova.cca(rda(roottraits_RDA_trans, env_tree)) #FT without controlling for soil
anova.cca(rda(roottraits_RDA_trans, env_FT, env_soil)) #FT while controlling for soil
anova.cca(rda(roottraits_RDA_trans, env_soil)) #soil without controlling for FT
anova.cca(rda(roottraits_RDA_trans, env_soil, env_FT)) #soil while controlling for FT
#all are significant
