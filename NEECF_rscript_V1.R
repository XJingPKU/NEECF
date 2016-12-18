#=======================================
# Data analyses for soil enzyme activities in NEECF
# 
# First created by Xin Jing, Oct. 20, 2016
# Last modified: December 16, 2016
# 
# R version 3.3.2 (2016-10-31)
# Email: jingxin0123@gmail.com
#====================================== 
rm(list = ls())

#======================================
# load library
library("plyr") # calls: ddply
library("dplyr") # calls: select
library("ggplot2") # calls: ggplot
library("reshape2") # calls: melt
library("asreml") # calls: asreml
## install 'pascal' package
# library(devtools)
# install_github("pascal-niklaus/pascal/pascal")
library("pascal") # calls: test.asreml, aov.ftest
library("ggbiplot") # calls: ggbiplot
library("gridExtra") # calls: grid.arrange

#======================================
# control the number of digits
options(digits = 4)

#======================================
# read table and tidy the data
enz.neecf <- read.csv("./data/data_processing/enz_forest_all.csv", header = T)

# transform variable types
enz.neecf$fN <- as.factor(enz.neecf$N)
enz.neecf$fb <- as.factor(enz.neecf$new.block)
enz.neecf$fp <- as.factor(enz.neecf$new.plot)
enz.neecf$fD <- as.factor(enz.neecf$depth)

# convert negative values of soil enzyme activity into zero
enz.vars <- c("BG", "CB", "NAG", "PHOS", "LAP", "POX", "PER")
enz.neecf[names(enz.neecf) %in% enz.vars][enz.neecf[names(enz.neecf) %in% enz.vars] < 0] <- 0

# data summary
enz.neecf.sum <- enz.neecf %>%
  select(site, forest.type, fD, fN, BG, CB, POX, PER, NAG, LAP, PHOS)
levels(enz.neecf.sum$fN) <- c("CK", "N20", "N25", "N50", "N100", "N150")
levels(enz.neecf.sum$fD) <- c("0cm", "0-10cm", "10-20cm", "20-40cm", "40-60cm")
enz.neecf.sum <- enz.neecf.sum[enz.neecf.sum$fN %in% c("CK", "N50", "N100"), ]
enz.neecf.sum <- enz.neecf.sum[enz.neecf.sum$fD %in% c("0-10cm", "10-20cm", "20-40cm", "40-60cm"), ]
enz.neecf.sum <- enz.neecf.sum[enz.neecf.sum$forest.type != "secondary_tropical_montane_rain_forest"
                               & enz.neecf.sum$forest.type != "birch_forest", ]
enz.neecf.sum <- droplevels(enz.neecf.sum)
df.enz.neecf.sum <- melt(enz.neecf.sum, 
                         id.vars = c("site", "forest.type", "fD", "fN"))

df.enz.neecf.summary <- ddply(df.enz.neecf.sum, c("site", "fD", "variable", "fN"),
      summarise,
      mean.x = mean(value, na.rm = T),
      sd.x = sd(value, na.rm = T))
# write.csv(df.enz.neecf.results, "./outputs/NEECF_data_summary.csv")

#============================

# C:N, C:P, and NP ratio
enz.neecf$CNratio <- log(enz.neecf$BG)/log(enz.neecf$NAG + enz.neecf$LAP)
enz.neecf$CPratio <- log(enz.neecf$BG)/log(enz.neecf$PHOS)
enz.neecf$NPratio <- log(enz.neecf$NAG + enz.neecf$LAP)/log(enz.neecf$PHOS)

# enzyme activity transformation
enz.neecf$BG <- log(enz.neecf$BG)
enz.neecf$CB <- log(enz.neecf$CB)
enz.neecf$NAG <- log(enz.neecf$NAG)
enz.neecf$LAP <- log(enz.neecf$LAP)
enz.neecf$PHOS <- log(enz.neecf$PHOS)
enz.neecf$POX <- log(enz.neecf$POX)
enz.neecf$PER <- log(enz.neecf$PER)

# convert "Inf" into "NA"
enz.neecf[mapply(is.infinite, enz.neecf)] <- NA

# convert "NaN" into "NA"
enz.neecf[mapply(is.nan, enz.neecf)] <- NA

# deal with the outliers of CNratio
enz.neecf$PER[enz.neecf$PER < -2.5] <- NA
enz.neecf$CNratio[enz.neecf$CNratio > 2.2 | enz.neecf$CNratio < -0.5] <- NA
enz.neecf$CPratio[enz.neecf$CPratio > 1.5 | enz.neecf$CPratio < -0.5] <- NA
enz.neecf$NPratio[enz.neecf$NPratio > 1.5 | enz.neecf$NPratio < -0.5] <- NA

# extract subset of enz.neecf CK, N50, N100
enz.neecf <- subset(enz.neecf, treat %in% c("CK", "N50", "N100"))

# extract subset of enz.neecf c("0-10cm", "10-20cm", "20-40cm", "40-60cm")
enz.neecf <- subset(enz.neecf, 
                    soil.layer %in% c("0-10cm", "10-20cm", "20-40cm", "40-60cm"))
# remove secondary forest
enz.neecf <- subset(enz.neecf, 
                    forest.type != "secondary_tropical_montane_rain_forest" &
                      forest.type != "birch_forest")
enz.neecf <- droplevels(enz.neecf)

###########################################################
# experiment design
xtabs(~ site + fD, data = enz.neecf)
xtabs(~ site + fb, data = enz.neecf)

###########################################################
# PCA analysis
pc.vars <- c("BG", "CB", "POX", "PER",
             "NAG", "LAP", "PHOS")
(enz.pca <- prcomp(na.omit(enz.neecf[pc.vars]), 
                  center = TRUE, scale = TRUE))
plot(enz.pca, type = "l")
summary(enz.pca)
grp <- na.omit(enz.neecf)$site
ggbiplot(enz.pca, choices = 1:2, obs.scale = 1, var.scale = 1, 
         varname.size = 4.5, group = grp, 
         ellipse = TRUE) +
  geom_point(aes(color = grp), size = 3) +
  geom_vline(xintercept = 0, color = "gray", size = 0.8) +
  geom_hline(yintercept = 0, color = "gray", size = 0.8) +
  scale_color_discrete(name = '') +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")
# ggsave("./outputs/NEECF_PCA.pdf", width = 6, height = 6.0)

###########################################################
# bivariate correlations between enzyme activities and latitude 
# select subdata
enz.df <- enz.neecf %>%
  select(site, latitude, longitude, altitude, MAT, MAP, 
         fN, fb, fp, fD, 
         BG, CB, POX, PER, NAG, LAP, PHOS,  
         CNratio, CPratio, NPratio)

# melt the subdata
enz.df <- melt(enz.df, 
               id.vars = c("site", "latitude", "longitude", "altitude", "MAT", "MAP", 
                           "fN", "fb", "fp", "fD"))
levels(enz.df$fD) <- c("0-10 cm", "10-20 cm", "20-40 cm", "40-60 cm")
levels(enz.df$fN) <- c("CK", "N50", "N100")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract subset of the data
# by act
enz.df.act <- enz.df %>%
  filter(variable %in% c("BG", "CB", "NAG", "LAP", "PHOS", "POX", "PER")) 
# by ratio
enz.df.ratio <- enz.df %>%
  filter(variable %in% c("CNratio", "CPratio", "NPratio")) 
enz.df.ratio <- droplevels(enz.df.ratio)
levels(enz.df.ratio$variable) <- c("C:N ratio", "C:P ratio", "N:P ratio")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# significance test for bivariate correlations
# number of values used for the significance test
ddply(enz.df.act, c("fD", "variable"), nrow)
# aov.ftest to corect F-test
aovFtest <- function(x) {
  mod <- aov(value ~ latitude + site, data = x)
  aov.ftest(mod, latitude ~ site, table = TRUE)[c(7, 8)]
}

enz.df.act.aovftest <- ddply(enz.df.act, c("fD", "variable"), aovFtest)
enz.df.act.aovftest$labels <- with(enz.df.act.aovftest,
                                   paste("F = ", F, "p =", P))
enz.df.act.aovftest$xpos <- rep(32.5, 28)
enz.df.act.aovftest$ypos <- rep(c(-1.0, 4.65, 3.0, -1.0, 6.0, -2.0, 2.2), each = 4)

# plot bivariate correlations with ggplot
# by depth
# enzyme activity
ggplot(enz.df.act, aes(x = latitude, y = value, color = fD)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fD, scales = "free_y") +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "Activity") +
  geom_text(data = enz.df.act.aovftest, aes(x = xpos, y = ypos, label = labels, 
                                     group = NULL, hjust = 0.5, vjust = 0.5)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12))
# ggsave("./outputs/NEECF_enz_lat_depth_act.pdf", height = 10, width = 8)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ratio
# number of values used for the significance test
ddply(enz.df.ratio, c("fD", "variable"), nrow)
# aov.ftest to corect F-test
enz.df.ratio.aovftest <- ddply(enz.df.ratio, c("fD", "variable"), aovFtest)
enz.df.ratio.aovftest$labels <- with(enz.df.ratio.aovftest,
                                   paste("F = ", F, "p =", P))
enz.df.ratio.aovftest$xpos <- rep(35, 12)
enz.df.ratio.aovftest$ypos <- rep(c(-0.2, 0.15, -0.25), 4)
p1 <- ggplot(enz.df.ratio, aes(x = latitude, y = value, color = fD)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fD, scales = "free_y") +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "") +
  geom_text(data = enz.df.ratio.aovftest, aes(x = xpos, y = ypos, label = labels, 
                                            group = NULL, hjust = 0.5, vjust = 0.5)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12)) 
# ggsave("./outputs/NEECF_enz_lat_depth_ratio.pdf", height = 4.5, width = 8)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CN vs CP ratio
ratio.fit <- dlply(enz.neecf, "fD", function(x) lm(CNratio ~ CPratio, data = x))
ratio.fit
ratio.labels <- data.frame(xpos = rep(0.4, 4),
                           ypos = rep(-0.25, 4),
                           fD = levels(enz.neecf$fD),
                           r2 = sapply(ratio.fit, function(x) summary(x)$r.squared),
                           p = sapply(ratio.fit, function(x) anova(x)[1, 5])
                           )
ratio.labels$labels <- with(ratio.labels,
                            paste(expression(R^2), "=", round(r2, 2), 
                                  "p =", round(p, 3), sep = " "))
p2 <- ggplot(data = enz.neecf, aes(x = CPratio, y = CNratio, color = fD)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(~ fD, scales = "free_y") +
  geom_text(data = ratio.labels, aes(x = xpos, y = ypos, label = labels, 
                                  group = NULL, hjust = 0.5, vjust = 0.5)) +
  labs(x = "C:P ratio [ln(BG):ln(PHOS)]", 
       y = "C:N ratio [ln(BG):ln(NAG+LAP)]") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.key = element_rect("NA"),
        axis.title = element_text(size = 12))
# pdf("./outputs/NEECF_enz_lat_depth_ratio.pdf", width = 8, height = 10)
gridExtra::grid.arrange(p1, p2, ncol = 1,
                        heights = c(0.73, 0.27))
# dev.off()
###########################################################
# treatment and depth effects
# fitting linear mixed effects model in asreml
vars <- c("BG", "CB", "POX", "PER",
          "NAG", "LAP", "PHOS",
          "CNratio", "CPratio", "NPratio")
for (i in vars) {
  print(i)
  df.test <- enz.df %>% filter(value != 0 & variable == i)
  fit <- asreml(fixed =  value ~ latitude + fN + fD + latitude:fN + latitude:fD + fN:fD, 
                random =  ~ site/(fb + fD + fN),
                keep.order = TRUE,
                control = asreml.control(maxiter = 500),
                data = df.test)
  test.asreml(fit)
  # AIC
  cat("\n----AIC:\n")
  print(asremlPlus::info.crit.asreml(fit)[2])
  # cat("\n---- Coefficients:\n")
  # print(summary(fit, all = T)$coef.fixed)
  cat("\n####################################################\n")
}
boxplot(value ~ fD*fN, data = enz.df %>% filter(value != 0 & variable == "PHOS"))


###########################################################
# END script                                              #
###########################################################