#==========================================================
# Data analyses for soil enzyme activities in NEECF
# 
# First created by Xin Jing, Oct. 20, 2016
# Last modified: March 5, 2017
# 
# R version 3.3.2 (2016-10-31)
# Email: jingxin0123@gmail.com
#==========================================================
rm(list = ls())

#==========================================================
# load library
library(plyr)  # calls: ddplyr
library("multcomp")  # tukey HSD Note 'multcomp' need to load before 'dplyr'
library("dplyr")  # calls: select, filter, mutate...
library("reshape2")  # calls: melt
library("asreml")  # calls: asreml
library("pascal")  # calls: test.asreml
library("ggplot2")  # calls: ggplot
library("gridExtra")  # calls: grid.arrange, arrangeGrob
library("grid")  # calls: textGrob, gpar
library("psych")  # calls: corr.test
require(ggrepel) # calls: geom_text_repel

#==========================================================
# read table
enz.neecf <- read.csv("./data/data_processing/enz_forest_all_V1.csv")

# clean the data
enz.neecf.clean <- enz.neecf %>%
  select(site, latitude, MAT, MAP, forest.type, treat, new.block, new.plot,
         soil.layer, STC, STN, STP, pH, BG, CB, NAG, PHOS, LAP, POX, PER) %>%
  filter(forest.type != "birch_forest"
         & forest.type != "secondary_tropical_montane_rain_forest"
         & soil.layer != "O-horizon"
         & treat != "N20" & treat != "N25" & treat != "N50" & treat != "P50" 
         & treat != "N50_P50" & treat != "N100_P50" & treat != "N150") %>%
  mutate(site = factor(site, levels = c("JFL", "WYS", "GNJ", "DLS", "WY", "GH")),
         lat = latitude,
         fb = factor(new.block),
         fp = factor(new.plot),
         fd = factor(soil.layer, levels = c("0-10cm", "10-20cm", "20-40cm", "40-60cm")),
         AP = PHOS) %>%
  select(site, lat, MAT, MAP, treat, fb, fp, fd, STC, STN, STP, pH,
         BG, CB, POX, PER, NAG, LAP, AP)

# drop levels
enz.neecf.clean <- droplevels(enz.neecf.clean)

# convert negative values of soil enzyme activity into zero
enz.vars <- c("BG", "CB", "POX", "PER", "NAG", "LAP", "AP")
enz.neecf.clean[names(enz.neecf.clean) 
                %in% enz.vars][enz.neecf.clean[names(enz.neecf.clean) 
                                               %in% enz.vars] < 0] <- 0

# convert "NaN" into "NA"
enz.neecf.clean[mapply(is.nan, enz.neecf.clean)] <- NA

#==========================================================
# supplementary table S2 -- summary of the enzyme activities
# reshape the data
df.summary <- melt(enz.neecf.clean,
                   id.vars = c("site", "lat", "MAT", "MAP", 
                               "treat", "fb", "fp", "fd", "STC", "STN", "STP", "pH"))
# relevels of variable
df.summary$variable <- factor(df.summary$variable,
                              levels = c("BG", "CB", "POX", "PER", "NAG", "LAP", "AP"))

# data summary
df.summary.res <- ddply(df.summary, c("site", "fd", "variable", "treat"),
                        summarise,
                        mean.x = mean(value, na.rm = T),
                        se.x = sd(value, na.rm = T)/length(value))
# write.csv(df.summary.res, "./outputs/NEECF_supplementary_table_data_summary.csv")

#==========================================================
# Table 2: linear mixed-effects models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# by enzyme activities
# data transformation
df.enz <- df.summary %>%
  mutate(logenz = log(value + 1))

# summary of the linear mixed-effects models
for (i in enz.vars) {
  print(i)
  df.test <- df.enz %>% filter(variable == i)
  fit <- asreml(fixed =  logenz ~ site + treat + fd +
                  site:treat + site:fd + treat:fd, 
                random =  ~ fb/fp,
                keep.order = TRUE,
                control = asreml.control(maxiter = 500),
                data = df.test)
  test.asreml(fit)
  # AIC
  # cat("\n----AIC:\n")
  # print(asremlPlus::info.crit.asreml(fit)[2])
  # cat("\n---- Coefficients:\n")
  # print(summary(fit, all = T)$coef.fixed)
  cat("\n####################################################\n")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# by ecoenzyme ratios
df.enz.ratio <- enz.neecf.clean %>%
  mutate(CNratio = log(BG + 1)/log(NAG + LAP + 1),
         CPratio = log(BG + 1)/log(AP + 1),
         NPratio = log(NAG + LAP + 1)/log(AP + 1)) %>%
  select(site, lat, MAT, MAP, treat, fb, fp, fd, CNratio, CPratio, NPratio)

# convert "Inf" into "NA"
df.enz.ratio[mapply(is.infinite, df.enz.ratio)] <- NA

# reshape the data
df.enz.ratio <- melt(df.enz.ratio,
                     id.vars = c("site", "lat", "MAT", "MAP", 
                                 "treat", "fb", "fp", "fd"))
# summary of the linear mixed-effects models
enz.ratio <- c("CNratio", "CPratio", "NPratio")
for (i in enz.ratio) {
  print(i)
  df.test <- df.enz.ratio %>% filter(variable == i)
  fit <- asreml(fixed =  value ~ site + treat + fd +
                  site:treat + site:fd + treat:fd, 
                random =  ~ fb/fp,
                keep.order = TRUE,
                control = asreml.control(maxiter = 500),
                data = df.test)
  test.asreml(fit)
  cat("\n####################################################\n")
}

#==========================================================
# Table 3 corr.test
# clean the data
df.corr.test <- enz.neecf.clean %>%
  select(lat, MAT, MAP, treat, fd, STC, STN, STP, pH,
         BG, CB, POX, PER, NAG, LAP, AP) %>%
  filter(treat == "CK" & fd == "0-10cm") %>%
  mutate(soilCN = STC/STN,
         soilCP = STC/STP,
         soilNP = STN/STP,
         enzCN = log(BG + 1)/log(NAG + LAP + 1),
         enzCP = log(BG + 1)/log(AP + 1),
         enzNP = log(NAG + LAP + 1)/log(AP + 1)) %>%
  select(lat, MAT, MAP, STC, STN, STP, soilCN, soilCP, soilNP, pH,
         BG, CB, POX, PER, NAG, LAP, AP, enzCN, enzCP, enzNP)
# corr.test
# by pearson
corr.test(df.corr.test, use = "pairwise", method = "pearson",
          adjust = "fdr", alpha = .05)
# by spearman
corr.test(df.corr.test, use = "pairwise", method = "spearman",
          adjust = "fdr", alpha = .05)

#==========================================================
# Figure 1: enzyme activities at each site*treatment*depth
# summarise the data
df.fig1 <- ddply(df.enz, c("site", "fd", "variable", "treat"),
                 summarise,
                 mean.x = mean(logenz, na.rm = T),
                 se.x = sd(logenz, na.rm = T)/length(logenz),
                 ci.x = qt(0.95, (length(is.na(logenz) == FALSE) - 1))*se.x)  # 90 percent confidence interval
# Tukey's HSD test
dep.vars <- c("0-10cm", "10-20cm", "20-40cm", "40-60cm")
out <- NULL
for (i in enz.vars) {
  for (j in dep.vars) {
    x <- df.enz %>% filter(variable == i, fd == j)
    fit <- aov(logenz ~ site, data = x)
    tuk.lett <- cld(glht(fit, linfct = mcp(site = "Tukey")))$mcletters$Letters
    tuk.lett <- as.data.frame(tuk.lett)
    res <- data.frame(variable = rep(i, length(tuk.lett)),
                      fd = rep(j, length(tuk.lett)),
                      site = row.names(tuk.lett),
                      tuk = tuk.lett[, 1],
                      treat = rep("CK", length(tuk.lett)))
    out <- rbind(out, res)
  }
}
# paired t test
out.pair <- ddply(df.enz, .(variable, fd, site), function(x) {
  g1 <- x %>% filter(treat == "CK")
  g2 <- x %>% filter(treat == "N100")
  mod <- t.test(g1$logenz, g2$logenz, paired = TRUE)
  res <- mod$p.value
})
out.pair <- data.frame(variable = out.pair$variable, 
                       fd = out.pair$fd,
                       site = out.pair$site,
                       p.val = round(out.pair$V1, 2),
                       treat = rep("CK", length(row.names(out.pair))))
out.pair$sig <- ifelse(out.pair$p.val < 0.05, "*", 
                       ifelse(out.pair$p.val < 0.10, "†", NA))

# plot the resut 
df.fig1 %>%
  ggplot(aes(x = site, y = mean.x, group = treat, shape = treat)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c(1, 19)) +
  geom_errorbar(aes(ymin = mean.x - ci.x, ymax = mean.x + ci.x), 
                width = 0.2, size = 0.25, 
                position = position_dodge(width = 0.7)) +
  facet_grid(variable ~ fd, scales = "free_y") +
  geom_text(data = out, aes(x = site, y = 10, label = tuk, vjust = 1.0), size = 3.5) +
  geom_text(data = out.pair, aes(x = site, y = 7.5, label = sig), size = 4.5) +
  labs(x = "Sites", y = "Enzyme activities [ln(x + 1)]") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))
# ggsave("./outputs/NEECF_enz_figure1.pdf", height = 11.5, width = 9)

#==========================================================
# Figure S1: enzyme activities at each site*treatment*depth
# summarise the data
df.figs1 <- ddply(df.enz, c("fd", "site", "variable", "treat"),
                  summarise,
                  mean.x = mean(logenz, na.rm = T),
                  se.x = sd(logenz, na.rm = T)/length(logenz),
                  ci.x = qt(0.95, (length(is.na(logenz) == FALSE) - 1))*se.x)  # 90 percent confidence interval
# Tukey's HSD test
site.vars <- c("JFL", "WYS", "GNJ", "DLS", "WY",  "GH")
out <- NULL
for (i in enz.vars) {
  for (j in site.vars) {
    x <- df.enz %>% filter(variable == i, site == j)
    fit <- aov(logenz ~ fd, data = x)
    tuk.lett <- cld(glht(fit, linfct = mcp(fd = "Tukey")))$mcletters$Letters
    tuk.lett <- as.data.frame(tuk.lett)
    res <- data.frame(variable = rep(i, length(tuk.lett)),
                      site = rep(j, length(tuk.lett)),
                      fd = row.names(tuk.lett),
                      tuk = tuk.lett[, 1],
                      treat = rep("CK", length(tuk.lett)))
    out <- rbind(out, res)
  }
}

# plot the resut 
df.figs1 %>%
  ggplot(aes(x = fd, y = mean.x, group = treat, shape = treat)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c(1, 19)) +
  geom_errorbar(aes(ymin = mean.x - ci.x, ymax = mean.x + ci.x), 
                width = 0.2, size = 0.25, 
                position = position_dodge(width = 0.7)) +
  facet_grid(variable ~ site, scales = "free_y") +
  geom_text(data = out, aes(x = fd, y = 10, label = tuk, vjust = 1.0), size = 4.5) +
  geom_text(data = out.pair, aes(x = fd, y = 7.5, label = sig), size = 5.5) +
  labs(x = "Soil depth", y = "Enzyme activities [ln(x + 1)]") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
# ggsave("./outputs/NEECF_enz_figureS1.pdf", height = 10.88, width = 13.34)

#==========================================================
# Figure 2: enzyme ratios at each site*treatment*depth
# summarise the data
df.fig2 <- ddply(df.enz.ratio, c("site", "fd", "variable", "treat"),
                 summarise,
                 mean.x = mean(value, na.rm = T),
                 se.x = sd(value, na.rm = T)/length(value),
                 ci.x = qt(0.95, (length(is.na(value) == FALSE) - 1))*se.x)  # 90 percent confidence interval
# Tukey's HSD test
out2 <- NULL
for (i in enz.ratio) {
  for (j in dep.vars) {
    x <- df.enz.ratio %>% filter(variable == i, fd == j)
    fit <- aov(value ~ site, data = x)
    tuk.lett <- cld(glht(fit, linfct = mcp(site = "Tukey")))$mcletters$Letters
    tuk.lett <- as.data.frame(tuk.lett)
    res <- data.frame(variable = rep(i, length(tuk.lett)),
                      fd = rep(j, length(tuk.lett)),
                      site = row.names(tuk.lett),
                      tuk = tuk.lett[, 1],
                      treat = rep("CK", length(tuk.lett)))
    out2 <- rbind(out2, res)
  }
}
# paired t test
df.enz.ratio[mapply(is.nan, df.enz.ratio)] <- NA
myTtest <- function(...) {
  # Silently return NA when there are not enough 'x' observations
  # I copied this t.test funciton from https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html
  obj <- try(t.test(...), silent = TRUE)
  if (is(obj, "try-error")) return(NA)
  else return(obj$p.value)
}
out.pair2 <- ddply(df.enz.ratio, .(variable, fd, site), function(x) {
  g1 <- x %>% filter(treat == "CK")
  g2 <- x %>% filter(treat == "N100")
  mod <- myTtest(g1$value, g2$value, paired = TRUE)
  # res <- mod$p.value
})
out.pair2 <- data.frame(variable = out.pair2$variable, 
                        fd = out.pair2$fd,
                        site = out.pair2$site,
                        p.val = round(out.pair2$V1, 2),
                        treat = rep("CK", length(row.names(out.pair2))))
out.pair2$sig <- ifelse(out.pair2$p.val < 0.05, "*", 
                        ifelse(out.pair2$p.val < 0.10, "†", NA))


# plot the resut 
df.fig2 %>%
  ggplot(aes(x = site, y = mean.x, group = treat, shape = treat)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c(1, 19)) +
  geom_errorbar(aes(ymin = mean.x - ci.x, ymax = mean.x + ci.x), 
                width = 0.2, size = 0.25, 
                position = position_dodge(width = 0.7)) +
  facet_grid(variable ~ fd, scales = "free_y") +
  geom_text(data = out2, aes(x = site, y = 5, label = tuk, hjust = 1.0), size = 3.5) +
  geom_text(data = out.pair2, aes(x = site, y = 3.5, label = sig), size = 4.5) +
  labs(x = "Sites", y = "Enzyme ratios") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))
# ggsave("./outputs/NEECF_enz_figure2.pdf", height = 5.5, width = 9)

#==========================================================
# Figure S2: enzyme ratios at each site*treatment*depth
# summarise the data
df.figs2 <- ddply(df.enz.ratio, c("fd", "site", "variable", "treat"),
                  summarise,
                  mean.x = mean(value, na.rm = T),
                  se.x = sd(value, na.rm = T)/length(value),
                  ci.x = qt(0.95, (length(is.na(value) == FALSE) - 1))*se.x)  # 90 percent confidence interval
# Tukey's HSD test
out2 <- NULL
for (i in enz.ratio) {
  for (j in site.vars) {
    x <- df.enz.ratio %>% filter(variable == i, site == j)
    fit <- aov(value ~ fd, data = x)
    tuk.lett <- cld(glht(fit, linfct = mcp(fd = "Tukey")))$mcletters$Letters
    tuk.lett <- as.data.frame(tuk.lett)
    res <- data.frame(variable = rep(i, length(tuk.lett)),
                      site = rep(j, length(tuk.lett)),
                      fd = row.names(tuk.lett),
                      tuk = tuk.lett[, 1],
                      treat = rep("CK", length(tuk.lett)))
    out2 <- rbind(out2, res)
  }
}

# plot the resut 
df.figs2 %>%
  ggplot(aes(x = fd, y = mean.x, group = treat, shape = treat)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c(1, 19)) +
  geom_errorbar(aes(ymin = mean.x - ci.x, ymax = mean.x + ci.x), 
                width = 0.2, size = 0.25, 
                position = position_dodge(width = 0.7)) +
  facet_grid(variable ~ site, scales = "free_y") +
  geom_text(data = out2, aes(x = fd, y = 5, label = tuk, vjust = 1.0), size = 4.5) +
  geom_text(data = out.pair2, aes(x = fd, y = 3.5, label = sig), size = 5.5) +
  labs(x = "Soil depth", y = "Enzyme ratios") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
# ggsave("./outputs/NEECF_enz_figureS2.pdf", height = 6.5, width = 13.34)

#==========================================================
# bivariate correlations between enzyme activities and latitude 

# significance test for bivariate correlations
# number of values used for the significance test
ddply(df.enz, c("fd", "variable"), nrow)
# aov.ftest to corect F-test
aovFtest <- function(x) {
  mod <- aov(logenz ~ lat + site, data = x)
  aov.ftest(mod, lat ~ site, table = TRUE)[c(7, 8)]
}

enz.df.act.aovftest <- ddply(df.enz, c("fd", "variable"), aovFtest)
enz.df.act.aovftest$labels <- with(enz.df.act.aovftest,
                                   paste("F = ", F, "P =", P))
enz.df.act.aovftest$xpos <- rep(32.5, 28)
enz.df.act.aovftest$ypos <- rep(c(8, 6, 6, 6, 8, 6, 12), 4)

# plot bivariate correlations with ggplot
# by depth
# enzyme activity
ggplot(df.enz, aes(x = lat, y = logenz, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "Activity") +
  geom_text(data = enz.df.act.aovftest, 
            aes(x = xpos, y = ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12))
# ggsave("./outputs/NEECF_enz_lat_depth_act.pdf", height = 10, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# significance test for bivariate correlations by MAT
# aov.ftest to corect F-test
aovFtest.mat <- function(x) {
  mod <- aov(logenz ~ MAT + site, data = x)
  aov.ftest(mod, MAT ~ site, table = TRUE)[c(7, 8)]
}

enz.df.act.aovftest <- ddply(df.enz, c("fd", "variable"), aovFtest.mat)
enz.df.act.aovftest$labels <- with(enz.df.act.aovftest,
                                   paste("F = ", F, "P =", P))
enz.df.act.aovftest$xpos <- rep(8, 28)
enz.df.act.aovftest$ypos <- rep(c(8, 6, 6, 6, 8, 6, 12), 4)

# plot bivariate correlations with ggplot
# by depth
# enzyme activity
ggplot(df.enz, aes(x = MAT, y = logenz, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = expression("Mean annual temperature ("*degree*C*")"), 
       y = "Activity") +
  geom_text(data = enz.df.act.aovftest, 
            aes(x = xpos, y = ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12))
# ggsave("./outputs/NEECF_enz_lat_depth_act_MAT.pdf", height = 10, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# significance test for bivariate correlations by MAP
# aov.ftest to corect F-test
aovFtest.map <- function(x) {
  mod <- aov(logenz ~ MAP + site, data = x)
  aov.ftest(mod, MAP ~ site, table = TRUE)[c(7, 8)]
}

enz.df.act.aovftest <- ddply(df.enz, c("fd", "variable"), aovFtest.map)
enz.df.act.aovftest$labels <- with(enz.df.act.aovftest,
                                   paste("F = ", F, "P =", P))
enz.df.act.aovftest$xpos <- rep(1400, 28)
enz.df.act.aovftest$ypos <- rep(c(8, 6, 6, 6, 8, 6, 12), 4)

# plot bivariate correlations with ggplot
# by depth
# enzyme activity
ggplot(df.enz, aes(x = MAP, y = logenz, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = "Mean annual precipitation (mm)", 
       y = "Activity") +
  geom_text(data = enz.df.act.aovftest, 
            aes(x = xpos, y = ypos, label = labels,
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12))
# ggsave("./outputs/NEECF_enz_lat_depth_act_MAP.pdf", height = 10, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ratio
# aov.ftest to corect F-test
aovFtest <- function(x) {
  mod <- aov(value ~ lat + site, data = x)
  aov.ftest(mod, lat ~ site, table = TRUE)[c(7, 8)]
}
enz.df.ratio.aovftest <- ddply(df.enz.ratio, c("fd", "variable"), aovFtest)
enz.df.ratio.aovftest$labels <- with(enz.df.ratio.aovftest,
                                     paste("F = ", F, "P =", P))
enz.df.ratio.aovftest$xpos <- rep(35, 12)
enz.df.ratio.aovftest$ypos <- rep(c(5.25, 1.5, 1.25), 4)
ggplot(df.enz.ratio, aes(x = lat, y = value, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = expression("Latitude ("*degree*N*")"), 
       y = "") +
  geom_text(data = enz.df.ratio.aovftest, 
            aes(x = xpos, y = ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12)) 
# ggsave("./outputs/NEECF_enz_lat_depth_ratio.pdf", height = 4.5, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# aov.ftest to corect F-test
aovFtest.mat <- function(x) {
  mod <- aov(value ~ MAT + site, data = x)
  aov.ftest(mod, MAT ~ site, table = TRUE)[c(7, 8)]
}
enz.df.ratio.aovftest <- ddply(df.enz.ratio, c("fd", "variable"), aovFtest.mat)
enz.df.ratio.aovftest$labels <- with(enz.df.ratio.aovftest,
                                     paste("F = ", F, "P =", P))
enz.df.ratio.aovftest$xpos <- rep(8, 12)
enz.df.ratio.aovftest$ypos <- rep(c(5.25, 1.5, 1.25), 4)
ggplot(df.enz.ratio, aes(x = MAT, y = value, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = expression("Mean annual temperature ("*degree*C*")"), 
       y = "") +
  geom_text(data = enz.df.ratio.aovftest, 
            aes(x = xpos, y = ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12)) 
# ggsave("./outputs/NEECF_enz_lat_depth_ratio_MAT.pdf", height = 4.5, width = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ratio
# aov.ftest to corect F-test
aovFtest.map <- function(x) {
  mod <- aov(value ~ MAP + site, data = x)
  aov.ftest(mod, MAP ~ site, table = TRUE)[c(7, 8)]
}
enz.df.ratio.aovftest <- ddply(df.enz.ratio, c("fd", "variable"), aovFtest.map)
enz.df.ratio.aovftest$labels <- with(enz.df.ratio.aovftest,
                                     paste("F = ", F, "P =", P))
enz.df.ratio.aovftest$xpos <- rep(1400, 12)
enz.df.ratio.aovftest$ypos <- rep(c(5.25, 1.5, 1.25), 4)
ggplot(df.enz.ratio, aes(x = MAP, y = value, color = fd)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = 'lm', color = "black", size = 0.4, se = FALSE) +
  facet_grid(variable ~ fd, scales = "free_y") +
  labs(x = "Mean annual precipitation (mm)", 
       y = "") +
  geom_text(data = enz.df.ratio.aovftest, 
            aes(x = xpos, y = ypos, label = labels, 
                group = NULL, hjust = 0.5, vjust = 1.0)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12)) 
# ggsave("./outputs/NEECF_enz_lat_depth_ratio_MAP.pdf", height = 4.5, width = 8)

#==========================================================
# Figure 3: PCA analysis
# select and transform the data for PCA
df.fig3 <- enz.neecf.clean %>%
  select(site, treat, fd, BG, CB, POX, PER, NAG, LAP, AP) %>%
  mutate(BG = log(BG + 1),
         CB = log(CB + 1),
         POX = log(POX + 1),
         PER = log(PER + 1),
         NAG = log(NAG + 1),
         LAP = log(LAP + 1), 
         AP = log(AP + 1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MyGgBiplot

MyggBiplot <- function(dat, pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                       varname.adjust = 2.5) {
  # extract scores of two PCs
  stopifnot(length(choices) == 2)
  df.u <- data.frame(pcobj$x[, choices], 
                     site = dat$site, 
                     treat = dat$treat)
  
  # calculate the scaled variables and observations 
  # code modified from biplot.prcomp by getAnywhere(biplot.prcomp)
  lam <- pcobj$sdev[choices]  # pc.biplot = TRUE
  n <- NROW(pcobj$x)
  lam <- lam * sqrt(n)
  if (scale < 0 || scale > 1)
    warning("'scale' is outside [0, 1]")
  if (scale != 0)
    lam <- lam^scale
  else lam <- 1
  if (pc.biplot)
    lam <- lam/sqrt(n)
  df.u[, choices] <- t(t(df.u[, choices])/lam)
  names(df.u)[1:2] <- c("x.var", "y.var")
  df.v <- t(t(pcobj$rotation[, choices])*lam)
  
  # calculate mean and se of the values of x and y
  # require("plyr")
  df.u <- ddply(df.u, .(site, treat),
                summarise,
                xvar = mean(x.var),
                xvar.se = mean(x.var)/sqrt(length(x.var)),
                yvar = mean(y.var),
                yvar.se = mean(y.var)/sqrt(length(y.var)))
  
  # caculate variable directions
  # code modified from biplot.default by getAnywhere(biplot.default)
  df.v <- df.v * 0.8
  df.v <- as.data.frame(df.v)
  names(df.v) <- c("xvar", "yvar")
  df.v$varname <- rownames(df.v)
  
  # Append the proportion of explained variance to the axis labels
  # code modified from ggbiplot
  u.axis.labs <- paste('PC', choices, sep='')
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[1:2]^2 / sum(pcobj$sdev^2)))
  # Variables for text label position
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust <- with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # site lable
  df.sit <- ddply(df.u, .(site, treat),
                  summarise,
                  xvar = mean(xvar),
                  yvar = mean(yvar))
  df.sit$hjust <- with(df.sit, (1 - varname.adjust * sign(xvar)) / 2)
  df.sit$vjust <- with(df.sit, (1 - varname.adjust * sign(yvar)) / 2)
  
  # plot the biplot with ggplot2
  # the main plot
  mainplot <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + 
    coord_equal() +
    geom_point(aes(shape = treat), size = 3.5) +
    scale_shape_manual(values = c(1, 19)) +
    geom_errorbar(aes(ymin = yvar - yvar.se, ymax = yvar + yvar.se), size = 0.5) +
    geom_errorbarh(aes(xmin = xvar - xvar.se, xmax = xvar + xvar.se), size = 0.5) +
    geom_text_repel(data = df.sit, 
                    aes(label = site, x = xvar, y = yvar), 
                    color = 'black', size = 3.5) +
    # geom_text(data = df.sit, 
    #           aes(label = site, x = xvar, y = yvar,
    #               hjust = hjust, vjust = vjust), 
    #           color = 'black', size = 2.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
    xlim(-3.2, 3.2) + ylim(-3.2, 3.2) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none',  # c(0.15, 0.90),
          legend.title = element_blank())
  
  # plot with arrows
  subplot <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) +
    coord_equal() +
    geom_segment(data = df.v,
                 aes(x = 0, y = 0, xend = xvar, yend = yvar),
                 arrow = arrow(length = unit(0.25, 'picas')), 
                 color = 'red') +
    geom_text_repel(data = df.v,  
                    aes(label = varname, x = xvar, y = yvar), 
                    color = 'black', size = 4.5) +
    # geom_text(data = df.v, 
    #           aes(label = varname, x = xvar, y = yvar, 
    #               angle = angle, hjust = hjust), 
    #           color = 'black', size = 2.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
    xlim(-1.2, 1.2) + ylim(-1.2, 1.2) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none',
          legend.title = element_blank())
  list(mainplot, subplot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA results

# by 0-10 cm soil depth
df.fig3A <- df.fig3[df.fig3$fd == "0-10cm", ]
df.fig3A <- na.omit(df.fig3A)
pca3A <- prcomp(df.fig3A[enz.vars], 
                center = TRUE, scale = TRUE)
p1 <- MyggBiplot(df.fig3A, pca3A, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[1]]
p5 <- MyggBiplot(df.fig3A, pca3A, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[2]]

# by 10-20 cm soil depth
df.fig3B <- df.fig3[df.fig3$fd == "10-20cm", ]
df.fig3B <- na.omit(df.fig3B)
pca3B <- prcomp(df.fig3B[enz.vars], 
                center = TRUE, scale = TRUE)
p2 <- MyggBiplot(df.fig3B, pca3B, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[1]]
p6 <- MyggBiplot(df.fig3B, pca3B, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[2]]

# by 20-40 cm soil depth
df.fig3C <- df.fig3[df.fig3$fd == "20-40cm", ]
df.fig3C <- na.omit(df.fig3C)
pca3C <- prcomp(df.fig3C[enz.vars], 
                center = TRUE, scale = TRUE)
p3 <- MyggBiplot(df.fig3C, pca3C, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[1]]
p7 <- MyggBiplot(df.fig3C, pca3C, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[2]]

# by 40-60 cm soil depth
df.fig3D <- df.fig3[df.fig3$fd == "40-60cm", ]
df.fig3D <- na.omit(df.fig3D)
pca3D <- prcomp(df.fig3D[enz.vars], 
                center = TRUE, scale = TRUE)
p4 <- MyggBiplot(df.fig3D, pca3D, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[1]]
p8 <- MyggBiplot(df.fig3D, pca3D, choices = 1:2, scale = 1, pc.biplot = TRUE,
                 varname.adjust = 1.5)[[2]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save the mainplot
p1 = arrangeGrob(p1, ncol=1, 
                 left = textGrob("(a)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p2 = arrangeGrob(p2, ncol=1, 
                 left = textGrob("(b)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p3 = arrangeGrob(p3, ncol=1, 
                 left = textGrob("(c)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p4 = arrangeGrob(p4, ncol=1, 
                 left = textGrob("(d)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))

# put the four graphs into one
final <- arrangeGrob(p1, p2, p3, p4)

# save the graph
# ggsave(filename = "./outputs/NEECF_enz_figure3_V1.pdf", 
#        plot = final, width = 10.27, height = 11.69)
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save the subplot
p5 = arrangeGrob(p5, ncol=1, 
                 left = textGrob("(a)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p6 = arrangeGrob(p6, ncol=1, 
                 left = textGrob("(b)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p7 = arrangeGrob(p7, ncol=1, 
                 left = textGrob("(c)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))
p8 = arrangeGrob(p8, ncol=1, 
                 left = textGrob("(d)", y = 0.9, vjust = 1, hjust = 0.5,
                                 gp = gpar(fontsize = 20)))

# put the four graphs into one
final <- arrangeGrob(p5, p6, p7, p8)

# save the graph
# ggsave(filename = "./outputs/NEECF_enz_figure3_V2.pdf", 
#        plot = final, width = 10.27, height = 11.69)

###########################################################
# END of the script                                       #
###########################################################