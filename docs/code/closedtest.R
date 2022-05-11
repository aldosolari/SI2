rm(list=ls())

#===================
# CLOSED TESTING
#===================

#-----------------------------------
# Comparing three groups
#-----------------------------------

require("dobson")
require("tidyr") 
require("DescTools") 
require("ggplot2") 
require("ggpubr") 
require("multcomp") 
require("multtest") 

fit <- lm(weight ~ group -1, plant.dried)
fit$coefficients

summary(fit)$sigma^2

fit_12 <- lm(weight ~ I(group=="TreatmentB"), plant.dried)
p_12 <- anova(fit_12,fit)[2,"Pr(>F)"]
fit_13 <- lm(weight ~ I(group=="TreatmentA"), plant.dried)
p_13 <- anova(fit_13,fit)[2,"Pr(>F)"]
fit_23 <- lm(weight ~ I(group=="Control"), plant.dried)
p_23 <- anova(fit_23,fit)[2,"Pr(>F)"]

p_tuk <- TukeyHSD(aov(weight ~ group, data = plant.dried))$group[,'p adj']
p_dun <- DunnettTest(weight ~ factor(group), plant.dried, control="Control")$Control[,'pval']

fit_123 <- lm(weight ~ 1, plant.dried)
p_A <- anova(fit_123,fit)[2,"Pr(>F)"]
p_B <- min(p_tuk)
p_C <- min(p_dun)
p_D <- p_12

adjp <- matrix(NA,nrow=4,ncol=4,
               dimnames = list(c("A", "B", "C", "D"),
                               c("H12", "H13", "H23", "H123")))
adjp["A",] <- pmax(c(p_12,p_13,p_23,p_A),p_A)
adjp["B",] <- pmax(c(p_12,p_13,p_23,p_B),p_B)
adjp["C",] <- pmax(c(p_12,p_13,p_23,p_C),p_C)
adjp["D",] <- pmax(c(p_12,p_13,p_23,p_D),p_D)
adjp

#pdf("Figure_plant_bp.pdf")
bp <- ggboxplot(plant.dried, x = "group", y = "weight")
stats <- compare_means(weight ~ group, data = plant.dried, method = "t.test")
stats$p.adj <- round(adjp["A",-4],3)
bp + stat_compare_means(method = "anova", label.y = 10) + 
  stat_pvalue_manual(stats, label = "p = {p.adj}", y.position =   c(9, 8, 7))
#dev.off()