## KI_Sym_Microbes Colors

# Clear working environment
rm(list=ls())

# Set colors for each disturbance level

distcols <- c("VeryLow"="#01665e",
              "Low"="#5ab4ac",
              "Medium"="#c7eae5",
              "High"="#d8b365",
              "VeryHigh"="#8c510a")

speccols <- c("#FF0000","#FFDB00","#49FF00","#00FF92","#0092FF","#4900FF","#FF00DB")


# Save
save.image("figures/KI_Sym_Microbes_colors.RData")
