#Enter parameters here!
#Enter file location
FileLocation <- "~/Library/CloudStorage/OneDrive-InternationalMedicalUniversity/Research Data - Mutant rice.xlsx"
#Enter sheet name
SheetName <- "Set 2"
#Enter variable name
VariableName <- "Tproline"
#Plot labbeling
TitleName <- "Proline content in mutant rice lines under control and stress conditions"
xaxis <- "Mutant rice lines"
yaxis <- "Proline content (mg/g)"
#Define y-axis limit for sig. labels, [yes/no]
#ylim <- "yes"
ylimit <- 0.5

#You are free to run the code from here!

#Code starts here

library(ggplot2)
library(dplyr)
library(readxl)
library(data.table)
library(ggsignif)
library(ggh4x)

Dataset <- read_excel(FileLocation, sheet = SheetName)

Dataset$Condition <- factor(Dataset$Condition,
                            levels = c(1, 2),
                            labels=c('Control', 'Stress'))

summary.data.frame(Dataset)

#summarises mean and sd of proline content
Group_proline <- Dataset %>%
  group_by(Condition, Line) %>%
  summarise(proline_mean = mean(Tproline), proline_sd = sd(Tproline))

#ANOVA analysis p-value for control and stress condition
anova_result <- aov(Proline ~ Condition * Line, data = Dataset)
summary(anova_result)

# Post-hoc test using Tukey's HSD
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
plot(tukey_result) +
  with (par(mar=c(1,8,0,1) + 3.5),{plot(TukeyHSD(anova_result), las=1, cex.axis = 0.7, cex.lab = 1, cex.main = 2)})

# Plotting the bar graph
plot <- ggplot(Group_proline, aes(fill=Condition, y=proline_mean, x=Line), y=y, group = group) + 
  geom_bar(position="dodge", stat="identity") +
  expand_limits(y = ylimit) +
  geom_errorbar( aes(x=Line, ymin=proline_mean-proline_sd, ymax=proline_mean+proline_sd), 
                 width=0.4, colour="orange", linewidth=0.5, position = position_dodge(0.9)) +
  #geom_signif(y_position = c(0.475,0.45,0.425,0.4), xmin = c(1.25,2.25,3.25,3.75), 
  #            xmax = c(4.225,4.225,4.225,4.225), annotation = c("*","***","***","***"),
  #            tip_length = 0)+
  theme_minimal() + 
  labs(title = TitleName, x = xaxis, y = yaxis)

print(plot)
