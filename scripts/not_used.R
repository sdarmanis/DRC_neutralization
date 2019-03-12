
# new dose levels as support for the line
newdata <- as.data.frame(seq(0.01,10000,10))
# predictions and confidence intervals
pm <- predict(list.reg[[i]], newdata=newdata, interval="confidence")
# new data with predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]
newdata$conc0 <-  newdata$`seq(0.01, 10000, 10)`
newdata$conc0log <- log10(newdata$conc0)
# Change 0 to something non-zero for log axis 
df <- as.data.frame(cbind(cor.means, conc))
df$conc0 <- df$conc
df$conc0[df$conc0 == 0] <- 0.001
df$conc0log <- log10(df$conc0)
# add stdev for each mean on df
df$lower <- cor.means-cor.st
df$upper <- cor.means+cor.st
# Plot and store plots 
list.reg.plots[[i]] <-  ggplot(df, aes(x = conc0, y = cor.means)) +
  geom_point() +
  coord_trans(x="log") + 
  # geom_errorbar(data=df, mapping=aes(x=conc0, ymin=upper, ymax=lower), width=0.2, size=1, color="blue") + 
  geom_line(data=newdata, aes(x=conc0, y=p, colour="red") ) +
  geom_ribbon(data=newdata, aes(x=conc0, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  theme_classic() + 
  coord_trans(x="log") + 
  xlab("Conc (ng/ml)") + ylab("%GFP") + ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + 
  scale_x_continuous(breaks=conc) + 
  ggtitle(abs.left[i])






