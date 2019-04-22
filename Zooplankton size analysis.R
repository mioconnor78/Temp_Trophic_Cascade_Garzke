### Section S3.4: Zooplankton size analysis
```{r Size analysis}
sizes <- read.csv("../data/GarzkeSizes.csv")
#View(sizes)

sizes$week <- as.factor(sizes$week)
sizes$Tank <- as.factor(sizes$Tank)

sizes$taxon <- sizes$species
levels(sizes$taxon)[levels(sizes$taxon)== "calanoid"] <- "Copepod"
levels(sizes$taxon)[levels(sizes$taxon)== "cyclopoid"] <- "Copepod"
levels(sizes$taxon)[levels(sizes$taxon)== "Bosmina"] <- "Daphnia"
```

```{r}
data.S <- as_tibble(sizes) %>%
  group_by(Tank, trophic.level, week) %>%
  filter(trophic.level != 'P') %>%
  arrange(trophic.level, week) 
#dplyr::select(., Tank, trophic.level, week, invTT, average.temp, abundance.Daphnia, abundance.copepods, Daphnia.Copepod.Ratio, total.zoo.abundance.liter) %>%
```

`r fig_nums("Fig_S3.7")`

```{r Fig_S3.7 Zooplankton size data, fig.width = 3, fig.height = 3}  
hist(data.S$size)
hist(log(data.S$size))
```

```{r size model}
## tried with untransformed data, and truncated models would not converge. So trying again with log transformed data. 
# first try: lm random effects models with normally distributed errors
m2 <- lmer(log(size) ~ invTT + trophic.level + (1|Tank), data = data.S)
#shapiro.test(resid(m2)) # very non-normal (both raw and log trasnformed)

# poisson regression
m1 <- glm(size ~ trophic.level + invTT, family="poisson", data=data.S)
#shapiro.test(resid(m1)) # very non-normal; probably need a zero-truncated poisson

#qqnorm(resid(m2)) # but not too bad; go ahead with this.
#qqline(resid(m2))

#m3 <- glmmadmb(size ~ invTT + trophic.level + (1|Tank), family = "poisson", data = data.S)
#shapiro.test(resid(m3)) # very non-normal

#m4 <- glmmadmb(size ~ invTT + trophic.level + (1|Tank), family = "truncnbinom1", data = data.S)
#shapiro.test(resid(m4)) # very non-normal
```


```{r size model comparison, echo = FALSE}
m2a <- lmer(log(size) ~ invTT*trophic.level*taxon + (1|Tank), REML = FALSE, data = data.S)
#summary(m2a)
#shapiro.test(resid(m2a))

m2b <- lmer(log(size) ~ I(-(invTT-mean(invTT)))*taxon + I(-(invTT-mean(invTT)))*trophic.level + (1|Tank), REML = FALSE, data = data.S)
m2j <- lmer(log(size) ~ I(-(invTT-mean(invTT)))*taxon + trophic.level + (1|Tank), REML = FALSE, data = data.S)
m2c <- lmer(log(size) ~ I(-(invTT-mean(invTT))) + trophic.level*taxon + (1|Tank), REML = FALSE, data = data.S)
m2d <- lmer(log(size) ~ I(-(invTT-mean(invTT))) + trophic.level + (1|Tank), REML = FALSE, data = data.S)
m2k <- lmer(log(size) ~ I(-(invTT-mean(invTT)))*trophic.level + (1|Tank), REML = FALSE, data = data.S)
m2m <- lmer(log(size) ~ I(-(invTT-mean(invTT)))*trophic.level + taxon + (1|Tank), REML = FALSE, data = data.S)
m2e <- lmer(log(size) ~ I(-(invTT-mean(invTT))) + taxon + (1|Tank), REML = FALSE, data = data.S)
m2l <- lmer(log(size) ~ I(-(invTT-mean(invTT)))*taxon + (1|Tank), REML = FALSE, data = data.S)
m2f <- lmer(log(size) ~ I(-(invTT-mean(invTT))) + (1|Tank), REML = FALSE, data = data.S)
m2g <- lmer(log(size) ~ trophic.level*taxon + (1|Tank), REML = FALSE, data = data.S)
m2h <- lmer(log(size) ~ taxon + (1|Tank), REML = FALSE, data = data.S)
m2i <- lmer(log(size) ~ trophic.level + (1|Tank), REML = FALSE, data = data.S)

resS <- model.sel(m2b, m2c, m2d, m2e, m2f, m2g, m2h, m2i, m2j, m2k, m2l, m2m)
#summary(m2g)
#m.avg <- model.avg(m2g)
```

```{r}
S.plot2 <- ggplot(data = data.S, aes(x = -invTT, y = log(size), xmin = -40, xmax = -38.5, group = week, shape = as.factor(trophic.level), alpha = week)) +
  theme_bw() +
  geom_point(aes(group = as.character(week))) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        legend.key = element_rect(fill = NA), 
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  #scale_alpha("week", guide = "none") +
  scale_shape_manual(guide = guide_legend(title = "Trophic Treatment"), values = c(1,2,22,23,8,15,16,17)) +
  ylab("length (ln(cm))") + 
  scale_x_continuous("Temperature (1/kTi)", sec.axis =sec_axis(~((1/(k*-.))-273), name = xlab))

S.plot3 <- S.plot2 +
  geom_smooth(data = data.S, aes(x = -invTT, y = S(invTT)), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5)


```

\pagebreak

`r table_nums("Table_S3.8")`

```{r Table_S3.8: Size model comparison}
kable(resS, digits=2, col.names = c('Int','T~wj~','Tx','Z~j~','T~wj~*Tx','T~wj~*Z~j~','Tx*Z~j~', 'df', 'logLik', 'AICc', 'd', 'w'))
```


```{r zooplankton size coefficients}
Sd.PZ <- fixef(m2g)[1]
Sc.PZ <- fixef(m2g)[1] + fixef(m2g)[3]
Sd.PZN <- fixef(m2g)[1] + fixef(m2g)[2]
Sc.PZN <- fixef(m2g)[1] + fixef(m2g)[3] + fixef(m2g)[4]
confint(m2b) -> SizeInts

#df <- length(data.S$size)-(length(coef(m2g))+1)-1
#t.stat <- qt(0.975, df = df)
#D.l2 <- fixef(m2g)[1] - t.stat * sqrt(vcov(m2g)[1,1])
#D.u2 <- fixef(m2g)[1] + t.stat * sqrt(vcov(m2g)[1,1])
#C.l2 <- fixef(m2g)[1] + fixef(m2g)[3] - t.stat * sqrt(vcov(m2g)[1,1] + vcov(m2g)[3,3] + #2*vcov(m2g)[1,3])
#C.u2 <- fixef(m2g)[1] + fixef(m2g)[3] + t.stat * sqrt(vcov(m2g)[1,1] + vcov(m2g)[3,3] + 2*vcov(m2g)[1,3])

#S.est <- cbind(Estimate = c(fixef(m2g)[1], fixef(m2g)[1] + fixef(m2g)[3]), Lower = c(D.l2, C.l2), Upper = c(D.u2, C.u2))
S.est <- cbind(Daphnia = c(exp(Sd.PZ), exp(Sd.PZN)), Copepods = c(exp(Sc.PZ), exp(Sc.PZN)))
rownames(S.est) <- c("- Predators", "+ Predators")
```
`r table_nums("Table_S3.9")`

```{r Table_S3.9}
kable(S.est, digits = 2)
```



```{r size_plot}
S.plot <- ggplot(data = data.S, aes(x = -invTT, y = log(size), xmin = -40.2, xmax = -38.2)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(group = as.character(Tank), color = as.character(trophic.level))) +
  facet_grid(.~taxon) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        legend.text = element_text(size = 3), 
        legend.title = element_text(size = 3), 
        legend.key = element_rect(fill = NA), 
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  ylab("Zooplankton body size (log(mm))") + 
  scale_x_continuous("Mean Tank Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273), name = "Temperature C"))

S.plot2 <- S.plot +
  geom_ribbon(data = subset(data.S[(data.S$taxon == 'Daphnia'),]), aes(ymin = Sd.PZ-SizeInts[1], ymax = Sd.PZ+SizeInts[1]), fill = "gray", alpha = .5) +
  geom_smooth(data = subset(data.S[(data.S$taxon == 'Daphnia'),]), aes(x = -invTT, y = Sd.PZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'blue', size = 1) +
  geom_ribbon(data = subset(data.S[(data.S$taxon == 'Copepod'),]), aes(ymin = Sc.PZ-SizeInts[1], ymax = Sc.PZ+SizeInts[1]), fill = "gray", alpha = .5) +  
  geom_smooth(data = subset(data.S[(data.S$taxon == 'Copepod'),]), aes(x = -invTT, y = Sc.PZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'blue', size = 1.5)
```

`r fig_nums("Fig_S3.8")`

```{r size plot}
S.plot2
```