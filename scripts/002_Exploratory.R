
# Get the total and rare variants per gene --------------------------------
info.expanded$FLAG <- 1
info.expanded <- info.expanded %>% mutate(
  RARE = case_when(
    AF_CTRL <= 0.001 ~ 1,
    TRUE ~ 0
  )
)

info.expanded <- info.expanded %>% mutate(
  ROR = case_when(
    ESTIMATES >= 2 & PVAL <= 0.05 & AF_CTRL <= 0.05 ~ 1,
    TRUE ~ 0
  )
)

var.count <- aggregate(FLAG ~ SYMBOL, info.expanded, sum)
rare.count <- aggregate(RARE ~ SYMBOL, info.expanded, sum)
risk.count <- aggregate(ROR ~ SYMBOL, info.expanded, sum)

cohort.info <- NULL
cohort.info <- inner_join(var.count,rare.count,by = "SYMBOL")
cohort.info <- inner_join(cohort.info,risk.count,by = "SYMBOL")


gene.length <- vroom("data/gene.lengths.txt")
colnames(gene.length) <- c("SYMBOL","GID","TID","LENGTH")
gene.bak <- gene.length
gene.length <- aggregate(LENGTH ~ SYMBOL, gene.length, max)



cohort.info <- inner_join(cohort.info,gene.length,by = "SYMBOL")

colnames(cohort.info) <- c("SYMBOL","TOTAL","RARE","RISK","LENGTH")


# # Correlation -------------------------------------------------------------
ggplot(cohort.info, aes(x = LENGTH,y = TOTAL, group = 1)) +
  geom_point(aes()) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5) +
  ggtitle("Total vs LENGTH") +
  xlab("Length of gene") +
  ylab("Total variants per Gene") +
  stat_cor(method = "pearson", cor.coef.name = "r", vjust = 1, size = 4)

ggplot(cohort.info, aes(x = LENGTH,y = RARE, group = 1)) +
  geom_point(aes()) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5) +
  ggtitle("Rare vs LENGTH") +
  xlab("Length of gene") +
  ylab("Rare variants per Gene") +
  stat_cor(method = "pearson", cor.coef.name = "r", vjust = 1, size = 4)

ggplot(cohort.info, aes(x = LENGTH,y = RISK, group = 1)) +
  geom_point(aes()) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5) +
  ggtitle("Risk vs LENGTH") +
  xlab("Length of gene") +
  ylab("Risky variants per Gene") +
  stat_cor(method = "pearson", cor.coef.name = "r", vjust = 1, size = 4)


ggplot(cohort.info, aes(y = RISK,x = TOTAL, group = 1)) +
  geom_point(aes()) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5) +
  ggtitle("Risk vs TOTAL") +
  xlab("Total variants per gene") +
  ylab("Risky variants per Gene") +
  stat_cor(method = "pearson", cor.coef.name = "r", vjust = 1, size = 4)



# Model -------------------------------------------------------------------

model <- lm(cohort.info$RISK~cohort.info$TOTAL)
model

cohort.info <- cohort.info %>% mutate(
  EXP.RISKY = 0.02234*TOTAL + 0.08667
)




group <- c("OBSERVED","PREDICTED")

png("graph/Rare_Linear_Model.png", width = 84, height = 42, units = 'cm', res = 900)
ggplot(cohort.info, aes(x = TOTAL, group = 1)) +
  geom_point(aes(y = RISK, color = group[1])) +
  geom_point(aes(y = EXP.RISKY, color = group[2])) +
  geom_smooth(linetype='dotted',aes(x = TOTAL,y = RISK),method = "lm", se = TRUE, size = 0.5) +
  ggtitle("Observed vs Predicted Variants per Gene") +
  xlab("Total Variants per Gene") +
  ylab("Rare variants (OR > 2) per Gene") +
  theme_bw() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2),
        axis.text.y = element_text(size = rel(0.8))) +
  scale_colour_discrete("Rare Variants") +
  stat_cor(aes(y = RISK, x = TOTAL),method = "pearson", cor.coef.name = "r", vjust = 1, size = 4)
dev.off()

qualified.genes <- cohort.info %>% filter(RISK > EXP.RISKY) %>% select(SYMBOL) %>% unlist
