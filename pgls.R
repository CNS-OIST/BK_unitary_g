rm(list = ls())
wd = "./"
setwd(wd)
################################################## Import the tree
library(ape)
library(treeio)
iq.tree <- read.iqtree("iqtree_LG+C20+F+G_20260624044521-renamed.tree")
phy <- as.phylo(iq.tree)
plot(phy)
################################################## Import data
traits <- read.table(
  "pS_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
pgls.data <- traits[, c("Animal", "pS", "phylum")]
rownames(pgls.data) <- pgls.data$Animal
head(pgls.data)
################################################## Matching rows to tips
setdiff(rownames(pgls.data), phy$tip.label)
identical(rownames(pgls.data), phy$tip.label)
pgls.data <- pgls.data[phy$tip.label, , drop = FALSE]
head(pgls.data)
identical(rownames(pgls.data), phy$tip.label)
################################################## Matching weights to rows
### Compute "weight" argument for non-ultrametric tree
# https://blog.phytools.org/2012/04/using-nlmegls-for-phylogenetic.html
w <- diag(vcv.phylo(phy))
w
data.frame(
  data = rownames(pgls.data),
  tree = phy$tip.label,
  w = w
)
################################################## MODEL
### taking correlation structure with form as tip labels (book of Revell & Harmon, on p.70)
# corBM<-corBrownian(phy=primate.tree,form=~spp)
# corLambda<-corPagel(value=1,phy=primate.tree,form=~spp)
labels_row_order <- rownames(pgls.data)
corP<-ape::corPagel(value=1, phy=phy, fixed = FALSE, form = ~labels_row_order)
### PGLS
# model_1 <- nlme::gls(pS ~ 1, data = pgls.data, na.action = na.omit, correlation=corP, weights=nlme::varFixed(~w))
# This is just testing the intercept

### Model 0
model_0 <- nlme::gls(pS ~ phylum, data = pgls.data, na.action = na.omit,
                     method="REML")
summary(model_0)
### Model 1
model_1 <- nlme::gls(pS ~ phylum, data = pgls.data, na.action = na.omit,
                     correlation=corP,
                     method="REML")
summary(model_1)
### Model in paper
model_2 <- nlme::gls(pS ~ phylum, data = pgls.data, na.action = na.omit,
                        correlation=corP,
                        weights=nlme::varFixed(~w), method="REML")
summary(model_2)

