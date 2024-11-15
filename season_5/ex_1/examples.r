library(tidyverse)
library(MVA)
colName <- c("country", "m0", "m25", "m50", "m75", 
"w0", "w25", "w50", "w75")
dat <- read_table(file = file.choose(), col_names = colName, 
col_types = cols("c", "d", "d", "d", "d", "d", "d", "d", "d"))
print(dat, n = 31)
################################
n <- nrow(dat)
sapply(1:3, function(f) factanal(dat[, -1], factors = f, 
method = "mle")$PVAL)
fact_model1 <- factanal(dat[, -1], factors = 3, method = "mle")
fact_model1
dat %>%
column_to_rownames(var = "country") %>%
factanal(x = ., factors = 3, method = "mle", 
scores = "regression") %>% .$scores %>% round(6) -> scores
print(scores)


library(ggrepel)
scores %>% as.data.frame  %>% 
mutate(country = dat$country) -> dat_score
dat_score %>%
ggplot() +
geom_label_repel(aes(x = Factor1, y = Factor2, label = country), 
fill = "yellow") + 
coord_cartesian(xlim = extendrange(scores[, 1]), 
ylim = extendrange(scores[, 2])) + 
theme_bw()

dat_score %>%
ggplot() +
geom_label_repel(aes(x = Factor1, y = Factor3, label = country), 
fill = "yellow") + 
coord_cartesian(xlim = extendrange(scores[, 1]), 
ylim = extendrange(scores[, 2])) + 
theme_bw()


dat_score %>%
ggplot() +
geom_label_repel(aes(x = Factor2, y = Factor3, label = country), 
fill = "yellow") + 
coord_cartesian(xlim = extendrange(scores[, 1]), 
ylim = extendrange(scores[, 2])) + 
theme_bw()

library(psych)
library(GPArotation)
fact_model3 <- fa(dat[, -1], nfactors = 3, 
rotate = "varimax")
fa.diagram(fact_model3)

fact_model4 <- fa(dat[, -1], nfactors = 3, 
rotate = "oblimin")

fa.diagram(fact_model4)


####################################

## ex2

colName2 <- read_table(file = file.choose(), 
col_names = F, col_type = cols("c")) %>% .$X1

cordat <- read_table(file = file.choose(), 
col_names = F) %>% setNames(colName2) %>% 
mutate(drug = colName2) %>%
column_to_rownames(var = "drug") %>%
mutate(drug = NULL) %>% .[, nrow(.):1]

f <- function(x) x * 1/100
library(corrplot)
library(contigencyTable2)
cordat %>% list_to_dataframe %>% 
mutate(nname = names(.)) %>%
column_to_rownames(var = "nname") %>%
mutate(nname = NULL) %>%
data.matrix -> cor_mat


cor_mat %>% f %>%
corrplot(method = "ellipse", type = "full", is.corr = T, 
col.lim = c(-1, 1), p.mat = cor_mat, insig = "p-value"
)

#########################
1/100 * cor_mat[nrow(cor_mat):1, ] -> cor_mat2
rownames(cor_mat2) <- rev(rownames(cor_mat2))

sapply(1:6, function(nf){
factanal(covmat = cor_mat2, factors = nf,
method = "mle", n.obs = 1634)$PVAL})

factanal(covmat = cor_mat2, factors = 6,
method = "mle", n.obs = 1634)

pfun <- function(nf) {
model_fa <- fa(r = cor_mat2, nfactors = nf,
rotate = "varimax", n.obs = 1634)
return(model_fa$residual)
}



pfun(6)

