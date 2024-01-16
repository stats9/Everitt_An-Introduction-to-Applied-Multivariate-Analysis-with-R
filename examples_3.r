library(tidyverse)
library(smacof)
dat1 <- read.table(file = file.choose(), header = F)
dat1 %>% mutate(V1 = NULL) %>%
setNames(paste0("x", 1:5)) -> dat1
dat1
D_x <- dist(dat1)
Model_cmd <- cmdscale(dist(dat1), k = 5, add = TRUE)
Model2_smacof <- smacofSym(dist(dat1), ndim = 5)
Model2_smacof$conf -> v1
v1
colSums(v1)
colSums(Model_cmd)

Model_cmd$GOF


#################################

library(tidyverse)
library(MVA)
tapply(1:nrow(skulls), skulls$epoch, function(x){
dat = skulls[x, -1]
return(cov(dat))
}) -> skulls_var
id_row <- tapply(1:nrow(skulls), skulls$epoch, function(x){
dat = skulls[x, -1]
return(nrow(dat))
})
id_row
S <- matrix(0, 4, 4)
for(i in 1:5){
    S = S + (id_row[[i]] - 1) * skulls_var[[i]]
}
Sp = S / (nrow(skulls) - length(skulls_var)) 
Sp %>% round(3) -> Sp
Sp
xbar1 <- tapply(1:nrow(skulls), skulls$epoch, function(x){
colMeans(skulls[x, -1])
}) 
xbar1
class(xbar1)
dim(xbar1)
xbar1 %>% transpose %>% as_tibble() %>% unnest(cols = c(mb, bh, bl, nh)) %>% 
as.data.frame
