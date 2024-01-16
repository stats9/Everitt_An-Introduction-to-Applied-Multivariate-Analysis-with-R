Teen_sex <- matrix(c(21, 21, 14, 13, 8, 
8, 9, 6, 8, 2, 
2, 3, 4, 10, 10), 3, 5, byrow = T)
cols = c("16", "16-17", "17-18", "18-19", "19-20")
rows = c("No boyfriend", "Boyfriend no sex", "Boyfriend sex")
dimnames(Teen_sex) = list(rows, cols)
dat <- Teen_sex
get_chi_measure <- function(dat, Graph = TRUE){
rr <- nrow(dat)
cc <- ncol(dat)
dcol = matrix(0, cc, cc)
drow = matrix(0, rr, rr)
n <- sum(dat)
cd1 <- colSums(dat)
rd1 <- rowSums(dat)

## for columns
for(i in 1:cc){
    for(j in 1:cc){
        for(k in 1:rr){
            pki = dat[k, i]/cd1[i]
            pkj = dat[k, j] / cd1[j]
           dcol[i, j] = dcol[i, j] +  n/rd1[k] * (pki - pkj)^2
        }
    }
}

## for rows
for(i in 1:rr){
    for(j in 1:rr){
        for(k in 1:cc){
            pik = dat[i, k]/rd1[i]
            pjk = dat[j, k] / rd1[j]
           drow[i, j] = drow[i, j] +  n/cd1[k] * (pik - pjk)^2
        }
    }
}
if(Graph){
c1 <- cmdscale(dcol, eig = TRUE)
r1 <- cmdscale(drow, eig = TRUE)
c1$points -> Dim_cols
r1$points -> Dim_rows
tot <- rbind(Dim_rows, Dim_cols)
coord_res <- extendrange(as.vector(tot))
Labels = unlist(dimnames(dat))
D <- as.data.frame(tot) %>%
mutate(Labels = Labels) %>%
mutate(Id = factor(rep(c("Rows", "Cols"), c(nrow(dat), ncol(dat))))) %>%
setNames(c("Coordinate:1", "Coordinate:2", "Labels", "Id"))
D %>% 
ggplot(aes(x = `Coordinate:1`, y = `Coordinate:2`, color = Id)) + 
geom_label(aes(label = Labels),  alpha = .4) + theme_bw() + 
geom_hline(yintercept = 0, linetype = 2) + 
geom_vline(xintercept = 0, linetype = 2) +
coord_cartesian(xlim = coord_res, ylim = coord_res) -> P
print(P)
}
return(list(dcols = dcol, drows = drow, Plot = P))
}
res <- get_chi_measure(dat)



library(MASS)
?geom_label_repel
