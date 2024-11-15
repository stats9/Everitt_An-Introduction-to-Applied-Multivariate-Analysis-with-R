library(tidyverse)
library(sem)
library(corrplot)

input3 <- textConnection("
1.00
-0.04 1.00
0.61 -0.07 1.00
0.45 -0.12 0.59 1.00
0.03 0.49 0.03 -0.08 1.00
-0.29 0.43 -0.13 -0.21 0.47 1.00
-0.30 0.30 -0.24 -0.19 0.41 0.63 1.00
0.45 -0.31 0.59 0.63 -0.14 -0.13 -0.26 1.00
0.30 -0.17 0.32 0.37 -0.24 -0.15 -0.29 0.40 1.00
")
dat_3 <- readLines(input3)
n <- 123
corr_mat_2 <- matrix(NA, 9, 9)
j <- 1
for(i in 2:10){
    dat_3[i] %>% strsplit(" ") %>% unlist %>% as.numeric-> temp
    corr_mat_2[j, 1:(i-1)] <- corr_mat_2[1:(i-1), j] <- temp
    j <- j + 1
}

name1 <- paste0("statement", 1:9)
dimnames(corr_mat_2) <- list(name1, name1)
corr_mat_2

corr_mat_2 %>%
corrplot(method = "ellipse", type = "full", is.corr = T, 
col.lim = c(-1, 1), p.mat = corr_mat_2, insig = "p-value", 
sig.level = -.4
)


patient_model <- specifyModel(file = file.choose())

cfa_patient <- sem(patient_model, S = corr_mat_2, N = 123)

summary(cfa_patient)
attributes(cfa_patient)

pathDiagram(cfa_patient)
capture.output(pathDiagram(cfa_patient), file = "sem.dot")
install.packages("semPlot")
library(semPlot)
semPaths(cfa_patient, intercept = TRUE, whatLabel = "est",
         residuals = TRUE, exoCov = TRUE)

pathDiagram(cfa_patient, edge.labels = "values", 
ignore.double = FALSE, output.type = "html")

