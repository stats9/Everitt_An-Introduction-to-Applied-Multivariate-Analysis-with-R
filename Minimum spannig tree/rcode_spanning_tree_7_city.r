library(tidyverse)
library(corrplot)
library(sf)
library(igraph)


city_per <- c("شیراز", 
"اصفهان", "آمل", 
"بیرجند", "بجنورد", "مشهد", "سیرجان")
mydat <- data.frame(
    City = c("Shiraz", "Isfahan", "Amol", "Birjand", 
    "Bojnurd", "Mashhad", "Sirjan"), 
    Lat = c(29.591768, 32.661343, 36.471546, 
    32.83333, 37.4745, 36.310699, 29.45137), 
    Lon = c(52.583698, 51.680374, 52.355087, 59.221375, 
    57.32903, 59.599457, 55.6809)
)
mydat2 <- data.frame(
    City = city_per, 
    Lat = c(29.591768, 32.661343, 36.471546, 
    32.83333, 37.4745, 36.310699, 29.45137), 
    Lon = c(52.583698, 51.680374, 52.355087, 59.221375, 
    57.32903, 59.599457, 55.6809)
)



iran <- st_read("D:\\habib_repository_github\\iran_data_map\\irn_admbnda_adm1_unhcr_20190514.shp")

ggplot(iran) + geom_sf(fill = NA)  + 
theme_bw() -> map_iran
map_iran +
geom_text(data = mydat2, 
aes(x = Lon, y = Lat, label = City), color = "darkblue") -> P1
P1
dist_mat <- matrix(
    c(0, 489, 1110, 1073, 1344, 1354, 387, 489, 0, 609, 870, 
    1019, 1147, 627, 1110, 609, 0, 1045, 518, 802, 1147, 
    1073, 870, 1045, 0, 659, 492, 740, 1344, 1019, 518, 659, 0, 
    274, 1183, 1354, 1147, 802, 492, 274, 0, 1085, 
    387, 627, 1147, 740, 1183, 1085, 0
    ), 7, 7, dimnames = list(mydat$City, mydat$City)
)

corrplot(dist_mat, type = "full",
method = "number", 
is.corr = FALSE, col = c("green", "red", "darkblue"), 
number.digits = 0)

G <- make_empty_graph(n = 7)

V(G)$color <- "yellow"
V(G)$size <- 20
coord <- cbind(lon = mydat2$Lon, lat = mydat$Lat)
col1 <- c(rep(1, 6), rep(2, 5), 
rep(3, 4), rep(4, 3), rep(5, 2), rep(6, 1))
col2 <- c(2:7, 3:7, 4:7, 5:7, 6:7, 7)
col3 <- c(dist_mat[2:7, 1], 
dist_mat[3:7, 2], dist_mat[4:7, 3], 
dist_mat[5:7, 4], dist_mat[6:7, 5], dist_mat[7, 6])
names(col1) <- names(col2) <- names(col3) <- NULL
Edge <- numeric(42)
Edge[seq(1, 42, 2)] <- col1
Edge[seq(2, 42, 2)] <- col2

G %>% 
  add_edges(Edge) %>%
  set_edge_attr("color", value = "red") -> G
E(G)$weight <- E(G)$label <- col3
V(G)$label = city_per
G <- set_graph_attr(G, "layout", coord)
G <- set_edge_attr(G, "color", value = 1:21)
G <- set_edge_attr(G, "label.color", value = 1:21)
G_mst <- mst(G, algorithm = "prim")
plot.igraph(G, edge.arrow.mode = 3)
plot.igraph(G_mst, edge.arrow.mode = 3)








