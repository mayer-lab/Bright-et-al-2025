## create source data for each figure ##

library(writexl)

## main figures:
df1 <- read.table("Results/source_data/fig1g.csv", sep = ",")
df2 <- read.table("Results/source_data/fig1h_left.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/fig1h_right.csv", sep = ",", h=T)
df4 <- read.table("Results/source_data/fig1m.csv", sep = ",", h=T)
df5 <- read.table("Results/source_data/fig1n.csv", sep = ",", h=T)

df0 <- data.frame(
  value = c(
    -65.3967,-60.3433,-64.9933,-68.4333,-68.5933,-66.2933,
    -57,-60.9,-66.09,-58.84,-61.1,-65.54,-61.5733,
    -69.7567,-57.5067,-69.3667,-56.18,-64.7767,-66.8833,-67.27,-66.3867,-67.6467,
    -46.83,-43.4367,-24.24,-42.6,-32.3,-39.9533,-35.41,
    -59.06,-67.1667,-66.5167,-64.04,-65.1267,-69.1433,-62.38,-70.23,-70.17,-69.58,
    -61.26,-63.2233,-62,
    -58.9533,-65.9067,
    -60.6167,-52.5733,-64.7333,-71.1767,-63.44,-72.2233,-69.66
  ),
  spatial = c(rep("LGE",6),rep("CGE",7),rep("MGE",9),rep("cortex",7), rep("LGE",10),rep("CGE",3),rep("MGE",2),rep("cortex",7)),
  stage = c(rep("e13.5",29), rep("e15.5",22))
)

df_list <- list(
  "Fig 1g" = df1,
  "Fig 1h-left" = df2,
  "Fig 1h-right" = df3,
  "Fig 1j" = df0,
  "Fig 1m" = df4,
  "Fig 1n" = df5
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/fig1_source_data.xlsx", format_headers = T)

## 2:
df1 <- read.table("Results/source_data/fig2b.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/fig2c.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/fig2d.csv", sep = ",", h=T)
df4 <- read.table("Results/source_data/fig2e.csv", sep = ",", h=T)
df5 <- read.table("Results/source_data/fig2g_upper.csv", sep = ",", h=T)
df6 <- read.table("Results/source_data/fig2g_bottom.csv", sep = ",", h=T)
df7 <- read.table("Results/source_data/fig2i.csv", sep = ",", h=T)
df8 <- read.table("Results/source_data/fig2k.csv", sep = ",", h=T)

df_list <- list(
  "Fig 2b" = df1,
  "Fig 2c" = df2,
  "Fig 2d" = df3,
  "Fig 2e" = df4,
  "Fig 2g-upper" = df5,
  "Fig 2g-bottom" = df6,
  "Fig 2i" = df7,
  "Fig 2k" = df8
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/fig2_source_data.xlsx", format_headers = T)

## 3:
df1 <- read.table("Results/source_data/fig3a_nodes.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/fig3a_edges.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/fig3b_nodes.csv", sep = ",", h=T)
df4 <- read.table("Results/source_data/fig3b_edges.csv", sep = ",", h=T)
df5 <- read.table("Results/source_data/fig3f.csv", sep = ",", h=T)

df_list <- list(
  "Fig 3a-nodes" = df1,
  "Fig 3a-edges" = df2,
  "Fig 3b-nodes" = df3,
  "Fig 3b-edges" = df4,
  "Fig 3f" = df5
)
write_xlsx(df_list, path = "Results/source_data/bundled/fig3_source_data.xlsx", format_headers = T)

## 4:
df1 <- read.table("Results/source_data/fig4b.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/fig4c.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/fig4h.csv", sep = ",", h=T)
df4 <- read.table("Results/source_data/fig4i.csv", sep = ",", h=T)
df5 <- read.table("Results/source_data/fig4j_upper.csv", sep = ",", h=T)
df6 <- read.table("Results/source_data/fig4j_lower.csv", sep = ",", h=T)
df7 <- read.table("Results/source_data/fig4k.csv", sep = ",", h=T)

df_list <- list(
  "Fig 4b" = df1,
  "Fig 4c" = df2,
  "Fig 4h" = df3,
  "Fig 4i" = df4,
  "Fig 4j-upper" = df5,
  "Fig 4j-bottom" = df6,
  "Fig 4k" = df7
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/fig4_source_data.xlsx", format_headers = T)

## Extended Figures ##

## EDF1:
df1 <- read.table("Results/source_data/edf1b.csv", sep = ",",h=T)

df_list <- list(
  "EDF 1b" = df1
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf1_source_data.xlsx", format_headers = T)


## EDF 2:
df1 <- read.table("Results/source_data/edf2a.csv", sep = ",",h=T)

df_list <- list(
  "EDF 2a" = df1
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf2_source_data.xlsx", format_headers = T)

## EDF 3:
df1 <- read.table("Results/source_data/edf3a.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/edf3e.csv", sep = ",", h=T)

df_list <- list(
  "EDF 3a" = df1,
  "EDF 3e" = df2
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf3_source_data.xlsx", format_headers = T)

## EDF 4:
df1 <- read.table("Results/source_data/edf4e.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/edf4f.csv", sep = ",", h=T)

df_list <- list(
  "EDF 4e" = df1,
  "EDF 4f" = df2
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf4_source_data.xlsx", format_headers = T)

## EDF 5:
df1 <- read.table("Results/source_data/fig3a_nodes.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/fig3a_edges.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/edf5b_nodes.csv", sep = ",",h=T)
df4 <- read.table("Results/source_data/edf5b_edges.csv", sep = ",", h=T)
df5 <- read.table("Results/source_data/edf5c_nodes.csv", sep = ",",h=T)
df6 <- read.table("Results/source_data/edf5c_edges.csv", sep = ",", h=T)
df7 <- read.table("Results/source_data/edf5d.csv", sep = ",",h=T)
df8 <- read.table("Results/source_data/edf5f_nodes.csv", sep = ",", h=T)
df9 <- read.table("Results/source_data/edf5f_edges.csv", sep = ",",h=T)
df10 <- read.table("Results/source_data/edf5g_nodes.csv", sep = ",", h=T)
df11 <- read.table("Results/source_data/edf5g_edges.csv", sep = ",",h=T)
df12 <- read.table("Results/source_data/edf5h.csv", sep = ",", h=T)

df_list <- list(
  "EDF 5a-nodes" = df1,
  "EDF 5a-edges" = df2,
  "EDF 5b-nodes" = df3,
  "EDF 5b-edges" = df4,
  "EDF 5c-nodes" = df5,
  "EDF 5c-edges" = df6,
  "EDF 5d" = df7,
  "EDF 5f-nodes" = df8,
  "EDF 5f-edges" = df9,
  "EDF 5g-nodes" = df10,
  "EDF 5g-edges" = df11,
  "EDF 5h" = df12
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf5_source_data.xlsx", format_headers = T)

## EDF 6:
df1 <- read.table("Results/source_data/edf6e_left.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/edf6e_right.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/edf6f.csv", sep = ",",h=T)
df4 <- read.table("Results/source_data/edf6g.csv", sep = ",", h=T)

df_list <- list(
  "EDF 6e-left" = df1,
  "EDF 6e-right" = df2,
  "EDF 6f" = df3,
  "EDF 6g" = df4
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf6_source_data.xlsx", format_headers = T)

## EDF 7:
df1 <- read.table("Results/source_data/edf7e.csv", sep = ",",h=T)
df2 <- read.table("Results/source_data/edf7f.csv", sep = ",", h=T)
df3 <- read.table("Results/source_data/edf7g.csv", sep = ",",h=T)
df4 <- read.table("Results/source_data/edf7h.csv", sep = ",", h=T)

df_list <- list(
  "EDF 7e" = df1,
  "EDF 7f" = df2,
  "EDF 7g" = df3,
  "EDF 7h" = df4
)
df_list <- lapply(df_list, function(el) {
  cbind(" " = rownames(el), el)
})
write_xlsx(df_list, path = "Results/source_data/bundled/edf7_source_data.xlsx", format_headers = T)
