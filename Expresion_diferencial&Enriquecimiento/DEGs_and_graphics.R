
#Limpiamos el entorno
rm(list = ls())

#Seleccionamos el directorio de trabajo
setwd("~/Máster en Bioinformática UNIR/Secuenciación y Ómicas de próxima generacion/Actividad 2/mubio03_act2")

#Cargamos las librerías necesarias
library("dplyr")
library("readr")
library("DESeq2")
library("ggplot2")
library("edgeR")
library("limma")
library("ggrepel")
library("pheatmap")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ReactomePA")
library("patchwork")

#Cargamos el archivo con la matriz de datos que contiene los read counts por transcrito y muestra
counts <- read.csv("matriz_counts_salmon.csv",
                   row.names = 1,
                   check.names = FALSE)

counts$transcript <- rownames(counts) # Esto añade una columna nueva con el nombre de cada transcrito

#Leemos el archivo que contiene la asociación entre transcrito y gen
map <- read_tsv("Transcrito_a_gen.tsv",
                col_names = c("transcript", "gene"))

#Añadimos la columna del gen correspondiente al lado de la del transcrito
df <- counts %>%
  left_join(map, by = "transcript")

#Filtramos transcritos sin gen (por si los hubiera)
df <- df %>% filter(!is.na(gene))

#Sumamos los read counts de cada transcrito y nos quedamos con el nombre de cada gen en lugar del transcrito
gene_counts <- df %>%
  dplyr::select(-transcript) %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()

#Preparamos la matriz para hacer el análisis por expresión diferencial mediante edgeR+limma (genes x individuo)
rownames(gene_counts) <- gene_counts$gene
gene_counts$gene <- NULL

#Quitamos todos los resultados que sean 0
gene_counts <- gene_counts[rowSums(gene_counts) > 0, ]

#Preparamos los metadatos para el análisis
coldata <- data.frame(
  row.names = c("AbrahamSimpson", "HomerSimpson",
                "BartSimpson", "LisaSimpson", "MaggieSimpson"),
  condition = c("Obeso", "Obeso", "Normopeso", "Normopeso", "Normopeso")
)

#Convertimos a factor para el análisis
coldata$condition <- factor(coldata$condition)

#Realizamos el análisis con edgeR
dge <- DGEList(counts = gene_counts, # Nuestra matriz
               group = coldata$condition) # Nuestros metadatos

#Filtramos para quitar posibles genes con baja expresión
keep <- filterByExpr(dge) # Determinamos que genes tienen suficientes cuentas para ser tenidos en cuenta en el análisis estadístico
dge <- dge[keep, , keep.lib.sizes = FALSE] # Los incluimos en el objeto dge

#Normalizamos 
dge <- calcNormFactors(dge)

#Transformación con limma
design <- model.matrix(~ coldata$condition)

v <- voom(dge, design, plot = TRUE)

#Ajustamos a modelo lineal
fit <- lmFit(v, design)
fit <- eBayes(fit)

#Resultados
res <- topTable(fit,
                coef = 2,
                number = Inf)

head(res)
write.csv(res, file = "Resultados EdgeR+Limma.csv")

#Preparamos los datos para representar el volcano plot
res$gene <- rownames(res)

res$significant <- ifelse(
  res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5,
  "sig",
  "no_sig"
) # Se consideran significativo aquellos genes cuyo valor de logFC > 0.5 y p.adj < 0.05

#Etiquetamos los genes más significativos
top_genes <- res[res$adj.P.Val < 0.05, ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 10)

#Representamos el volcano plot
ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 3.5) +
  scale_color_manual(values = c("no_sig" = "black",
                                "sig" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  labs(
    title = "Volcano plot: Obeso vs Normopeso",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

# Volcano plot

library(EnhancedVolcano)

EnhancedVolcano(res, 
                lab=res$gene, 
                x="logFC", 
                y="adj.P.Val", 
                labSize = 3, 
                axisLabSize = 10,
                title = "Volcano Plot",
                subtitle = "Obesos vs Normopeso",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 5.0,
                maxoverlapsConnectors = Inf
                )

#Generamos un objeto que almacene el nombre de los genes significativos
#Se consideran significativo aquellos genes cuyo valor de logFC > 0.5 y p.adj < 0.05
genes <- top_genes$gene

#Generamos la matriz que se va a emplear para representar el heatmap
mat <- gene_counts[rownames(gene_counts) %in% genes, ]

#Normalizamos para que sea comparable
mat_log <- log2(mat + 1)
gene_counts_log <- log2(gene_counts + 1)

#Escalamos por gen
mat_scaled <- t(scale(t(mat_log)))
gene_counts_scaled <- t(scale(t(gene_counts_log)))
gene_counts_scaled <- na.omit(gene_counts_scaled)

#Anotamos los grupos
annotation_col <- data.frame(
  Condition = coldata$condition
)
rownames(annotation_col) <- colnames(gene_counts_scaled)

#Representamos el heatmap
pheatmap(
  gene_counts_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Heatmap de genes sin estadístico"
)

#Preparamos el ranking para hacer un GSEA (filtramos por el valor t del resultado del análisis)
gene_list <- res$t
names(gene_list) <- res$gene
gene_list <- sort(gene_list, decreasing = TRUE) # Ordenamos descendente

#Convertimos genes a su valor ENTREZ ID
gene_df <- bitr(names(gene_list),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

#Alineamos ranking con ENTREZ
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

#Primero hacemos GSEA para enriquecimiento en KEGG
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

df_gsk <- gsea_kegg@result
df_gsk <- df_gsk[order(df_gsk$NES), ]
df_gsk$Description <- factor(df_gsk$Description, levels = df_gsk$Description)

p1 <- ggplot(df_gsk, aes(x = NES, y = Description)) +
      geom_point(aes(fill = NES),
                shape = 21,
                size = 4,
                color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      labs(y = NULL) +
      theme_bw() # No hay enriquecimiento significativo en valores KEGG

#Ahora hacemos GSEA para GO
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

df_gsg <- gsea_go@result

# Top 10 positivos (mayor NES)
top_pos <- df_gsg[df_gsg$NES > 0, ]
top_pos <- top_pos[order(top_pos$NES, decreasing = TRUE), ]
top_pos <- head(top_pos, 10)

# Top 10 negativos (menor NES)
top_neg <- df_gsg[df_gsg$NES < 0, ]
top_neg <- top_neg[order(top_neg$NES), ]
top_neg <- head(top_neg, 10)

# Unir ambos
df20 <- rbind(top_pos, top_neg)

# Orden para el gráfico
df20 <- df20[order(df20$NES), ]
df20$Description <- factor(df20$Description, levels = df20$Description)

p2 <- ggplot(df20, aes(x = NES, y = Description)) +
      geom_point(aes(fill = NES),
                shape = 21,
                size = 4,
                color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      labs(y = NULL) +
      theme_bw() # No hay enriquecimiento significativo en valores GO

#Por último, hacemos GSEA en Reactome
gsea_reactome <- gsePathway(
  geneList = gene_list,
  organism = "human",
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

df_gsr <- gsea_reactome@result
df_gsr <- df_gsr[order(df_gsr$NES), ]
df_gsr$Description <- factor(df_gsr$Description, levels = df_gsr$Description)

p3 <- ggplot(df_gsr, aes(x = NES, y = Description)) +
      geom_point(aes(fill = NES),
                shape = 21,
                size = 4,
                color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      labs(y = NULL) +
      theme_bw() # No hay enriquecimiento significativo en valores Reactome

#Juntamos todos los dotplots en uno
all_nes <- c(df_gsk$NES, df_gsg$NES, df_gsr$NES)
x_lim <- range(all_nes, na.rm = TRUE)

p1 <- p1 + coord_cartesian(xlim = x_lim)
p2 <- p2 + coord_cartesian(xlim = x_lim)
p3 <- p3 + coord_cartesian(xlim = x_lim)

(p1 / p2 / p3)
write.csv(top_genes, file = "Genes Significativos.csv")




library(glmGamPoi)

ego2 <- setReadable(
  gsea_go,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
  )

cnetplot(ego2)

