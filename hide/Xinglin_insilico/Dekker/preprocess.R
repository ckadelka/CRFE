go = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/Fazhir_GO.txt", header = F)
gsub( "\t.*$", "", go[1,1])
rn = c()
for (i in 1:nrow(go)) {
  rn[i] = gsub( "\t.*$", "", go[i,1])
}
go$rnames = rn

ws = c()
for (i in 1:nrow(go)) {
  ws[i] = sub(".*?\\t", "", go[i,1])
}
for (i in 1:nrow(go)) {
  ws[i] = sub(".*?\\t", "", ws[i])
}
for (i in 1:nrow(go)) {
  ws[i] = gsub("\t", "s", ws[i])
}

for (i in 1:nrow(go)) {
  ws[i] = gsub("s", "G ", ws[i])
}

go_to_save = data.frame(rn, ws)
write.table(go_to_save, file = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/Fazhir_GO_cleanG.txt", sep ="\t", row.names = F, quote = F, col.names = F)

genes = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/Fazhir.txt", header = F, sep = "\t")
genes$V1 = paste0(genes$V1, "G")
write.table(genes, file = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/Fazhir_G.txt", sep ="\t", row.names = F, quote = F, col.names = F)

genes_full = read.csv("/Users/JasonJia/Downloads/rankedlistcortiso_full.rnk", header = F, sep = "\t")
temp_g = as.vector(genes_full$V1[!genes_full$V1 %in% genes$V1])
temp_score = rep(0, length(temp_g))

g1 = c(genes$V1, temp_g)
g2 = c(genes$V2, temp_score)
genes1 = data.frame(g1,g2)
genes1$g1 = paste0(genes1$g1, "G")
write.table(genes1, file = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/Fazhir_G_full.txt", sep ="\t", row.names = F, quote = F, col.names = F)


genes = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/fazhir_real.txt", header = F, sep = "\t")
#genes$V1 = paste0(genes$V1, "G")
write.table(genes, file = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/Dekker/fazhir_real_G.txt", sep ="\t", row.names = F, quote = F, col.names = F)

