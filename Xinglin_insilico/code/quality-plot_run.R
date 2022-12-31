library(plyr)


## load ontology data
setwd("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/")

d = scan("GOecoli.txt", what = "", sep = "\n")
ecoli = strsplit(d, "[[:space:]]+")
names(ecoli) <- sapply(ecoli, `[[`, 1)
ecoli <- lapply(ecoli, `[`, -1)
ecoli = ecoli[sapply(ecoli, length) %in% c(20:200)]

d = scan("GOhuman.txt", what = "", sep = "\n")
human = strsplit(d, "[[:space:]]+")
names(human) <- sapply(human, `[[`, 1)
human <- lapply(human, `[`, -1)
human_anno = d[sapply(human, length) %in% c(20:200)]
human = human[sapply(human, length) %in% c(20:200)]

par(mfrow=c(4,4)) 
### test for belief
design = data.frame(c(0.4),c(0.25),c(2),c(10))
colnames(design) = c("alpha", "beta", "belif","n")

design = data.frame(c(rep(0.1, 8), rep(0.3, 8)), 
                    c(rep(0.3, 8), rep(0.1, 8)),
                    c(rep(c(rep(2, 4),rep(10, 4)),2)),
                    c(rep(c(5,10,50,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")
## testing 
design = data.frame(c(rep(0.25, 4), rep(0.4, 4)), 
                    c(rep(0.4, 4), rep(0.25, 4)),
                    c(rep(c(rep(2, 2),rep(10, 2)),2)),
                    c(rep(c(10,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")

par(mfrow=c(1,2)) 
design = data.frame(.4,.25,2,10)
colnames(design) = c("alpha", "beta", "belif", "n")

dd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/human4_20_50/"
kd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/human4_20_50/"
seeds = 1 ## number of replicates
tp = 50

for (k in 1:nrow(design)) {
  des = design[k,]
  b1 = preplot_v1(data_dir = paste0(dd,"b1_n10rep/temp"),
                  key_dir = kd,
                  anno = human,
                  top_percent = tp,
                  design = des,
                  e_or_h = "human",
                  n_points = 1000,
                  m = seeds)
  
  b2 = preplot_v1(data_dir = paste0(dd,"b2_n10rep/temp"),
                  key_dir = kd,
                  anno = human,
                  top_percent = tp,
                  design = des,
                  e_or_h = "human",
                  n_points = 1000,
                  m = seeds)
  
  b5 = preplot_v1(data_dir = paste0(dd,"b5_n10rep/temp"),
                  key_dir = kd,
                  anno = human,
                  top_percent = tp,
                  design = des,
                  e_or_h = "human",
                  n_points = 1000,
                  m = seeds)
  
  b10 = preplot_v1(data_dir = paste0(dd,"b10_n10rep/temp"),
                   key_dir = kd,
                   anno = human,
                   top_percent = tp,
                   design = des,
                   e_or_h = "human",
                   n_points = 1000,
                   m = seeds)
  
  ## truth alpha = 0.4, beta = 0.25, b =10
  x = seq(0,0.999,1/1000)
  color = c("red", "blue", "black","purple", "green")
  plot(x, c(1,b1[2:1000]+3), pch = " ", main = paste0("human 30 rep, 10 repeats, genes = 20-50, α=",paste(des[1])," β=",paste(des[2]), " b=",paste(des[3])," n=", paste(des[4])),
       xlab = "Proportion of Top 50% Perturbed Genes that are Explained", ylab = "Quality score")
  
  lines(x, b1, pch = ".", col = color[1])
  lines(x, b2, pch = ".", col = color[2])
  lines(x, b5, pch = ".", col = color[3])
  lines(x, b10, pch = ".", col = color[4])
  
  #legend(0.7, 8, legend=c("CRFE b=1","CRFE b=2", "CRFE b=5", "CRFE b=10"), fill = color[1:4])
}



## run real data
setwd("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/")
d = scan("GOhuman_ft_named.txt", what = "", sep = "\n")
human = strsplit(d,"\t")
names(human) <- sapply(human, `[[`,1)
human <- lapply(human, `[`, -1)
for (i in 1:length(human)) {
  human[[i]] = strsplit(human[[i]], "[[:space:]]+")
}

dd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/"
gd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/adeno_symbol.txt"

b1 = preplot_real_overexpress(crfe_dir = paste0(dd,"b1/temp/adeno_symbol.txt"),
                              data_dir = gd,
                              anno = human,
                              top_percent = 50,
                              n_points = 1000,
                              threshold = 2)

b2 = preplot_real_overexpress(crfe_dir = paste0(dd,"b2/temp/adeno_symbol.txt"),
                              data_dir = gd,
                              anno = human,
                              top_percent = 50,
                              n_points = 1000,
                              threshold = 2)

b5 = preplot_real_overexpress(crfe_dir = paste0(dd,"b5/temp/adeno_symbol.txt"),
                              data_dir = gd,
                              anno = human,
                              top_percent = 50,
                              n_points = 1000,
                              threshold = 2)
b10 = preplot_real_overexpress(crfe_dir = paste0(dd,"b10/temp/adeno_symbol.txt"),
                              data_dir = gd,
                              anno = human,
                              top_percent = 50,
                              n_points = 1000,
                              threshold = 2)


color = c("red", "blue", "black","purple", "green")
plot(b5$percent_explained, seq((length(b5$percent_explained)+1)/length(b5$percent_explained),max(c(b5$quality, b1$quality, b2$quality, b10$quality)), (max(c(b5$quality, b1$quality, b2$quality, b10$quality))-1)/length(b5$percent_explained)), 
     pch = " ", main = "adeno_symbol, overexpress, threshold = 2",
     xlab = "Proportion of Top 50% Perturbed Genes that are Explained", ylab = "Quality score", xlim = c(0,1))

lines(b1$percent_explained, b1$quality, pch = ".", col = color[1])
lines(b2$percent_explained, b2$quality, pch = ".", col = color[2])
lines(b5$percent_explained, b5$quality, pch = ".", col = color[3])
lines(b10$percent_explained, b10$quality, pch = ".", col = color[4])

legend(0.4, 4, legend=c("CRFE b=1","CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:4]
)


#run qr
design = data.frame(c(rep(0.1, 32), rep(0.3, 32)), 
                    c(rep(0.3, 32), rep(0.1, 32)),
                    c(rep(c(rep(2, 16),rep(10, 16)),2)),
                    c(rep(c(rep("b1",4),rep("b2",4),rep("b5",4),rep("b10",4)),4)),
                    c(rep(c(5,10,50,100),16)))
colnames(design) = c("alpha", "beta", "belif", "n","b")
dd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/human4/"
kd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human4/"
seeds = 30 ## number of replicates

rec_result = c()
pre_result = c()
for (k in 1:nrow(design)) {
  des = design[k,]
  result = qr(data_dir = paste0(dd,paste(des[4]),"/cut"),
                  key_dir = kd,
                  anno = human,
                  design = des,
                  e_or_h = "human",
                  m = seeds)
  rec_result = rbind(rec_result, result$rec)
  pre_result = rbind(pre_result, result$pre)
}

par(mfrow=c(4,4)) 
color = c("red", "blue", "black","purple")
k=1
for (k in c(1,2,3,4,17,18,19,20,33,34,35,36,49,50,51,52)){
  plot(c(0,1), c(0,1), pch = " ", main = paste0("human α=",paste(design[k,1])," β=",paste(design[k,2]), " b=",paste(design[k,3])," n=", paste(design[k,5])),
       xlab = "Recall", ylab = "Precision")
  j = k
  for (i in 1:4) {
    lines(rec_result[j+4*(i-1),], pre_result[j+4*(i-1),], pch = ".", col = color[i])
  }
}



