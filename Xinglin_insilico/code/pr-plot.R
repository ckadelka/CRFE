library(plyr)
## PR plot
setwd("/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/output/ecoli/filtered_ecoli")
key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data/ecoli/"

design = data.frame(c(rep(0.05, 9), rep(0.3, 9)), 
                    c(rep(0.05, 18)),
                    c(rep(c(2,5,10), 6)),
                    c(rep(c(5,5,5,20,20,20,100,100,100), 2)))
colnames(design) = c("alpha", "beta", "belif","n_activated_terms")

for (j in 1:nrow(design)) {
  d = unlist(design[j,])
  filename = paste0("ecoli", paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_",paste(d[4]), ".txt")
  key_name = paste0("ecoli", paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_",paste(d[4]), "key.txt")
  data = read.csv(filename, sep = "\t")
  term = data[,8]
  key = read.csv(paste0(key_dir, key_name), header = F)
  # PR plot
  pr = c()
  re = c()
  p = 0
  for (i in 1:length(term)) {
    if (term[i] %in% unlist(key)){
      p = p +1}
    ppr = p/i
    rre = p/length(unlist(key))
    pr[i] = ppr
    re[i] = rre
  }
  plot(re, pr, main = filename)
}


## PR plot average
setwd("/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/output/ecoli/average/filtered_ecoli/filtered_ecoli_average")
key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data/ecoli/average"
design = data.frame(c(rep(0.05, 9), rep(0.3, 9)), 
                    c(rep(0.05, 18)),
                    c(rep(c(2,5,10), 6)),
                    c(rep(c(5,5,5,20,20,20,100,100,100), 2)))
colnames(design) = c("alpha", "beta", "belif","n_activated_terms")

## gen data for average PR
pr_lst = list()
re_lst = list()
for (i in c(7,8,9,16,17,18)){
  d = unlist(design[i,])
  for (seed in 1:100) {
    filename = paste0("ecoli", paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_",paste(d[4]), "seed_", seed, ".txt")
    key_name = paste0("ecoli", paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_",paste(d[4]),"seed_", seed,"key.txt")
    data = read.csv(filename, sep = "\t")
    term = data[,8]
    key = read.csv(paste0(key_dir, key_name), header = F)
    # PR plot
    pr = c()
    re = c()
    p = 0
    for (i in 1:length(term)) {
      if (term[i] %in% unlist(key)){
        p = p +1}
      ppr = p/i
      rre = p/length(unlist(key))
      pr[i] = ppr
      re[i] = rre
    }
    pr_lst[seed] = pr
    re_lst[seed] = re
    pr_mat = plyr::ldply(pr_lst, rbind)
    re_mat = plyr::ldply(re_lst, rbind)
    ## compute average
    pr_ave = c()
    re_ave = c()
    for (ind in 1:ncol(pr_mat)) {
      pr_ave[ind] = mean(pr_mat[,ind])
      re_ave[ind] = mean(re_mat[,ind])
    }
    plot(re_ave, pr_ave, main = paste0("ecoli", paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_",paste(d[4])_"average"))
  }
}
