# gen model for CRFE testing using normal distribution 
# Xinglin Jia
# 7/8/2022 

setwd("/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/")

d = scan("GOecoli.txt", what = "", sep = "\n")
ecoli = strsplit(d, "[[:space:]]+")
names(ecoli) <- sapply(ecoli, `[[`, 1)
ecoli <- lapply(ecoli, `[`, -1)
ecoli_anno = d[sapply(ecoli, length) %in% c(20:200)]
ecoli = ecoli[sapply(ecoli, length) %in% c(20:200)]



d = scan("GOhuman.txt", what = "", sep = "\n")
human = strsplit(d, "[[:space:]]+")
names(human) <- sapply(human, `[[`, 1)
human <- lapply(human, `[`, -1)
human_anno = d[sapply(human, length) %in% c(20:200)]
human = human[sapply(human, length) %in% c(20:200)]

#lapply(human, cat, "\n", file="ecoli_anno.txt", append=TRUE)
#write.table(names(human), file = "ecolinames.txt", quote = F, col.names = F, row.names = F)


letsgen = function(data, alpha, beta, belif, n){
  # data: list of terms, each term is a sub-list with genes 
  # false_pos_rate: false positive for unperturbed genes
  # false_neg_rate: false negative for perturbed genes
  # belif: how much false negative rate grow
  # n: number of perturbed terms
  erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
  
  std = 1
  anchor = 0 # anchor point
  beta_i = c() ## beta_i
  #alpha_i = c() ## alpha_i
  means = c()
  for (i in 1:n) {
    deno=0
    for (j in 1:n) {
      deno = deno + (n-j+ (j-1) * belif)
      #print(deno)
    }
    beta_i[i] = beta*(n-i + (i-1)*belif)*n/deno
    #alpha_i[i] = alpha*(n-i + (i-1)*belif)*n/deno
    means[i] = -sqrt(2)*erfinv(2*beta_i[i]-1)
  }
  
  mean_up = -sqrt(2)*erfinv(2*(1-alpha)-1)
  
  # randomly select n terms
  p = sample(length(data), n)
  
  perturbed_gene = c()
  perturbed_gene_score = c()

  for (j in 1:n) {
    term = data[[p[j]]]
    clean_no_dup = setdiff(unlist(term), perturbed_gene)
    perturbed_gene = c(perturbed_gene, clean_no_dup)
    perturbed_gene_score = c(perturbed_gene_score, rnorm(length(clean_no_dup), mean = means[j], sd = 1))
  }
  
  unperturbed_genes = unique(unlist(data) [! unlist(data) %in% perturbed_gene])
  unperturbed_gene_score = rnorm(length(unperturbed_genes), mean = mean_up, sd = 1)
  
  genes = c(perturbed_gene, unperturbed_genes)
  score = c(perturbed_gene_score, unperturbed_gene_score)
  gene_rank = data.frame(genes, score)
  gene_rank = gene_rank[order(gene_rank$score, decreasing = TRUE),]
  
  activated_terms = names(data)[p]
  
  perturbed_genes = data.frame(perturbed_gene, perturbed_gene_score)
  perturbed_genes = perturbed_genes[order(perturbed_genes$perturbed_gene_score, decreasing = TRUE),]
  
  up = data.frame(unperturbed_genes, unperturbed_gene_score)
  up = up[order(up$unperturbed_gene_score, decreasing = TRUE),]
  
  betas = c(beta_i)
  meanss = c(mean_up, means)
  
  result = list("activated_terms" = activated_terms, "gene_rank" = gene_rank, "perturbed_genes" = perturbed_gene,
                "unperturbed_genes" = up, "anchor" = anchor, "betas" = betas, "means" = meanss)
  return(result)
}

result = letsgen(data = ecoli, alpha = .1, beta = .1, belif = 10, n =5)


## testing 
design = data.frame(c(rep(0.1, 3), rep(0.4, 3)), 
                    c(rep(0.4, 3), rep(0.25, 3)),
                    c(rep(c(2,5,10), 2)))
colnames(design) = c("alpha", "beta", "belif")

design = data.frame(c(.25),c(.25),c(5))
colnames(design) = c("alpha", "beta", "belif")

## gen data for average PR

## ecoli

write.table(ecoli_anno, file = "data_small/ecoli/ecoli_anno.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:30) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen(data = ecoli, alpha = d[1], beta = d[2], belif = d[3], n = 30)
    write.table(result$gene_rank, file =paste0("data_small/ecoli/30ecoli", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]),"seed_",paste(j)  ,".txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("data_small/ecoli/30ecoli", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"seed_",paste(j), "key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}


## human
write.table(human_anno, file = "data_small/human/human_anno.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:30) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen(data = human, alpha = d[1], beta = d[2], belif = d[3], n = 100)
    write.table(result$gene_rank, file =paste0("data_small/human/human", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]),"seed_",paste(j),"n100.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("data_small/human/human", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"seed_",paste(j), "n20key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}

