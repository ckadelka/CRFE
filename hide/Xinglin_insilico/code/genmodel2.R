### Gen-model for CRFE v2
### put genes in P and U buckets based on beta_i and alpha then shuffle
##Xinglin Jia
##10/14/2022
##decision tree structure 

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

letsgen2 = function(data, alpha, beta, belif, n){
  beta_i = c()
  alpha_i = c()
  for (i in 1:n) {
    deno=0
    for (j in 1:n) {
      deno = deno + (n-j+ (j-1) * belif)
      #print(deno)
    }
    beta_i[i] = beta*(n-i + (i-1)*belif)*n/deno
    alpha_i[i] = alpha*(n-i + (i-1)*belif)*n/deno
  }
  # randomly select n terms
  p = sample(length(data), n)
  
  perturbed_gene = c()
  unperturbed_gene = c()
  for (j in 1:n) {
    term = data[[p[j]]]
    no_dup = unique(setdiff(unlist(term), perturbed_gene)) ### check 
    no_dup_p = sample(no_dup, (1-beta_i[j])*length(no_dup))
    no_dup_up = no_dup[!no_dup %in% no_dup_p]
    print(length(no_dup))
    print(length(no_dup_p))
    print(length(no_dup_up))
    print("next")
    perturbed_gene = c(perturbed_gene, no_dup_p)
    unperturbed_gene = c(unperturbed_gene, no_dup_up)
  }
  
  unperturbed_genes = unique(unlist(data) [! unlist(data) %in% perturbed_gene])
  fp_genes = sample(unperturbed_genes, alpha*length(unperturbed_genes))
  unperturbed_genes = unperturbed_genes[!unperturbed_genes %in% fp_genes]
  unperturbed_genes = unique(c(unperturbed_gene, unperturbed_genes))
  
  perturbed_gene = unique(c(perturbed_gene, fp_genes))
  unperturbed_genes = setdiff(unperturbed_genes, perturbed_gene)
  print("unperturbed_genes number")
  print(length(unperturbed_genes))
  print("unperturbed_gene number")
  print(length(unperturbed_gene))
  print("perturbed_genes number")
  print(length(perturbed_gene))
  unperturbed_gene_score = runif(length(unperturbed_genes), -10, -.001)
  perturbed_gene_score = runif(length(perturbed_gene), 0.001, 10)
  
  genes = c(perturbed_gene, unperturbed_genes)
  score = c(perturbed_gene_score, unperturbed_gene_score)
  gene_rank = data.frame(genes, score)
  gene_rank = gene_rank[order(gene_rank$score, decreasing = TRUE),]
  
  activated_terms = names(data)[p]
  
  up = data.frame(unperturbed_genes, unperturbed_gene_score)
  up = up[order(up$unperturbed_gene_score, decreasing = TRUE),]
  
  betas = c(beta_i)
  
  perturbed_genes = data.frame(perturbed_gene, perturbed_gene_score)
  perturbed_genes = perturbed_genes[order(perturbed_genes$perturbed_gene_score, decreasing = TRUE),]
  result = list("activated_terms" = activated_terms, "gene_rank" = gene_rank, "perturbed_genes" = perturbed_gene,
                "unperturbed_genes" = up, "betas" = betas)
  return(result)
}

## testing 
set.seed(1)
result = letsgen2(data = ecoli, alpha = .3, beta = .3, belif = 10, n =5)


## testing 
design = data.frame(c(rep(0.1, 3), rep(0.33, 3)), 
                    c(rep(0.33, 3), rep(0.25, 3)),
                    c(rep(c(2,5,10), 2)))
colnames(design) = c("alpha", "beta", "belif")

design = data.frame(c(.25),c(.25),c(5))
colnames(design) = c("alpha", "beta", "belif")

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

