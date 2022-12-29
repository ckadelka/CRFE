### Gen-model for CRFE v3
### put genes in P and U buckets based on beta_i and alpha_i
##Xinglin Jia
##10/24/2022
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
 
letsgen3 = function(data, alpha, beta, belif, n){ ## U < 0, 0 < P < 10 
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
  
  tp_g = c() ## true postie gene
  fn_g = c() ## false negative gene 
  C = c()
  
  for (j in 1:n) {
    term = data[[p[j]]] # pick a term
    no_dup = unique(setdiff(unlist(term), C)) # take out unique different genes in that term
    #print(length(no_dup))
    I = runif(length(no_dup), 0,1) # bernoulli identifier series
    no_dup_p = no_dup[I > beta_i[j]] # if bernulli>beta_i, perturbe, true postive gene
    #no_dup_up = no_dup[I < beta_i[j]]
    #print(length(no_dup_p))
    no_dup_up = no_dup[!no_dup %in% no_dup_p]# rest genes are false negative
    tp_g = c(tp_g, no_dup_p)
    fn_g = c(fn_g, no_dup_up)
    C = c(tp_g, fn_g)
  }
  
  ## uniformly assign score to tp genes
  tp_gs = runif(length(tp_g), 0, 1) ## true perturbed gene score
  
  ### now deal with unperturbed groups 
  up_g = unique(unlist(data))
  up_g = up_g[! up_g %in% C]
  I = runif(length(up_g), 0,1)
  fp_g = up_g[I<alpha] ## false postive gene
  tn_g = up_g[!up_g %in% fp_g] ## true negative gene
  
  up_g = c(fn_g, tn_g)
  up_gs = runif(length(up_g), -1, 0) ## unperturbed gene score
  
  #inters = 2*belif/(1+belif)
  #slope = 2-2*inters
  #cdf = function(x) -x*((b-1)*x-2*b)/(b+1)
  inv1 = function(x) (-belif+sqrt(x*(belif+1)*(-belif+1)+belif^2))/(-belif+1)
  #inv2 = function(x) -(b+sqrt(x*(b+1)*(-b+1)+b^2))/(-b+1)
  fp_gs = c()
  for (i in 1:length(fp_g)) {
    fp_gs[i] = inv1(runif(1,0,1))
  }
  
  p_g = c(tp_g, fp_g)
  p_gs = c(tp_gs, fp_gs)
  
  genes = c(p_g, up_g)
  score = c(p_gs, up_gs)
  
  gene_rank = data.frame(genes, score)
  gene_rank = gene_rank[order(gene_rank$score, decreasing = TRUE),]
  
  activated_terms = names(data)[p]
  
  up = data.frame(up_g, up_gs)
  up = up[order(up$up_gs, decreasing = TRUE),]
  P = data.frame(p_g, p_gs)
  P = P[order(P$p_gs, decreasing = TRUE),]
  
  betas = c(beta_i)
  alphas = c(alpha_i)
  
  result = list("activated_terms" = activated_terms, "gene_rank" = gene_rank, "perturbed_genes" = P,
                "unperturbed_genes" = up, "betas" = betas, "alphas" = alphas)
  return(result)
}

## testing 
set.seed(1)
result = letsgen3(data = human, alpha = .1, beta = 1/3, belif = 10, n =200)


## testing 
design = data.frame(c(rep(0.1, 8), rep(0.3, 8)), 
                    c(rep(0.3, 8), rep(0.1, 8)),
                    c(rep(c(rep(2, 4),rep(10, 4)),2)),
                    c(rep(c(5,10,50,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")

design = data.frame(c(rep(0.4, 1)), 
                    c(rep(0.25, 1)),
                    c(rep(c(5), 1)))
colnames(design) = c("alpha", "beta", "belif")

## ecoli

write.table(ecoli_anno, file = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli3/ecoli_anno.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:30) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen3(data = ecoli, alpha = d[1], beta = d[2], belif = d[3], n = d[4])
    write.table(result$gene_rank, file =paste0("data_small/ecoli3/ecoli", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]), "_", paste(d[4]), "_seed_",paste(j)  ,".txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("data_small/ecoli3/ecoli", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"_", paste(d[4]), "_seed_",paste(j), "key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}

## human
write.table(human_anno, file = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human4/human_anno.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:1000) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen3(data = human, alpha = d[1], beta = d[2], belif = d[3], n = 200)
    write.table(result$gene_rank, file =paste0("data_small/human3/human3/human", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]),"_", paste(d[4]), "_seed_",paste(j)  ,".txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("data_small/human3/human3/human", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"_", paste(d[4]), "_seed_",paste(j), "key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}

