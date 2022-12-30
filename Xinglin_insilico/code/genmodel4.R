### Gen-model for CRFE v4
### v3: put genes in P and U buckets based on beta_i and alpha_i
### v4: slopes for true postive terms
##Xinglin Jia
##10/28/2022
##decision tree structure 
## load ontology data
setwd("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico")

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
human_anno = d[sapply(human, length) %in% c(20:50)]
human = human[sapply(human, length) %in% c(20:50)]

letsgen4 = function(data, alpha, beta, belif, n){ ## U < 0, 0 < P < 10 
  inters = 2*belif/(1+belif)
  slope = 2-2*inters
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
  tp_gs = c()
  
  tp_uppers = seq(0.01,1,0.01)
  tp_lowers = tp_uppers - 1/100
  
  tp_n = c() ##debug
  
  b1 = 4*belif/(1+belif)-2
  b0 = 1-b1/2
  ratio = b1*2/n
  
  for (j in 1:n) {
    term = data[[p[j]]] # pick a term
    no_dup = unique(setdiff(unlist(term), C)) # take out unique different genes in that term
    #print(length(no_dup))
    I = runif(length(no_dup), 0,1) # bernoulli identifier series
    no_dup_p = no_dup[I > beta_i[j]] # if bernulli>beta_i, perturbe, true postive gene
    no_dup_p_gs = c()
    
    density_set =c()
    for (ct in 1:100) {
      density_set[ct] = ct/100*b1+b0
    }
    density_set = density_set/sum(density_set)
    
    b1 = b1-ratio
    b0 = 1-b1/2
    
    density_in_line = c()
    sum_temp = 0
    for (i in 1:length(density_set)) {
      sum_temp = sum(sum_temp+density_set[i])
      density_in_line[i] = sum_temp
    }
    
    for (kct in 1:length(no_dup_p)) {
      loc = sum(runif(1,0,1)>(density_in_line))
      no_dup_p_gs[kct] = runif(1, tp_lowers[loc+1], tp_uppers[loc+1])
    }
    
    tp_n = c(tp_n, length(no_dup_p))
    #no_dup_up = no_dup[I < beta_i[j]]
    #print(length(no_dup_p))
    no_dup_up = no_dup[!no_dup %in% no_dup_p]# rest genes are false negative
    tp_g = c(tp_g, no_dup_p)
    tp_gs = c(tp_gs, no_dup_p_gs)
    fn_g = c(fn_g, no_dup_up)
    C = c(tp_g, fn_g)
  }
  
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
  fp_gs = c()
  for (i in 1:length(fp_g)) {
    fp_gs[i] = inv1(runif(1,0,1))
  }
  
  tp_gs = tp_gs[1:length(tp_g)]
  
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
                "unperturbed_genes" = up, "betas" = betas, "alphas" = alphas, "tp_n" = tp_n, "tp_gs" = tp_gs)
  return(result)
}

## testing 
set.seed(1)
result = letsgen4(data = human, alpha = .1, beta = 1/3, belif = 2, n =10)

result$tp_n
result$tp_gs
## testing 
design = data.frame(c(rep(0.1, 8), rep(0.3, 8)), 
                    c(rep(0.3, 8), rep(0.1, 8)),
                    c(rep(c(rep(2, 4),rep(10, 4)),2)),
                    c(rep(c(5,10,50,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")

## ecoli

write.table(ecoli_anno, file = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli4/ecoli_anno.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:30) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen4(data = ecoli, alpha = d[1], beta = d[2], belif = d[3], n = 5)
    write.table(result$gene_rank, file =paste0("data_small/ecoli4/ecoli", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]),"_seed_",paste(j),".txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("data_small/ecoli4/ecoli", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"_seed_",paste(j), "key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}

## testing 
design = data.frame(c(rep(0.25, 4), rep(0.4, 4)), 
                    c(rep(0.4, 4), rep(0.25, 4)),
                    c(rep(c(rep(2, 2),rep(10, 2)),2)),
                    c(rep(c(10,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")
design = data.frame(.4,.25,2,100)
colnames(design) = c("alpha", "beta", "belif", "n")
## human
write.table(human_anno, file = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/human4_20_50/human_anno_20_50.txt", col.names = FALSE, row.names = F, quote = F)
for (i in 1:nrow(design)) {
  for (j in 1:30) {
    set.seed(j)
    d = unlist(design[i,])
    result = letsgen4(data = human, alpha = d[1], beta = d[2], belif = d[3], n = d[4])
    write.table(result$gene_rank, file =paste0("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/human4_20_50/human", paste(d[1]), "_", paste(d[2]), "_",
                                               paste(d[3]),"_", paste(d[4]), "_seed_",paste(j)  ,".txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
    write.table(result$activated_terms, file = paste0("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/human4_20_50/human", paste(d[1]), "_", paste(d[2]), 
                                                      "_",paste(d[3]),"_", paste(d[4]), "_seed_",paste(j), "key.txt"), 
                sep ="\t", row.names = F, quote = F, col.names = F)
  }
}
