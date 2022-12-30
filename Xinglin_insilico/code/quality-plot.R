## Quality measure PR plot for CRFE testing
## Xinglin Jia
## 8/3/2022

library(plyr)


## load ontology data
setwd("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico")

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


## combination design for alpha, beta, belief
design = data.frame(c(rep(0.1, 3), rep(0.4, 3)), 
                    c(rep(0.4, 3), rep(0.25, 3)),
                    c(rep(c(2,5,10), 2)))
colnames(design) = c("alpha", "beta", "belif")

preplot_v1 = function(data_dir, key_dir, anno, top_percent, design, e_or_h, n_points, m){
  ## pr-plot  
  setwd(data_dir)
  result = c()
  for (i in 1:nrow(design)) { #loop through each design
    d = unlist(design[i,])
    ready_to_ave = c()
    quality_all = c()
    proportion_all = c()
    for (j in 1:m) {  # loop through each seed, we did 30 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, ".txt") 
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, "key.txt")
      data = read.csv(filename, sep = "\t") # table from CRFE
      term = data[,9] # CRFE ranked terms
      key = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
      all_genes = read.csv(paste0(key_dir, filename), sep = "\t",header = F) # expression file, ranked genes
      
      ## put perturbed genes (gene expression > threshold) in a vector
      P = all_genes[all_genes[,2] >0,1]
      U = all_genes[all_genes[,2] <=0,1]
      P = unique(P)
      U = unique(U)
      
      n_X = floor(length(P)*top_percent/100) # number of top X percent of the perturbed genes
      X = all_genes[all_genes$V1 %in% P,][1:n_X,1]# top X percent perturbed genes
      
      ## quality and proportion
      quality = c()
      proportion = c()
      proportion_2 = c()
      C = c()
      for (q in 1:length(term)) { # for each term in CRFE result calculate quality and proportion of perturbed genes explained up to this term
        C_temp = unlist(anno[term[q]])
        C = unique(c(C_temp, C))
        nu = sum(C %in% X)/n_X
        deno = sum(C %in% U)/length(U)
        quality[q] = nu/deno
        proportion[q] = nu
        proportion_2[q] = deno
      }
      ## imputation
      f = tail(proportion,1)
      t = 1
      proportion_imputated = c(proportion, seq(from = f, to = t, by = ((t - f)/(n_points - floor(tail(proportion,1)*1000)))))
      f = tail(quality,1)
      quality_imputated =c(quality, seq(from = f, to = t, by = ((t - f)/(n_points - floor(tail(proportion,1)*1000)))))
      
      
      ## use nearest proportion to make 1000 points
      ## the 1000 points can be seen as 0.1%, 0.2%, ... 100% of the top 50% perturbed genes
      ## start form first point and first term's proportion  
      ## if first point (0.1%) <= first term's proportion, we know this term covered more % of the perturbed genes
      ## store quality score for first point, 1/1000 finished, move to next point
      ## until there is a point e.g. 5.3% > first term's proportion, we know term 1 no longer covers it, we move to next term and check 
      ## so that we can assign 1000 points for each stimulate data
      temp = c()
      ct = 1
      for (pt in 1:n_points) {
        if (pt/n_points <= proportion_imputated[ct]) {
          temp[pt] = quality_imputated[ct]
        } else {
          ct = ct +1
          temp[pt] = quality_imputated[ct]
        }
      }
      ready_to_ave = rbind(ready_to_ave, temp)
      #quality_all = rbind(quality_all, quality)
      #proportion_all = rbind(proportion_all, quality)
      #print(nrow(ready_to_ave))
    }
    ready_to_ave[sapply(ready_to_ave, is.infinite)] <- NA ## replaced Inf to NA and excluded for mean
    ave_qua = apply(ready_to_ave, 2, median)
    #ave_qua = colMeans(ready_to_ave, na.rm = T)
    result = rbind(result, as.vector(ave_qua))
    #print(nrow(result))
  }
  #result_all = list("result" = result, "ready_to_ave" = ready_to_ave, "quality_all" = quality_all, "proportion_all" = proportion_all)
  return(result)
} # now median # now P = pass threshold 


preplot_real_overexpress = function(data_dir, crfe_dir, anno, top_percent, n_points, threshold){
  crfe = read.csv(crfe_dir, sep = "\t") # table from CRFE
  term = crfe[,8] # CRFE ranked terms
  all_genes = read.csv(data_dir, sep = "\t",header = F) # expression file, ranked genes
  all_genes = all_genes[order(all_genes[,2], decreasing = TRUE),]
  U = all_genes[all_genes[,2] <=threshold,1]
  U = unique(U)
  P = all_genes[all_genes[,2] > threshold, 1]
  P = unique(P)
  
  n_X = floor(length(P)*top_percent/100) # number of top X percent of the perturbed genes
  X = all_genes[all_genes$V1 %in% P,][1:n_X,1]# top X percent perturbed genes
  
  ## quality and proportion
  quality = c()
  proportion = c()
  proportion_2 = c()
  C = c()
  for (q in 1:length(term)) { # for each term in CRFE result calculate quality and proportion of perturbed genes explained up to this term
    C_temp = unlist(anno[term[q]])
    C = unique(c(C_temp, C))
    nu = sum(C %in% X)/n_X
    deno = sum(C %in% U)/length(U)
    quality[q] = nu/deno
    proportion[q] = nu
    proportion_2[q] = deno
  }
  
  # temp = c()
  # ct = 1
  # for (ave in 1:n_points) {
  #   if (ave/n_points*100 <= proportion[ct]) {
  #     temp[ave] = quality[ct]
  #   } else {
  #     temp[ave] = quality[ct]
  #     ct = ct +1
  #   }
  # }
  # 
  result = list("quality" = quality, "percent_explained" = proportion)
  return(result)
} # now P = pass threshold


preplot_pr = function(data_dir, key_dir, design, e_or_h, n_points, m){
  ## pr-plot  
  setwd(data_dir)
  result = c()
  for (i in 1:nrow(design)) { #loop through each design
    d = unlist(design[i,])
    ready_to_ave = c()
    quality_all = c()
    proportion_all = c()
    for (j in 1:m) {  # loop through each seed, we did 30 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_seed_", j, ".txt") 
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]), "_seed_", j, "key.txt")
      data = read.csv(filename, sep = "\t") # table from CRFE
      term = data[,8] # CRFE ranked terms
      key = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
      all_genes = read.csv(paste0(key_dir, filename), sep = "\t",header = F) # expression file, ranked genes
      
      ## put perturbed genes (gene expression > threshold) in a vector
      P = all_genes[all_genes[,2] >0,1]
      U = all_genes[all_genes[,2] <=0,1]
      P = unique(P)
      U = unique(U)
      
      n_X = floor(length(P)*top_percent/100) # number of top X percent of the perturbed genes
      X = all_genes[all_genes$V1 %in% P,][1:n_X,1]# top X percent perturbed genes
      
      ## quality and proportion
      quality = c()
      proportion = c()
      proportion_2 = c()
      C = c()
      for (q in 1:length(term)) { # for each term in CRFE result calculate quality and proportion of perturbed genes explained up to this term
        C_temp = unlist(anno[term[q]])
        C = unique(c(C_temp, C))
        nu = sum(C %in% X)/n_X
        deno = sum(C %in% U)/length(U)
        quality[q] = nu/deno
        proportion[q] = nu
        proportion_2[q] = deno
      }
      ## imputation
      f = tail(proportion,1)
      t = 1
      proportion_imputated = c(proportion, seq(from = f, to = t, by = ((t - f)/(n_points - floor(tail(proportion,1)*1000)))))
      f = tail(quality,1)
      quality_imputated =c(quality, seq(from = f, to = t, by = ((t - f)/(n_points - floor(tail(proportion,1)*1000)))))
      
      
      ## use nearest proportion to make 1000 points
      ## the 1000 points can be seen as 0.1%, 0.2%, ... 100% of the top 50% perturbed genes
      ## start form first point and first term's proportion  
      ## if first point (0.1%) <= first term's proportion, we know this term covered more % of the perturbed genes
      ## store quality score for first point, 1/1000 finished, move to next point
      ## until there is a point e.g. 5.3% > first term's proportion, we know term 1 no longer covers it, we move to next term and check 
      ## so that we can assign 1000 points for each stimulate data
      temp = c()
      ct = 1
      for (pt in 1:n_points) {
        if (pt/n_points <= proportion_imputated[ct]) {
          temp[pt] = quality_imputated[ct]
        } else {
          ct = ct +1
          temp[pt] = quality_imputated[ct]
        }
      }
      ready_to_ave = rbind(ready_to_ave, temp)
      quality_all = rbind(quality_all, quality)
      proportion_all = rbind(proportion_all, quality)
      #print(nrow(ready_to_ave))
    }
    ready_to_ave[sapply(ready_to_ave, is.infinite)] <- NA ## replaced Inf to NA and excluded for mean
    ave_qua = apply(ready_to_ave, 2, median)
    #ave_qua = colMeans(ready_to_ave, na.rm = T)
    result = rbind(result, as.vector(ave_qua))
    print(nrow(result))
  }
  #result_all = list("result" = result, "ready_to_ave" = ready_to_ave, "quality_all" = quality_all, "proportion_all" = proportion_all)
  return(result)
} # now median # now P = pass threshold 

qr = function(data_dir, key_dir, anno, design, e_or_h, m){
  ## pr-plot  
  setwd(data_dir)
  terms = sample(names(anno))
  
  result = c()
  d = unlist(design[1,])
  pre_all = c()
  rec_all = c()
  for (j in 1:30) {  # loop through each seed, we did 30 stimulation for each design
    filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[5]), "_seed_", j, ".txt") 
    key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[5]), "_seed_", j, "key.txt")
    data = read.csv(filename, sep = "\t") # table from CRFE
    term = data[,9] # CRFE ranked terms
    term = unique(c(term, terms))
    key = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
    pre = c()
    rec = c()
    for (i in 1:length(term)) {
      pre[i] = sum(term[1:i]%in%key)/i
      rec[i] = sum(term[1:i]%in%key)/length(key)
    }
    pre_all = rbind(pre_all, pre)
    rec_all = rbind(rec_all, rec)
  }
  pre_final = colSums(pre_all)/m
  rec_final = colSums(rec_all)/m
  return(list("pre" = pre_final, "rec" = rec_final))
}

result = preplot(data_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/human/cut", 
         key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human/",
         anno = human,
         top_percent = 50,
         design = design,
         e_or_h = "human",
         n_points = 100)


par(mfrow=c(1,2))
## plot
x = c(1:100)
## alpha = 0.1, beta = 0.4
plot(x, result[3,], pch = " ", main = "Use real belief: human alpha=0.1 beta=0.4", 
     xlab = "Proportion of Top 50% Perturbed Genes that are Explained", ylab = "Quality score")
color = c("red", "blue", "black","purple", "green")
for (i in 1:3) {
  lines(x, result[i,], pch = ".", col = color[i])
}
legend(70, 60, legend=c("CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:3]
)

## alpha = 0.4, beta = 0.25
plot(x, result[4,], pch = " ", main = "Use real belief: human alpha=0.4 beta=0.25",
     xlab = "Proportion of Top 50% Perturbed Genes that are Explained", ylab = "Quality score")
for (i in 4:6) {
  lines(x, result[i,], pch = ".", col = color[i-3])
}
legend(70, 6, legend=c("CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:3]
)

### test for belief
design = data.frame(c(0.1),c(0.4),c(2))
colnames(design) = c("alpha", "beta", "belif")

b1 = preplot(data_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/ecoli_test_belief/ecoli_4_25_10/b1/cut", 
                 key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli/",
                 anno = ecoli,
                 top_percent = 50,
                 design = design,
                 e_or_h = "ecoli",
                 n_points = 100)

b2 = preplot(data_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/ecoli_test_belief/ecoli_4_25_10/b2/cut", 
             key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli/",
             anno = ecoli,
             top_percent = 50,
             design = design,
             e_or_h = "ecoli",
             n_points = 100)

b5 = preplot(data_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/ecoli_test_belief/ecoli_4_25_10/b5/cut", 
             key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli/",
             anno = ecoli,
             top_percent = 50,
             design = design,
             e_or_h = "ecoli",
             n_points = 100)

b10 = preplot(data_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/ecoli_test_belief/ecoli_4_25_10/b10/cut", 
             key_dir = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/ecoli/",
             anno = ecoli,
             top_percent = 50,
             design = design,
             e_or_h = "ecoli",
             n_points = 100)

## truth alpha = 0.4, beta = 0.25, b =10
x = c(1:100)
color = c("red", "blue", "black","purple", "green")
plot(x, b1, pch = " ", main = "Use diff b on same data. Truth: ecoli alpha=0.4 beta=0.25 b=10",
     xlab = "Proportion of Top 50% Perturbed Genes that are Explained", ylab = "Quality score")

lines(x, b1, pch = ".", col = color[1])
lines(x, b2, pch = ".", col = color[2])
lines(x, b5, pch = ".", col = color[3])
lines(x, b10, pch = ".", col = color[4])

legend(70, 20, legend=c("CRFE b=1","CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:4]
)


