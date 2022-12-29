design = data.frame(c(rep(0.1, 32), rep(0.3, 32)), 
                    c(rep(0.3, 32), rep(0.1, 32)),
                    c(rep(c(rep(2, 16),rep(10, 16)),2)),
                    c(rep(c(rep(5,4),rep(10,4),rep(50,4),rep(100,4)),4)),
                    c(rep(c("b1","b2","b5","b10"),16)))
colnames(design) = c("alpha", "beta", "belif","n","b")

dd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/human4/"
kd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human4/"

data_dir = paste0(dd,"b1/cut")
key_dir = kd
anno = human
design = design
e_or_h = "human"
plot_dim = c(2,4)
nsets = 16

ct1 = 1
for (sets in 1:nsets) {
  sub_d = design[ct1:(ct1+3),]
  par(mfrow=c(plot_dim[1],plot_dim[2])) 
  pre_all = list()
  rec_all = list()
  for(j in 1:(plot_dim[1]*plot_dim[2])){
    for (ct in 1:4) {
  ## pr-plot  
  data_dir = c(paste0(dd,"b1/cut"),paste0(dd,"b2/cut"),paste0(dd,"b5/cut"),paste0(dd,"b10/cut"))
  setwd(data_dir[ct])
  result = c()
    d = unlist(sub_d[ct,])
  #for (j in 1:1) {  # loop through each seed, we did 30 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, ".txt") 
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, "key.txt")
      data = read.csv(filename, sep = "\t") # table from CRFE
      term = data[,8] # CRFE ranked terms
      g_exp = read.csv(paste0(key_dir, filename),header = F,sep = "\t") # gene expression file
      true_terms = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
      key = c()
      for (tts in 1:length(true_terms)) {
        g_temp = unlist(anno[[true_terms[tts]]])
        key = c(key, g_temp)
      }
      #key = unique(key)## true postive genes + false negative genes
      #key = unique(key[g_exp[match(key, g_exp[,1]),2]>0])## true postive genes
      key = unique(g_exp[g_exp[,2]>0,1]) ## all postive genes
      
      test = c()
      for (crfe_ts in 1:length(term)) {
        t_temp = unlist(anno[[term[crfe_ts]]])
        t_temp = t_temp[order(g_exp[match(t_temp,g_exp[,1]),2], decreasing = T)]
        test = c(test, t_temp)
      }
      test = unique(test)
      
      pre = c()
      rec = c()
      for (i in 1:length(test)) {
        pre[i] = sum(test[1:i]%in%key)/i
        rec[i] = sum(test[1:i]%in%key)/length(key)
      }
      pre_all[[ct]] = pre
      rec_all[[ct]] = rec
}

  color = c("red", "blue", "black","purple")
  plot(c(0,1), c(0,1), pch = " ", xlab = "Recall", ylab = "Precision", main = paste0("human α=",paste(d[1])," β=",paste(d[2]), " b=",paste(d[3])," n=", paste(d[4])))
    for (i in 1:4) {
      lines(rec_all[[i]], pre_all[[i]], pch = ".", col = color[i], ylim = c(0,1), xlim = c(0,1))
      }
  } 
  ct1= ct1+4
}

#####pr-plot term base
## testing 
design = data.frame(c(rep(0.25, 4), rep(0.4, 4)), 
                    c(rep(0.4, 4), rep(0.25, 4)),
                    c(rep(c(rep(2, 2),rep(10, 2)),2)),
                    c(rep(c(10,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")

dd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/human4/"
kd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human4/"

data_dir = paste0(dd,"b1/temp")
key_dir = kd
anno = human
design = design
e_or_h = "human"
plot_dim = c(3,5)
nsets = nrow(design)

for (sets in 1:nsets) {
  d = unlist(design[sets,])
  par(mfrow=c(plot_dim[1],plot_dim[2])) 
  pre_all = list()
  rec_all = list()
  for(j in 1:(plot_dim[1]*plot_dim[2])){
    for (ct in 1:4) {
      ## pr-plot  
      data_dir = c(paste0(dd,"b1/temp"),paste0(dd,"b2/temp"),paste0(dd,"b5/temp"),paste0(dd,"b10/temp"))
      setwd(data_dir[ct])
      result = c()
      #for (j in 1:1) {  # loop through each seed, we did 30 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, ".txt") 
      data = read.csv(filename, sep = "\t") # table from CRFE
      test = data[,8] # CRFE ranked terms
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, "key.txt")
      key = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
      
      pre = c()
      rec = c()
      for (i in 1:length(test)) {
        pre[i] = sum(test[1:i]%in%key)/i
        rec[i] = sum(test[1:i]%in%key)/length(key)
      }
      pre_all[[ct]] = pre
      rec_all[[ct]] = rec
    }
    
    color = c("red", "blue", "black","purple")
    plot(c(0,1), c(0,1), pch = " ", xlab = "Recall", ylab = "Precision", main = paste0("human α=",paste(d[1])," β=",paste(d[2]), " b=",paste(d[3])," n=", paste(d[4])))
    for (i in 1:4) {
      lines(rec_all[[i]], pre_all[[i]], pch = ".", col = color[i], ylim = c(0,1), xlim = c(0,1))
    }
  } 
}

