library(plyr)


## load ontology data
setwd("/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/")

d = scan("GOecoli.txt", what = "", sep = "\n")
ecoli = strsplit(d, "[[:space:]]+")
names(ecoli) <- sapply(ecoli, `[[`, 1)
ecoli <- lapply(ecoli, `[`, -1)
ecoli = ecoli[sapply(ecoli, length) %in% c(20:200)]

### test for belief
prebox = function(data_dir, key_dir, anno, design, e_or_h, n_points,m){
  ## pr-plot  
  setwd(data_dir)
  result = c()
  post = c()
  for (i in 1:nrow(design)) { #loop through each design
    d = unlist(design[i,])
    ready_to_ave = c()
    for (j in 1:m) {  # loop through each seed, we did 1000 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, ".txt") 
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, "key.txt")
      data = read.csv(filename, sep = "\t") # table from CRFE
      term = data[,9] # CRFE ranked terms
      p = data[,4]
      key = unlist(read.csv(paste0(key_dir, key_name), header = F)) # true activated terms
      result = rbind(result, match(key, term))
      post = rbind(post, data[,4][match(key,term)])
    }
    #post[is.na(post)] = 0
    #mean_p = colMeans(post)
    #result[is.na(result)] = 0
    #mean_rank = colMeans(result)
  }
  #return(mean_p)
  #return(mean_rank)
  return(result)
}

design = data.frame(c(rep(0.25, 4), rep(0.4, 4)), 
                    c(rep(0.4, 4), rep(0.25, 4)),
                    c(rep(c(rep(2, 2),rep(10, 2)),2)),
                    c(rep(c(10,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")
design = design[c(1,3,5,7),]


design = data.frame(c(rep(0.25, 4), rep(0.4, 4)), 
                    c(rep(0.4, 4), rep(0.25, 4)),
                    c(rep(c(rep(2, 2),rep(10, 2)),2)),
                    c(rep(c(10,100),4)))
colnames(design) = c("alpha", "beta", "belif", "n")
design = design[c(2,4,6,8),]

dd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/test_results_small/CRFE/human4/"
kd = "/Users/JasonJia/Desktop/GOEnrichment_with_Xinglin_copy/Xinglin/data_small/human4/"
seeds = 30 ## number of replicates

b1 = prebox(data_dir = paste0(dd,"b1_n10rep/temp"), 
             key_dir = kd,
             anno = human,
             design = design,
             e_or_h = "human",
            m = seeds)

b2 = prebox(data_dir = paste0(dd,"b2_n10rep/temp"), 
             key_dir = kd,
             anno = human,
             design = design,
             e_or_h = "human",m = seeds)

b5 = prebox(data_dir = paste0(dd,"b5_n10rep/temp"), 
             key_dir = kd,
             anno = human,
             design = design,
             e_or_h = "human",m = seeds)

b10 = prebox(data_dir = paste0(dd,"b10_n10rep/temp"), 
              key_dir = kd,
              anno = human,
              design = design,
              e_or_h = "human",m = seeds)

## for n = 10
rows = c(1:30)
par(mfrow=c(3,4))
boxplot(b1[rows,1], b2[rows,1],b5[rows,1],b10[rows,1], main = ".4, .25, 2, 10, term 1", ylim = c(0,100))
boxplot(b1[rows,2], b2[rows,2],b5[rows,2],b10[rows,2], main = "term 2", ylim = c(0,100))
boxplot(b1[rows,3], b2[rows,3],b5[rows,3],b10[rows,3], main = "term 3", ylim = c(0,100))
boxplot(b1[rows,4], b2[rows,4],b5[rows,4],b10[rows,4], main = "term 4", ylim = c(0,100))
boxplot(b1[rows,5], b2[rows,5],b5[rows,5],b10[rows,5], main = "term 5", ylim = c(0,100))
boxplot(b1[rows,6], b2[rows,6],b5[rows,6],b10[rows,6], main = "term 6", ylim = c(0,100))
boxplot(b1[rows,7], b2[rows,7],b5[rows,7],b10[rows,7], main = "term 7", ylim = c(0,100))
boxplot(b1[rows,8], b2[rows,8],b5[rows,8],b10[rows,8], main = "term 8", ylim = c(0,100))
boxplot(b1[rows,9], b2[rows,9],b5[rows,9],b10[rows,9], main = "term 9", ylim = c(0,100))
boxplot(b1[rows,10], b2[rows,10],b5[rows,10],b10[rows,10], main = "term 10", ylim = c(0,100))

## for n = 100
par(mfrow=c(2,2))
for (para in 1:4) {
  sta = (para-1)*30+1
  en = sta+29
  rows = c(sta:en)
  des = design[para,]
  b1_me = c()
  b2_me = c()
  b5_me = c()
  b10_me = c()
  for (i in 1:ncol(b1)) {
    b1_me[i] = median(b1[rows,i], na.rm = T)
    b2_me[i] = median(b2[rows,i], na.rm = T)
    b5_me[i] = median(b5[rows,i], na.rm = T)
    b10_me[i] = median(b10[rows,i], na.rm = T)
  }
  
  x = seq(1:100)
  color = c("red", "blue", "black","purple", "green")
  plot(x, c(11:110), pch = " ", main = paste0("human α=",paste(des[1])," β=",paste(des[2]), " b=",paste(des[3])," n=", paste(des[4])),
       xlab = "True term rank", ylab = "Median rank in CRFE result, rep = 30")
  
  lines(x, b1_me, pch = ".", col = color[1])
  lines(x, b2_me, pch = ".", col = color[2])
  lines(x, b5_me, pch = ".", col = color[3])
  lines(x, b10_me, pch = ".", col = color[4])
}


    