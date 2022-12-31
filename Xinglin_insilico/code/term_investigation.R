prebox_sum = function(data_dir, key_dir, anno, design, e_or_h, n_points,m){
  ## pr-plot  
  setwd(data_dir)
  result_size = c()
  result_np = c()
  result_ratio = c()
  for (i in 1:nrow(design)) { #loop through each design
    d = unlist(design[i,])
    for (j in 1:m) {  # loop through each seed, we did 1000 stimulation for each design
      filename = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, ".txt") 
      key_name = paste0(e_or_h, paste(d[1]), "_", paste(d[2]), "_",paste(d[3]),"_",paste(d[4]), "_seed_", j, "key.txt")
      data = read.csv(filename, sep = "\t") # table from CRFE
      size = data[1:50,2] # CRFE ranked terms
      n_p = data[1:50,3]
      ratio = n_p/size
      result_size = rbind(result_size, size)
      result_np = rbind(result_np, n_p)
      result_ratio =rbind(result_ratio, ratio)
    }
    #post[is.na(post)] = 0
    #mean_p = colMeans(post)
    #result[is.na(result)] = 0
    #mean_rank = colMeans(result)
  }
  #return(mean_p)
  #return(mean_rank)
  result = list("size" = result_size, "np" = result_np, "ratio" = result_ratio)
  return(result)
}
dd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/human4/"
kd = "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/data/human4/"
seeds = 30 ## number of replicates

b1 = prebox_sum(data_dir = paste0(dd,"b1_n10rep/temp"), 
            key_dir = kd,
            anno = human,
            design = design,
            e_or_h = "human",
            m = seeds)

b2 = prebox_sum(data_dir = paste0(dd,"b2_n10rep/temp"), 
            key_dir = kd,
            anno = human,
            design = design,
            e_or_h = "human",m = seeds)

b5 = prebox_sum(data_dir = paste0(dd,"b5_n10rep/temp"), 
            key_dir = kd,
            anno = human,
            design = design,
            e_or_h = "human",m = seeds)

b10 = prebox_sum(data_dir = paste0(dd,"b10_n10rep/temp"), 
             key_dir = kd,
             anno = human,
             design = design,
             e_or_h = "human",m = seeds)

rows = c(1:30)
des = design
b1_me = c()
b2_me = c()
b5_me = c()
b10_me = c()
for (i in 1:ncol(b1$ratio)) {
  b1_me[i] = median(b1$ratio[rows,i], na.rm = T)
  b2_me[i] = median(b2$ratio[rows,i], na.rm = T)
  b5_me[i] = median(b5$ratio[rows,i], na.rm = T)
  b10_me[i] = median(b10$ratio[rows,i], na.rm = T)
}

x = seq(1:50)
color = c("red", "blue", "black","purple", "green")
plot(x, seq(0.02,1,0.02), pch = " ", main = paste0("human α=",paste(des[1])," β=",paste(des[2]), " b=",paste(des[3])," n=", paste(des[4])),
     xlab = "CRFE/MGSA term rank", ylab = "Median ratio of perturbed genes in the term, rep=30")
color = c("red", "blue", "black","purple", "green")
lines(x, b1_me, pch = ".", col = color[1])
lines(x, b2_me, pch = ".", col = color[2])
lines(x, b5_me, pch = ".", col = color[3])
lines(x, b10_me, pch = ".", col = color[4])
legend(30, 0.5, legend=c("CRFE b=1","CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:4])


b1_s = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b1/temp/adeno_symbol.txt", sep = "\t")[1:50,2]
b1_r = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b1/temp/adeno_symbol.txt", sep = "\t")[1:50,3]
b1_ra = b1_r/b1_s
b2_s = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b2/temp/adeno_symbol.txt", sep = "\t")[1:50,2]
b2_r = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b2/temp/adeno_symbol.txt", sep = "\t")[1:50,3]
b2_ra = b2_r/b2_s
b5_s = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b5/temp/adeno_symbol.txt", sep = "\t")[1:50,2]
b5_r = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b5/temp/adeno_symbol.txt", sep = "\t")[1:50,3]
b5_ra = b5_r/b5_s
b10_s = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b10/temp/adeno_symbol.txt", sep = "\t")[1:50,2]
b10_r = read.csv("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/result/real/b10/temp/adeno_symbol.txt", sep = "\t")[1:50,3]
b10_ra = b10_r/b10_s



x = seq(1:50)
color = c("red", "blue", "black","purple", "green")
plot(x, seq(2,100,2), pch = " ", main = "adeno_raymol",
     xlab = "CRFE/MGSA term rank", ylab = "# of perturbed genes")
color = c("red", "blue", "black","purple", "green")
lines(x, b1_r, pch = ".", col = color[1])
lines(x, b2_r, pch = ".", col = color[2])
lines(x, b5_r, pch = ".", col = color[3])
lines(x, b10_r, pch = ".", col = color[4])
legend("topright",legend=c("CRFE b=1","CRFE b=2", "CRFE b=5", "CRFE b=10"), 
       fill = color[1:4],cex = 0.75)





       