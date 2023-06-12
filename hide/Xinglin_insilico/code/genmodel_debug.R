n = 10
alpha = 0.3
beta = 0.3
belif = 5
data = human

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
  
  for (k in 1:length(no_dup_p)) {
    loc = sum(runif(1,0,1)>(density_in_line))
    loc_db[k] = loc
    no_dup_p_gs[k] = runif(1, tp_lowers[loc+1], tp_uppers[loc+1])
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
p_g = c(tp_g, fp_g)
p_gs = c(tp_gs, fp_gs)

genes = c(p_g, up_g)
score = c(p_gs, up_gs)

gene_rank = data.frame(genes, score)
gene_rank = gene_rank[order(gene_rank$score, decreasing = TRUE),]

###debugs
tp_nn = c()
sum_temp = 0
for (i in 1:length(result$tp_n)) {
  sum_temp = sum(sum_temp+result$tp_n[i])
  tp_nn[i] = sum_temp
}

par(mfrow=c(2,5)) 
hist(result$tp_gs[1:tp_nn[1]], freq = F)
for (i in 1:(length(result$tp_n)-1)) {
  hist(result$tp_gs[tp_nn[i]:tp_nn[i+1]], freq = F) 
}




