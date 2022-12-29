###CRFE/MGSA test for likelihood 
## Xinglin Jia
## 11/1/2022
alpha_i = c()
beta_i = c()
alpha = 0.3
beta=  0.3
belif = 200000
n = 10
deno=0
for (i in 1:n) {
  deno=0
  for (j in 1:n) {
    deno = deno + (n-j+ (j-1) * belif)
    #print(deno)
  }
  beta_i[i] = beta*(n-i + (i-1)*belif)*n/deno
  alpha_i[i] = alpha*(n-i + (i-1)*belif)*n/deno
}

## gen Big/Small and top/mid/bot terms. Bt, Bm, Bb, St, Sm, Sb
term_gen = function(alpha, beta, b, n_terms, Big, Small){
  alpha_i = c()
  beta_i = c()
  Big = 200
  Small = 20
  deno=0
  for (i in 1:n_terms) {
    deno=0
    for (j in 1:n_terms) {
      deno = deno + (n-j+ (j-1) * belif)
      #print(deno)
    }
    beta_i[i] = beta*(n_terms-i + (i-1)*belif)*n_terms/deno
    alpha_i[i] = alpha*(n_terms-i + (i-1)*belif)*n_terms/deno
  }
  Bt = Big * (1-beta_i[1])
  Bm = Big * (1-beta)
  Bb = Big * (1-beta_i[n_terms])
  St = Small * (1-beta_i[1])
  Sm = Small * (1-beta)
  Sb = Small * (1-beta_i[n_terms])
  out = c(Bt, Bm, Bb, St, Sm, Sb)
  betas = c(beta_i[1],beta,beta_i[n_terms])
  return(list("out" = out, "betas" = betas, "alpha_i" = alpha_i))
}

t = term_gen(0.1, 0.3, 5, 100, 200, 20)

  

l_Big = c()
l_Small = c()
for (i in 1:3) {
  l_Big[i] = log((1-t$betas[i])^t$out[i] * t$betas[i]^(200-t$out[i]))
  l_Small[i] = log((1-t$betas[i])^t$out[i+3] * t$betas[i]^(20-t$out[i+3]))
}


l =function(p,a){
  b = .3
  u = a-p
  lh = log(b^u*(1-b)^p)
  return(lh)
}
c1 =c()
c2 = c()
c3 = c()
c4 = c()
for (i in 1:3) {
  v = c(180,140,100)
  c1[i] = l(v[i], 200)
}
 
for (i in 1:3) {
  v = c(18,14,10)
  c2[i] =l(v[i], 20)
}
 
 l =function(p,a){
   b = .9
   u = a-p
   lh = log(b^u*(1-b)^p)
   return(lh)
 }
 
 for (i in 1:3) {
   v = c(180,140,100)
   c3[i] = l(v[i], 200)
 }
 
 for (i in 1:3) {
   v = c(18,14,10)
   c4[i] =l(v[i], 20)
 }
 
 
 