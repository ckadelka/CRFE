belif = 5
fp_g = c(rep(1,100000))
inv1 = function(x) (-belif+sqrt(x*(belif+1)*(-belif+1)+belif^2))/(-belif+1)
fp_gs = c()
for (i in 1:length(fp_g)) {
  fp_gs[i] = inv1(runif(1,0,1))
}
hist(fp_gs, freq = F)






belif = 5
k = 5
fp_g = c(rep(1,100000))
b0 = 2*belif/(1+belif)
b1 = 2-2*b0

b1=2-4*belif/(1+belif)

ratio = b1*2/k
b1 = b1 -ratio 

inv1 = function(x) {(b1-2+sqrt((-b1+2)^2+4*b1^2*x))/b1^2}
inv1 = function(x) (belif-sqrt(x*(belif+1)*(belif-1)+belif^2))/(belif-1)

fp_gs = c()
for (i in 1:length(fp_g)) {
  fp_gs[i] = inv1(runif(1,0,1))
}
hist(fp_gs, freq = F)



b = 5
b1 = 4*b/(1+b)-2
b0 = 1-b1/2
ratio = b1*2/k

par(mfrow=c(2,5)) 

y = c()
k = 10
for (i in 1:k) {
  for (j in 1:100) {
    y[j] = j/100*b1+b0
  }
  b1 = b1 -ratio 
  b0 = 1-b1/2
  plot(c(seq(0.01,1,0.01)), y, ylim = c(0,2), type = "l")
}


