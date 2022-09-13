g2 = runif(10000,0,2.5)
b2 = runif(10000,0,1.5)
g3 = runif(10000,0,2.5)
b3 = runif(10000,-1.5,0)

r = 0.01
R = 1+r
sig = 0.04
beta = 10
g4 = 1.01
g1 = 0
b1 = 0
b4 = 0

generator_y = function(theta1,theta2,theta3,theta4,n){
  Y = matrix(NA, nrow=n, ncol=100)
  
  for(i in 1:n){
    Y[i,1] = 0
    Y[i,2] = 0
    Y[i,3] = 0
    for(j in 4:100){
      U1 = (Y[i,j-1]-R*Y[i,j-2])*(g1*Y[i,j-3]+b1-R*Y[i,j-2])
      U2 = (Y[i,j-1]-R*Y[i,j-2])*(theta1[i]*Y[i,j-3]+theta2[i]-R*Y[i,j-2])
      U3 = (Y[i,j-1]-R*Y[i,j-2])*(theta3[i]*Y[i,j-3]+theta4[i]-R*Y[i,j-2])
      U4 = (Y[i,j-1]-R*Y[i,j-2])*(g4*Y[i,j-3]+b4-R*Y[i,j-2])
      
      l1 = beta*U1
      l2 =  beta*U2
      l3 = beta*U3
      l4 =  beta*U4
      
      N1 = 1+exp(l2-l1)+exp(l3-l1)+exp(l4-l1)
      N2 = exp(l1-l2)+1+exp(l3-l2)+exp(l4-l2) 
      N3 = exp(l1-l3)+exp(l2-l3)+1+exp(l4-l3)
      N4 = exp(l1-l4)+exp(l2-l4)+exp(l3-l4)+1
      
      n1 = 1/N1
      n2 = 1/N2
      n3 = 1/N3
      n4 = 1/N4
      
      Y[i,j] = 1/R*(n1*(g1*Y[i,j-1]+b1)+n2*(theta1[i]*Y[i,j-1]+theta2[i])+n3*(theta3[i]*Y[i,j-1]+theta4[i])+n4*(g4*Y[i,j-1]+b4)+rcauchy(1,scale=sig))
    }
  }
  return(Y)
}

Y2 = generator_y(g2,b2,g3,b3,10000)
theta2 = cbind(g2,b2,g3,b3)
y2 = generator_y(rep(0.6,500),rep(0.65,500),rep(0.7,500),rep(-0.55,500),500)


write.table(Y2, 'Y2.csv', sep = ',')
write.table(theta2, 'theta2.csv', sep = ',')
write.table(y2, 'ytrue2.csv', sep = ',')

y11 = generator_y(rep(-0.7,500),rep(-0.4,500),rep(0.5,500),rep(0.3,500),500)
y22 = generator_y(rep(0.6,500),rep(0.65,500),rep(0.7,500),rep(-0.55,500),500)

write.table(y11, 'ytrue11.csv', sep = ',')
write.table(y22, 'ytrue22.csv', sep = ',')