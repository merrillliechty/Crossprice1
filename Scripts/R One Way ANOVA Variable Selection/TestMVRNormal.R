## Test Chol and MV Random Normal Genearation ##
nSample = 1000000
B = diag(K)
B[1,2] = 0.25
B[2,1] = B[1,2]
B[1,5] = -0.5
B[5,1] = B[1,5]
B[2,3] = 0.75
B[3,2] = B[2,3]
B[2,6] = 0.25
B[6,2] = B[2,6]
B[4,5] = -0.65
B[5,4] = B[4,5]
cholB = chol(B)
tCholB = t(cholB)
ss = matrix(0,nrow = K, ncol = K)
sst = matrix(0,nrow = K, ncol = K)
for(t in 1:nSample)
{
  tb = cholB %*% t(t(rnorm(K)))
  ss = ss + tb %*% t(tb)
  tb = tCholB %*% t(t(rnorm(K)))
  sst = sst + tb %*% t(tb)
}

ss = ss*(1/nSample)
sst = sst*(1/nSample)
absSS_m_B = abs(ss - B)
absSSt_m_B = abs(sst - B)
sum(absSS_m_B)
sum(absSSt_m_B)