##########################################
# Name: Ayesha Kumbhare                  #
# ID: 604936823                          #
# Homework 9 .R file                     #
##########################################

X = read.csv("cleaned_data.csv")

# cleaned_data.csv is a file of 400 observations and 45 variables from the
# Family Income and Expenditure in the Phillipines dataset

X.cs = scale(X,scale=T) # centered and scaled 

##########################################
####### K-MEANS ANALYSIS #################

D = dist(X.cs) # distance matrix
class(D) # distance matrix has class dist

max(D) # the largest distance between rows is 28.0622

# find the two observations with the max distance
D1 = as.matrix(D) 
rows.most.distant = matrix(c(0,0), ncol=1)
rows.most.distant
for (i in 1:400){
  for (j in 1:400){
    if (D1[i,j] == max(D)){
      rows.most.distant = matrix(c(i,j), ncol=1)
    }
  }
}
rows.most.distant

# initial means for variables in cluster 1
c1 = as.numeric(as.vector(X.cs[210,]))

# initial means for variables in cluster 2
c2 = as.numeric(as.vector(X.cs[31,]))

# add c3 and c4

X.cs[210, ] # view initial vector for cluster 1
X.cs[31, ] # view initial vector for cluster 2

pastIndicator=400:1   # initial value for z 
indicator=1:400  

### note: we initialize this way to get the algorithm started

## While the two indicator vectors are different, keep going. 
while(sum(pastIndicator!=indicator)!=0) 
{
  pastIndicator=indicator; 
  
  #distance to current cluster centers
  dc1 =colSums((t(X.cs)-c1)^2)
  dc2=colSums((t(X.cs)-c2)^2)
  # add dc3 and dc4
  dMat=matrix(c(dc1,dc2),,2)
  
  #decide which cluster each point belongs to 
  indicator = max.col(-dMat)
  
  # update the cluster centers
  c1=colMeans(X.cs[indicator==1,])
  c2=colMeans(X.cs[indicator==2,])
  # update colmeans for c3 and c4
  
  # Make plot  	
  
}

## See the means of each cluster. 
c1   
c2   

indicator
pastIndicator

Z = pastIndicator 
X = cbind(X, Z) # add Z vector to original data

table(Z) # table indicating the proportion of each type in each cluster

# Mean vectors for each cluster
cluster1 = X[Z==1,]
cluster2 = X[Z==2,]
meanvector_c1 = colMeans(cluster1)
meanvector_c2 = colMeans(cluster2)
meanvector_c1
meanvector_c2

# Median vectors for each cluster
median_c1 = apply(cluster1, 2, median)
median_c2 = apply(cluster2, 2, median)
median_c1
median_c2

# Standard deviation vectors for each cluster
sigmahat_c1 = (346/347)*var(cluster1)
sigmahat_c2 = (52/53)*var(cluster2)
std_dev_c1 = sqrt(diag(sigmahat_c1))
std_dev_c2 = sqrt(diag(sigmahat_c2))
std_dev_c1
std_dev_c2

##########################################
############ PC ANALYSIS #################

## Sample variance-covariance matrix of scaled data is the 
## correlation matrix of the original unscaled and uncentered data 

S=var(X.cs)

eigen(S)
eigenvalues=eigen(S)$values
round(eigenvalues,3)  #for presentation purposes only 
sum(eigenvalues)

## compute what percentage of the sum of eigenvalues is each eigenvalue
round(100*eigenvalues/sum(eigenvalues),3)


cumsum(round(100*eigenvalues/sum(eigenvalues),3))


V=eigen(S)$vectors 
V.r = round(V,3)  # for presentation purposes 
V.r

## check that eigenvectors are orthonormal
## it should be 0 in the off diagonals and 1 in the 
## diagonals when you do cross-product 
#Off-diagonals represent pairwise inner products 
# and diagonals represent norm^2. 

round((t(V))%*%V,3)
# alternatively
round(crossprod(V,V))

## See relation of each data variable to the PC
## I fo with rounded values for presentation purposes only 

#rownames(V.r)=colnames(X)
V.r
round(crossprod(V.r))


##########################
# Obtain the principal components 
###########################
## Note: X.tilde is the PC matrix 
X.tilde=X.cs%*%V  # X.tilde are the PC
round(var(X.tilde),2)

head(X.tilde)

#X.tilde  # too long to print. Only to illustrate the first time


#######Check that the PC are orthogonal

round(t(X.tilde)%*%X.tilde,3)

#### Check the variance-covariance  of the PC 
round(var(X.tilde),3)

sx = scale(X.tilde,scale=T)

# Find the correlation between the components of your PC matrix 
X = X[,1:45] # X without Z vector
correlation = cor(X, X.tilde)

##########################################
####### NEWTON ALGORITHM #################

####### Create contour plot of the function
### to see if possible initial values 
### suggest themselves by looking at the 
## contour plots 

# log likelihool of lognormal
# L(x_1,x_2) = -sum(log(data)) - n*log(sqrt(x2[j])) - n/2*log(2*pi) 
# - 1/(2*x2[j])*sum((log(data) - x1[i])^2)
data=X$Total.Household.Income
n=400
x1=seq(10,15,by=0.01)
length(x1)
x2=seq(0,3,by=0.01)
f=matrix(0,nrow=length(x1),ncol=length(x2))
for(i in 1:length(x1)){ 
  for(j in 1:length(x2)) {
    f[i,j]=-sum(log(data)) - n*log(sqrt(x2[j])) - n/2*log(2*pi) - 1/(2*x2[j])*sum((log(data) - x1[i])^2)
  }}

contour(x1,x2,f,nlevels=300,xlab="mu",ylab="sigma^2")  


#### Now numerical optimization

### start defining gradient and hessian

xt=c(100,0)   # this helps us get started 
eps=0.00000001  # tolerance for xtp1-xt 
xtp1=c(12,1.5)    # vector with initial values for x1(mu) and x2 (sigma^2) 

xHist=matrix(xtp1,2,1) # save history of xt (Stores newton output)
fHist=c()
xHist
# objective function to maximize. The log likelihood of the normal 
# distribution. Notice that xtp1[1]= mean mu and xtp1[2] is the 
# variance. 
f=-sum(log(data)) - n*log(sqrt(xtp1[2])) - n/2*log(2*pi) - 1/(2*xtp1[2])*sum((log(data) - xtp1[1])^2)

# History of the objective function
#fHist=f

while(sum((xtp1-xt)^2)>eps){
  xt=xtp1
  xt
  #compute first and second derivatives
  gradient=as.vector(c((-1/xt[2]*sum(xt[1] - log(data))) , 
                       -n/sqrt((xt[2])) + sum((log(data)-xt[1])^2/ sqrt(xt[2])*xt[2]))) 
  gradient                                
  
  hessian=matrix(c(-n/xt[2],
                   2*sum(xt[1]-log(data))/sqrt(xt[2])*xt[2],
                   2*sum(xt[1]-log(data))/sqrt(xt[2])*xt[2], 
                   n/(xt[2])-((3*sum(log(data)-xt[1])^2)/xt[2]^2)),
                 ncol=2,nrow=2)
  hessian                                
  ###  
  # compute xtp1 solve(hessian)*gradient=hessian^{-1}*gradient)
  ###   
  xtp1=xt-solve(hessian)%*%gradient # 
  xtp1
  ###   
  #save history
  
  xHist=matrix(c(xHist,xtp1),2)
  xHist
  f=-sum(log(data)) - n*log(sqrt(xtp1[2])) - n/2*log(2*pi) - 1/(2*xtp1[2])*sum((log(data) - xtp1[1])^2)
  
  fHist=c(fHist,f)
}

xHist
fHist  
### Compute, after all this, the hessian 
### at the stationary point 
gradient                
hessian

### and find out whether the solution is a minimum, 
### a maximum or a saddle point. 

eigen(hessian) # the sign of the eigenvalues will tell you

### Now  find the standard errors of the estimators

## First the variances 
diag(solve(-hessian)) 

## Then the standard errors 

#sqrt(diag(solve(-hessian)))  

# confidence intervals for parameters

conf_int_sigma2 <- xt[2]+c(-1,1)*0.06568165*qnorm(.975,lower.tail=TRUE)
conf_int_sigma2 #95% conf interval for sigma^2
conf_int_mu <- xt[1]+c(-1,1)*0.06568165*qnorm(.975,lower.tail=TRUE)
conf_int_mu #95% conf interval for mu

