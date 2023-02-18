library(matlib) ##-- plotEqn
###-- plotEqn3d for 3d plot
###-- Inverse 

library(MASS) ##for fractions in the matrix, for NULL function
library(pracma) ###This is for RREF
#### diag(8) ===  eye(8)

library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(rgl)
library(gMOIP)
library(igraph)
library(igraphdata)
library(lpSolve)
library(wCorr)


A<-matrix(c(1,-2,-1,2,1,-1,1,2,-2,1),2,5)
B<-matrix(c(3,1,2,4,-3,-4,0,2,-2,1),2,5)
C<-rbind(A,B)
C
nullspace(C)
null(C) ##same as nullspace(C)
Null(t(C)) ###nullspace(C)

### Matrix-Formation
M <- matrix(c(1,1,1,2,1,3,1,4),nrow=2,ncol=4)
M
##the entries are filled vertical wise.

M[4] ## the order in which we have filled the matrix, 4th-entry will be returned.

## Repeation
rep(1,10) ##repeating 1, 10 times.

diag(4) ## identity matrix of dim. 4x4


##To solve matrix equations.
library(matlib) ##install.packages("matlib)
##This is for many matrix functions.

A <- matrix(c(1,2,1,1), nrow = 2, ncol = 2)
b <- c(4, 5)
solve(A, b) ## Ax=b (solve() function will give us the matrix x.)
plotEqn(A,b)
## to plot equations obtained from AX=b.

## Example 42
A <- matrix(c(1,-2,-1,2,3,2,3,-2,1), nrow = 3, ncol = 3)
b <- c(6, -1, 2)
plotEqn3d(A,b, xlim=c(0,4), ylim=c(0,4))


## To get Reduced Echelon Form
## Example 44
A <- matrix(c(0, -1, 1, 0, 1, 1, 0, 1, 3, -4, 2, 0, -1, 0, 4, -4), 4, 4)
b <- c(1, 1, 5, -2)
showEqn(A, b)
echelon(A, b, verbose=TRUE, fractions=TRUE)####---Augmented matrix.





#########IMAGE-PROCESSING --- Uses Matrix Arithmetic

### Library
library(magick)

inp_img <- image_read("http://polytopes.net/Tora_color.png")
image_info(inp_img)
plot(inp_img)
inp_img

mod_img <- image_modulate(inp_img, brightness = 120, saturation = 100)
plot(mod_img)

##############

##Matrix Multiplication
A <- matrix(c(3, 0, -5, -1, -3, 4), nrow = 2, ncol = 3, byrow = TRUE)
B <- matrix(c(-5, 5, 2, 1, -2, 0), nrow = 3, ncol = 2, byrow = TRUE)
A %*% B ###---Dot Product

A <- matrix(c(4, -1, -5, 0, 1, -2), 2, 3, byrow = TRUE)
t(A)


#######################
library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(matlib)

## Standard deviation
sigma <- matrix(c(4,2,2,3), ncol = 2, nrow = 2) 

## Mean
mu <- c(1, 2)
n <- 1000
set.seed(123)
x <- rmvnorm(n = n, mean = mu, sigma = sigma) ##multivariate normal
d <- data.frame(x)
p2 <- ggplot(d, aes(x = X1, y = X2)) + geom_point(alpha = .5) + geom_density_2d()
p2

set.seed(1)
x <- rmvnorm(n = n, mean = mu, sigma = sigma) ##multivariate normal
mu <- c(1, 2)
y <- x - mu
dd <- data.frame(y)
p3 <- ggplot(dd, aes(x = X1, y = X2)) + geom_point(alpha = .5) + geom_density_2d()
p3

dd <- data.frame(y)
head(dd)
tran<-matrix(c(1/sqrt(2),-1/sqrt(2),1/sqrt(2),1/sqrt(2)),2,2)
dd<-as.matrix(dd)%*%tran
dd<-as.data.frame(dd)
head(dd)
p4<- ggplot(dd, aes(x = V1, y = V2)) + geom_point(alpha = .5) + geom_density_2d()
p4


dd <- data.frame(y)
head(dd)
tran<-matrix(c(0,-1,0,1),2,2)
dd<-as.matrix(dd)%*%tran
dd<-as.data.frame(dd)
head(dd)
p5<- ggplot(dd, aes(x = V1, y = V2)) + geom_point(alpha = .5) + geom_density_2d()
p5

dd <- data.frame(y)
head(dd)
tran<-matrix(c(0,1,0,-1),2,2)
dd<-as.matrix(dd)%*%tran
dd<-as.data.frame(dd)
head(dd)
p6<- ggplot(dd, aes(x = V1, y = V2)) + geom_point(alpha = .5) + geom_density_2d()
p6

dd <- data.frame(y)
head(dd)
tran<-matrix(c(1/sqrt(2),+1/sqrt(2),1/sqrt(2),-1/sqrt(2)),2,2)
dd<-as.matrix(dd)%*%tran
dd<-as.data.frame(dd)
head(dd)
p7<- ggplot(dd, aes(x = V1, y = V2)) + geom_point(alpha = .5) + geom_density_2d()
p7

ggarrange(p3,p4,p5,p6,p7,nrow=2,ncol=3) ## two plots in one area.
##################





############------Inverse
library(matlib)
A <- matrix(c(1,-2,-1,2,3,2,3,-2,1), nrow = 3, ncol = 3)
inv(A) ## solve(A)

library(MASS) ##for fractions in the inverse
fractions(inv(A))

t(A)
B <- matrix(c(0,0,0), nrow = 3, ncol = 1)
fractions(solve(A))####Here solve(A) will give us X such that AX=I.
## this will give inverse. of A
fractions(inv(A))
#########################


### Practical Applications
A <- matrix(c(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 
              0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 
              0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 
              0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1), nrow = 11, ncol = 18, byrow = TRUE)
b <- c(1298, 1948, 465 , 605 , 451 ,  338 , 260 ,  183 ,  282 ,  127 ,  535)
echelon(A, b, verbose=TRUE, fractions=TRUE) ####---reduced Echelon form of Augmented matrix.
rref(cbind(A,b))


G1 <- eye(8) ##Identity Matrix - 8X8 dim.
G2 <- matrix(rep(0, 80), 8, 10)
b2 <- c(266, 223, 140, 264, 137, 67, 130, 24)
G <- cbind(G1, G2, b2)
M <- rbind(E, G)
inv(M[,-19]) 



#####################----------New Way to find the inverse.
A <- matrix(c(1, 3, 2, 1, 2, 5, -2, -3, 2, 2, -3, -3, -4, -2, -1, 4), nrow = 4, ncol = 4, byrow=TRUE)
inv(A) ##method-1
M <- cbind(A, diag(4))
M
library(pracma)
R <- rref(M)
R
Ainv <- R[, 5:8]
Ainv ##method-2


###
E <- E[-11,]
G1 <- eye(8)
G2 <- matrix(1:80, 8, 10)
b2 <- c(266, 223, 140, 264, 137, 67, 130, 24)
G <- cbind(G1, G2, b2)
M <- rbind(E, G)
M2 <- M[,-19]
M3 <- cbind(M2, diag(18))
M4 <- rref(M3)
M4[,19:36]


####CONICS==>>
##a.x^2+b.x.y+c.y^2+d.x+e.y+f=0
library(conics) ##not able to download

##coefficients
v <- c(2, 2, 2, -20, -28, 10)
conicPlot(v)
##########


## Example 87
library(pracma)
A <- matrix(c(0, 1, 3, -1, -1, 1, -4, 0, 1, 0, 2, 4, 0, 1, 0, -4), 
            nrow = 4, ncol = 4, byrow = TRUE)
A
AI <- rowadd(A, 1, 2, -2) ##row2=row2+(-2)*row1
AI
det(AI)

##
dets<-numeric(length=500)
for(i in 1:500){
  k=i
  A <- matrix(c(k, k^2, k^3, k^2, k^4, k^6, k^3,k^6,k^9), 
              nrow = 3, ncol = 3, byrow = TRUE)
  dets[k]<-det(A)
}
mean(dets)

x<-1:500
plot(x,dets)

#################----Cramer's Rule
A <- matrix(c(0, 1, 3, -1, -1, 1, -4, 0, 1, 0, 2, 4, 0, 1, 0, -4), 
            nrow = 4, ncol = 4, byrow = TRUE)
b <- c(1, 1, 5, -2)
A1 <- A
A1[, 1] <- b
A2 <- A
A2[ ,2] <- b
A3 <- A
A3[ ,3] <- b
A4 <- A
A4[ ,4] <- b
x1 <- det(A1)/det(A)
x2 <- det(A2)/det(A)
x3 <- det(A3)/det(A)
x4 <- det(A4)/det(A)

solve(A,b)
####################


######vector spaces.
df<-USArrests
dim(df)
head(df)

df <- na.omit(df) ##clean the data by deleting missing values
df <- scale(df)
d <- dist(df, method = "euclidean")
d


########-------PLot vector space on the graph
## Example 98
library(rgl)
# Create some dummy data
dat <- replicate(2, 1:3)
# Initialize the scene, no data plotted
plot3d(dat, type = 'n', xlim = c(-1, 8), ylim = c(-1, 8), zlim = c(-10, 20), xlab = '', ylab = '', zlab = '')
# Define the linear plane
planes3d(2, 3, -1, 0, col = 'red', alpha = 0.6)
# Define the origin
points3d(x=0, y=0, z=0)


####-Example x+y+z=0
library(rgl)
# Create some dummy data
dat <- replicate(2, 1:3)
# Initialize the scene, no data plotted
plot3d(dat, type = 'n', xlim = c(-1, 8), ylim = c(-1, 8), zlim = c(-10, 20), xlab = '', ylab = '', zlab = '')
# Define the linear plane
planes3d(1,1,1, 0, col = 'red', alpha = 0.6)
# Define the origin
points3d(x=0, y=0, z=0)


######---NULL SPACE
######---Unique Solution if Null(A)=0
###-solution of exists if Ax=b b is the colspace(A)

library(pracma)
A <- matrix(c(1, -1, 4, 2, 0, -1, -1, -1, 5), nrow=3, ncol=3, byrow=TRUE)
b<-matrix(c(0,0,0),nrow = 3,ncol=1)
rref(A) ##row space
Solve(A,b)
plotEqn3d(A,b)


library(MASS) ####--Null matrix.
A <- matrix(c(1, -1, 4, 2, 0, -1, -1, -1, 5), nrow=3, ncol=3, byrow=TRUE)
Null(t(A)) #### method-1
nullspace(A) ## method-2

####---orth() gives the orthonormal basis of the column space.


t(rref(t(A))) ###the columns are the basis of the COLUMN SPACE
rref(A) ##the rows are the basis of the ROW SPACE
## an alternative method
orth(t(A)) # row space (as columns)
orth(A) # column space


## Non-zero Rows of rref form - are the basis of the matrix

A<-matrix(c(-5,5,-1,-7,-2,-4,1,3,4),nrow = 3,ncol=3,byrow=T)
A
Null(t(A))
rref(A) ## now zero rows- row space of A
rref(t(A)) ## now zero rows- row space of t(A) - column space of A
Rank(A)
ncol(Null(t(A))) ## dimension of null space
det(A)

####----CHange of BASIS

library(pracma)
B <- matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), 3, 3, byrow = TRUE)
B
v <- c(2, 4, 0)
v
inv(B) %*% v ## this is the change of basis. ([v]b=([B]b)^(-1) * v)  
solve(B,t(t(v))) ##this is also change of basis.

########heirarchy clustering
df <- USArrests
head(df)
df <- na.omit(df)
df <- scale(df)
d <- dist(df, method = "euclidean")
hc <- hclust(d)
sub_grp <- cutree(hc, k = 4)
library(factoextra)
plot(hc, cex = 0.6)


###-----K-Means Clustering
df <- USArrests
df <- na.omit(df)
df <- scale(df)
# nstart is the number of random starting points
clusters <- kmeans(df, 5, nstart = 10)
library(factoextra)
fviz_cluster(clusters, df)
### Practical Applications
library(animation)
kmeans.ani(df, 5)

fviz_nbclust(df, kmeans, method = "wss") ###no. of clusters


######3-------EIGEN
A <- matrix(c(2, 1, 1, 2), 2, 2)
eigen(A)

#############------Linear Programming
library(gMOIP)
obj <- c(2, 7)

A <- matrix(c(1, 1, 1, 7, -1, 1), nrow = 3, ncol = 2, byrow = TRUE)
b <- c(7, 35, 3)
plotPolytope(
  A,
  b,
  obj,
  type = rep("c", ncol(A)),
  crit = "max",
  faces = rep("c", ncol(A)),
  plotFaces = TRUE,
  plotFeasible = TRUE,
  plotOptimum = TRUE,
  labels = "coord"
)


### Practical Applications
A <- matrix(c(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1), nrow = 11, ncol = 18, byrow = TRUE)

f.con <- rbind(A, diag(18))
f.rhs <- c(1298, 1948, 465, 605, 451, 338, 260, 183, 282, 127, 535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

f.obj <- c(39.99, 126.27, 102.70, 81.68, 38.81, 71.99, 31.21, 22.28, 321.04, 145.36, 33.82, 154.05, 64.19, 87.90, 107.98, 65.45, 39.08, 167.38)
f.dir <- c("=", "=", "=", "=", "=", "=", "=", "=", "=", "=", "=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=")
library(lpSolve)
Sol <- lp ("min", f.obj, f.con, f.dir, f.rhs)
Sol$objval

##example----
A <- matrix(c(2,3,1,1,1,0,0,1,1,0), nrow = 5, ncol = 2, byrow = TRUE)
A

f.con <-A
f.rhs <- c(19,8,0,0,16)
f.rhs

f.obj <- c(5,7)
f.dir <- c( "<=", "<=", ">=", ">=", "<=")
library(lpSolve)
Sol <- lp ("max", f.obj, f.con, f.dir, f.rhs)
Sol$objval

################-------------Network
library(igraph)
library(igraphdata)

#####-------UNDIRECTED gRAPH
library(igraph)
G <- sample_gnp(10, 0.5) ###-Number of edges = nC2 * p
#here number is 23
plot(G)
as_adjacency_matrix(G) ###Symmetric


#####-------DIRECTED gRAPH
library(igraph)
N <- barabasi.game(10) ###-Number of edges = 
#here number is 9
plot(N)
as_adjacency_matrix(N) ###Not-Symmetric


##########------Forming A HISTOGRAM
Carbon<-read.table("http://stat4ds.rwth-aachen.de/data/Carbon.dat", header=TRUE)
head(Carbon)

breaks<-seq(2.0,10.0,by=1.0) ##gives a sequency from 2-10  -> 2,3,4,5,6,7,8,9,10
freq<-table(cut(Carbon$CO2,breaks,right=FALSE))
cbind(freq,freq/nrow(Carbon))
hist(Carbon$CO2,xlab="CO2",ylab="Proportion",freq=FALSE)# histogram
plot(density(Carbon$CO2))



###########------BoxPlot
Carbon<-read.table("http://stat4ds.rwth-aachen.de/data/Carbon.dat", header=TRUE)

summary(Carbon$CO2) # 1stQu=lowerquartile, 3rd Qu=upperquartile
c(mean(Carbon$CO2),sd(Carbon$CO2),quantile(Carbon$CO2,0.90))
boxplot(Carbon$CO2,xlab="CO2values",horizontal=TRUE)

Crime<-read.table("http://stat4ds.rwth-aachen.de/data/Murder2.dat", header=TRUE)
head(Crime)
boxplot(Crime$murder~Crime$nation,xlab="Murder Rate",horizontal=TRUE)
tapply(Crime$murder,Crime$nation,summary)


#########-----Bivariate Data
GS<-read.table("http://stat4ds.rwth-aachen.de/data/Guns_Suicide.dat", header=TRUE)
Guns<-GS$guns
Guns
Suicides<-GS$suicide
Suicides
plot(Guns,Suicides) # scatterplot with arguments x, y
abline(lm(Suicides~Guns))

cor(Guns,Suicides) # correlation

summary(lm(Suicides ~ Guns))


##########-----Contingency Table [MosaicPlot]
PID<-read.table("http://stat4ds.rwth-aachen.de/data/PartyID.dat", header=TRUE)
table(PID$race,PID$id) #forms contingency table (not shown here;see Table 1.1)
options(digits=2)
prop.table(table(PID$race,PID$id),margin=1) #Formargin=1, proportions
mosaicplot(table(PID$race,PID$id)) #graphical portrayalof cell sizes





######-------Probabiility_Distributions

###---Expectation
#--To illustrate, let’s randomly generate ten million values of Y = number of successes in n = 3 binary 
#“trials,” with probability 0.50 of a success on each trial, and then average the ten million observations.
#In R, the function rbinom(r, n, π) simulates r experiments, in each case finding the number of successes 
#in n binary trials with probability π of a success on each.
y<-rbinom(10000000,3,0.50) # simulate 10000000 times the no. of successes in 3 trials
mean(y)
sd(y)


###---BINOMIAL DISTRIBUTION
mu<-function(n,pi){n*pi} # function: binomial mean
sigma<-function(n,pi){sqrt(n*pi*1-pi)} # function: binomial standard deviation
Psd<-function(n,pi,k){  pbinom(mu(n,pi)+k*sigma(n,pi),n, pi) - pbinom(mu(n,pi)-k*sigma(n,pi),n,pi)  } 
## function:prob.within k std.dev.of mean

n=1500;pi=0.60;k=2
Psd(n, pi, 2) # probability within k=2 standard deviations of mean
Psd(n, pi, 3)

###---POISSON DISTRIBUTION
ppois(130, 100)-ppois(69, 100)
y<-seq(60, 140, 1)
plot(y, dpois(y, 100),type = "l")

###########------Binomial with large n and small p
y=seq(0, 1000, 1) # y values between 0 and 14 with increment of 1
dbinom(y,50000000,0.0000001)
plot(y,dbinom(y,50000000,0.0000001),type="l")
plot(y,dpois(y,5),type="l")


###----GAMMA
y<-seq(0, 1000, 1)
plot(y, dgamma(y,shape =1 ,scale=100),type = "l")



#############-------Joint and Conditional Distributions
GSS <- read.table("http://stat4ds.rwth-aachen.de/data/GSS2018.dat", header=T)
head(GSS)
gender <- factor(GSS$SEX, levels=c(1,2), labels = c("Male","Female"))
smallgap <- factor(GSS$SMALLGAP, levels=c(1:5), labels = c("strongly agree", "agree","neutral","disagree","strongly disagree"))
fairsociety <- table(gender,smallgap) # frequency table
fairsociety

joint.prob <- fairsociety/sum(fairsociety) # joint cell proportions (total = 1)
## add to joint prob. the row and column marginal probabilities:
round(addmargins(joint.prob),3)
barplot(joint.prob, density=30, main="In a fair society, differences
in people's standard of living should be small",xlab="", ylab="proportions",ylim=c(0,0.30),legend=rownames(joint.prob))

#######-------Correlation
probabilities <- c(0.2,0.1,0.0,0.1,0.2,0.1,0.0,0.1,0.2)
x <- c(1,1,1,2,2,2,3,3,3)
y <- c(1,2,3,1,2,3,1,2,3) #scores

library(wCorr)
weightedCorr(x, y, weights=probabilities, method="polyserial")

########################################----------------------------------------------------------

x<-rnorm(1000,70, 10) # simulate n=1000 from normal with mean=70 and standard dev=10
y<-rnorm(1000,70 + 0.6*(x-70),6) #E(Y|x)=70+0.6*(x-70)
plot(x,y)



runif(10,1,2)
sample()


########-------CLT on binomial
CLT_binom <- function(B,n,p) {
  # B: number of iterations used to approximate the distribution of Xmean
  # n: sample size
  # p: success probability pi
  Y <- rbinom(B,n,p)
  Ymean <- Y/n # vector (length B) with the p-estimates: Algorithm 1 (2)
  var.mean <-p*(1-p)/n # variance of the estimator of p
  p.MC <- mean(Ymean) # Monte Carlo estimate of p
  varp.MC <- var(Ymean) # MC variance estimate of var.mean
  h <- hist(Ymean, col = "gray", probability=TRUE, main=paste("n=",n))
  xfit<-seq(0, 1,length=5000)
  yfit<-dnorm(xfit,mean=p,sd=sqrt(p*(1-p)/n))
  gr <- lines(xfit, yfit, col="blue",lwd=2)
  list(var.mean=var.mean, p.MC=p.MC, varp.MC=varp.MC) 
}

par(mfrow=c(2,2)) # multiple graphs layout in a 2x2 table format
CLT_binom(100000, 10, 0.3)
CLT_binom(100000, 30, 0.3)
CLT_binom(100000, 100, 0.3)
CLT_binom(100000, 1000, 0.3)


#########--------CLT ON POISSOn
pois_CLT <- function(n, mu, B) {
  # n: vector of 2 sample sizes [e.g. n <- c(10, 100)]
  # mu: mean parameter of Poisson distribution
  # B: number of simulated random samples from the Poisson
  par(mfrow = c(2, 2))
  for (i in 1:2){
    Y <- numeric(length=n[i]*B)
    Y <- matrix(rpois(n[i]*B, mu), ncol=n[i])
    Ymean <- apply(Y, 1, mean) # or, can do this with rowMeans(Y)
    barplot(table(Y[1,]), main=paste("n=", n[i]), xlab="y",
            col="lightsteelblue") # sample data dist. for first sample
    hist(Ymean, main=paste("n=",n[i]), xlab=expression(bar(y)),
         col="lightsteelblue") # histogram of B sample mean values
    
  } 
}

# implement:with 100000 random sample sizes of 10 and 100, mean = 0.7
n <- c(10, 100)
pois_CLT(n, 0.7, 100000)




Y<-matrix(rpois(10*100000,0.7), ncol=10)
ncol(Y)
nrow(Y)
# simulate Poisson rv's with mean 0.7
Ymed<-apply(Y,1, median) # find median for each sample of size 10
hist(Ymed,freq=FALSE) # histogram of the 100,000 sample medians

