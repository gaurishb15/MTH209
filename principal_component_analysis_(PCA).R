#### Used to reduce dimensionality

library(ISLR)
data(Auto)
summary(Auto)
head(Auto)

auto<-Auto[,1:7]
summary(auto)

pc <- prcomp(auto, scale=TRUE)
plot(pc) ###Shows Importance of different variables

library(ggfortify)
autoplot(pc, data=Auto, colour = "origin")
