
library(plotrix)

# the evaluation results of setting1 is stored in data1; the evaluation results of setting2 is stored in data2; ...
load("~/Documents/papers/groupwise mixture/JASA_casestudy/simulation/to_plot.RData")

layout(matrix(c(1,2,1,2,3,4,3,4,5,6,5,6, 7,7), 7, 2, byrow = TRUE))
par(mai=c(0.3, 0.5, 0.3, 0.3))

gap.plot(c(0:6), data1[,1], gap=c(0.2, 0.9), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1,0.2, 0.9, 1), ylim=c(0,1), xlab="iteration", ylab="", main="Setting 1", cex=0.3) 
gap.plot(c(0:6), data1[,2], gap=c(0.2, 0.9), pch=8,  type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:6), data1[,3], gap=c(0.2, 0.9), pch=12,  type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:6), data1[,4], gap=c(0.2, 0.9), pch=23,  type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:6), data1[,5], gap=c(0.2, 0.9), pch=24,  type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

gap.plot(c(0:8), data2[,1], gap=c(0.2, 0.9), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1, 0.2, 0.9, 1), ylim=c(0,1), xlab="iteration", ylab="", main="Setting 2", cex=0.3)
gap.plot(c(0:8), data2[,2], gap=c(0.2, 0.9), pch=8, type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:8), data2[,3], gap=c(0.2, 0.9), pch=12, type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:8), data2[,4], gap=c(0.2, 0.9), pch=23, type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:8), data2[,5], gap=c(0.2, 0.9), pch=24, type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

gap.plot(c(0:8), data3[,1], gap=c(0.2, 0.9), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1, 0.2, 0.9, 1), ylim=c(0,1), xlab="iteration", ylab="", main="Setting 3", cex=0.3)
gap.plot(c(0:8), data3[,2], gap=c(0.2, 0.9), pch=8, type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:8), data3[,3], gap=c(0.2, 0.9), pch=12, type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:8), data3[,4], gap=c(0.2, 0.9), pch=23, type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:8), data3[,5], gap=c(0.2, 0.9), pch=24, type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

gap.plot(c(0:7), data4[,1], gap=c(0.2, 0.9), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1, 0.2, 0.9, 1), ylim=c(0,1), xlab="iteration", ylab="", main="Setting 4", cex=0.3)
gap.plot(c(0:7), data4[,2], gap=c(0.2, 0.9), pch=8, type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:7), data4[,3], gap=c(0.2, 0.9), pch=12, type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:7), data4[,4], gap=c(0.2, 0.9), pch=23, type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:7), data4[,5], gap=c(0.2, 0.9), pch=24, type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

gap.plot(c(0:13), data5[,1], gap=c(0.25, 0.9), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 1), ylim=c(0,1), xlab="iteration", ylab="", main="Setting 5", cex=0.3)
gap.plot(c(0:13), data5[,2], gap=c(0.25, 0.9), pch=8, type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:13), data5[,3], gap=c(0.25, 0.9), pch=12, type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:13), data5[,4], gap=c(0.25, 0.9), pch=23, type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:13), data5[,5], gap=c(0.25, 0.9), pch=24, type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

gap.plot(c(0:13), data6[,1], gap=c(0.2, 0.85), pch=21, type="b", bgcol = "white", bg="#FF0000", col = "#FF0000", 
         ytics = c(0.1, 0.2, 0.9, 1), xtics = c(0, 5, 10, 15, 20),
         ylim=c(0,1), xlab="iteration", ylab="", main="Setting 6", cex=0.3)
gap.plot(c(0:13), data6[,2], gap=c(0.2, 0.85), pch=8, type="b", bgcol = "white", bg="#FFBF00", col = "#FFBF00", add=TRUE)
gap.plot(c(0:13), data6[,3], gap=c(0.2, 0.85), pch=12, type="b", bgcol = "white", bg="#33FF00", col = "#33FF00", add=TRUE)
gap.plot(c(0:13), data6[,4], gap=c(0.2, 0.85), pch=23, type="b", bgcol = "white", bg="#3300FF", col = "#3300FF", add=TRUE)
gap.plot(c(0:13), data6[,5], gap=c(0.2, 0.85), pch=24, type="b", bgcol = "white", bg="#CC00FF", col = "#CC00FF", add=TRUE)

par(mai=c(0,0,0,0))
plot.new()
legend(x="center", bty="n", c("MSE", "accuracy", "classification error rate", "cluster number error rate", "Gaussian center MSE"),
       pch = c(21,8,12,23,24), pt.bg=c("#FF0000", "#FFBF00", "#33FF00", "#3300FF", "#CC00FF"),
      col = c("#FF0000", "#FFBF00", "#33FF00", "#3300FF", "#CC00FF"), text.font=2, ncol = 2) #horiz = TRUE)

dev.off()



