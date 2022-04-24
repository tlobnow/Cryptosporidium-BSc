library(colorRamps)
library(RColorBrewer)
library(dplyr)
#display.brewer.all()

Blue2purple2red <- c("#0000FF", "#400CF9", "#8019F3", "#B319BF", "#D90C5F", "#FF0000")
r <- c(0, 64, 128, 179, 217, 255)
g <- c(0, 12, 25, 25, 12,  0)
b <- c(255, 249, 243, 191,  95,   0)
BlPrpRd <- rgb(r, g, b, maxColorValue = 255)
 function (n, name = c("BlPrpRd")) 
 {
   BlPrpRd = rgb(c(0, 64, 128, 179, 217, 255),
                 c(0, 12, 25, 25, 12,  0),
                 c(255, 249, 243, 191,  95,   0),
                 maxColorValue = 255)
     name = match.arg(name)
     orig = eval(parse(text = name))
     rgb = t(col2rgb(orig))
     temp = matrix(NA, ncol = 3, nrow = n)
     x = seq(0, 1, , length(orig))
     xg = seq(0, 1, , n)
     for (k in 1:3) {
         hold = spline(x, rgb[, k], n = n)$y
         hold[hold < 0] = 0
         hold[hold > 255] = 255
         temp[, k] = round(hold)
     }
     palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
     palette
 }
 <environment: namespace:fBasics>

  beach <- function (n, name = c("beach.colors")) 
  {
    beach.colors = rgb(r,g,b,maxColorValue = 255)
    name = match.arg(name)
    orig = eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) {
      hold = spline(x, rgb[, k], n = n)$y
      hold[hold < 0] = 0
      hold[hold > 255] = 255
      temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    palette
  }
rm(BlPrpRd)
rm(Blue2purple2red)
#pal2 <- beach(n=50)
#par(mar = rep(0, 4))
#pie(rep(1, length(pal2)), col = pal2)
