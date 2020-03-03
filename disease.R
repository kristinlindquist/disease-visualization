library('plot.matrix');
library('animation');
library('spatstat');
library('sp');
library('ggplot2');

R0 <- 3;
totalExposed <- 10;
CFR <- 0.1;

getColor = function(status) {
  # 0 - unexposed
  # 1 - infected
  # 2 - non-infected
  # 3 - undead
  # 4 - dead
  return(
    if (status == 0) 'ghostwhite'
    else if (status == 1) 'coral'
    else if (status == 2) 'darkseagreen2'
    else if (status == 3) 'azure4'
    else if (status == 4) 'firebrick3'
    else 'purple'
  );
}

getColors = function(statuses) lapply(statuses, function(s) getColor(s));

pandemic <- function(size, steps, filename) {
  df <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size, status = 0))));
  df[df$x == ceiling(size/2) & df$y == ceiling(size/2),]$status = 1;
  coordMap <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size))));
  coordMap$distance <- spDistsN1(as.matrix(coordMap[c("x", "y")]), c(ceiling(size/2), ceiling(size/2)));
  coordMap <- coordMap[order(coordMap$distance),];
  saveGIF({
    for (i in 1:steps) {
      print(ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
        geom_tile(color="black", size = 0.25) +
        labs(title="Pandemic") +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank()
        )
      );
      copy <- df;
      for (x1 in 1:size) {
        for (y1 in 1:size) {
          if (df[df$x == x1 & df$y == y1,]$status == 1) {
            exposed <- 0;
            for(j in 1:nrow(coordMap)) {
              v <- coordMap[j,];
              if (copy[copy$x == v$x & copy$y == v$y,]$status == 0) {
                copy[copy$x == v$x & copy$y == v$y,]$status =
                  if (exposed <= R0) 1
                  else if (exposed < totalExposed) 2
                  else 0;
                
                exposed = exposed + 1;
              }
            }
            coordMap <- coordMap[coordMap$x != x1 | coordMap$y != y1,];
            copy[copy$x == x1 & copy$y == y1,]$status = sample(3:4, 1, prob=c(1-CFR, CFR));
          }
        }
      }
      df <- copy;
    }
  }, movie.name=filename, interval = 0.5, ani.width = 600);
}

pandemic(size=15, steps=15, file="pandemic.gif")
