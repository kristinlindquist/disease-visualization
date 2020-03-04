library('plot.matrix');
library('animation');
library('spatstat');
library('sp');
library('ggplot2');
library('R0');

# https://mathematicsinindustry.springeropen.com/track/pdf/10.1186/s13362-019-0058-7
getFinalUninfectedSusceptible <- function(R0) multiroot(
  f = function(R0, S) return(S - exp(-R0 * (1 - S))),
  start = 0,
  positive = TRUE,
  R0 = R0
)$root;
getR <- function(R0) (R0 * (1 - getFinalUninfectedSusceptible(R0))); 

getColor = function(status) {
  # 0 - unexposed
  # 1 - infected
  # 2 - non-infected
  # 3 - undead
  # 4 - dead
  return(
    if (status == 0) 'ghostwhite'
    else if (status == 1) 'firebrick3'
    else if (status == 2) 'ghostwhite'
    else if (status == 3) 'azure3'
    else if (status == 4) 'gray10'
    else 'purple'
  );
}

getColors = function(statuses) lapply(statuses, function(s) getColor(s));

# http://coleoguy.blogspot.com/2016/04/stochasticprobabilistic-rounding.html
getStochRound <- function(x) {
  q = abs(x - trunc(x));
  adj = sample(0:1, size = 1, prob = c(1 - q, q));
  if(x < 0) adj = adj * -1;
  trunc(x) + adj;
}

pandemic <- function(size, steps, filename, name, R0, CFR) {
  df <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size, status = 0))));
  df[df$x == ceiling(size/2) & df$y == ceiling(size/2),]$status = 1;
  coordMap <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size))));
  coordMap$distance <- spDistsN1(as.matrix(coordMap[c("x", "y")]), c(ceiling(size/2), ceiling(size/2)));
  coordMap <- coordMap[order(coordMap$distance),];
  saveGIF({
    for (i in 1:steps) {
      R <- getR(R0);
      cfrString = paste(CFR*100, '%', sep='');
      print(ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
        geom_point(color=getColors(df$status), size = 5) +
        labs(caption=paste(name, ' (R0=', R0, ', CFR=', cfrString, ')', sep='')) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.caption = element_text(colour = 'gray10', family='mono', hjust = 0.5, size=30, vjust = 0),
          plot.margin=unit(c(1,1,1,1.2), "cm")
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
                  if (exposed <= getStochRound(R)) 1
                  else if (exposed <= getStochRound(R0)) 2
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
  }, movie.name=filename, interval = 2, ani.width = 600, ani.height = 600);
}

# R0s from https://en.wikipedia.org/wiki/Basic_reproduction_number
# seasonal flu R0 https://www.ncbi.nlm.nih.gov/pubmed/19545404
# CFR Covid-19 https://www.nytimes.com/2020/03/04/world/coronavirus-news.html
# Ebola CFR https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciz678/5536742
# Flu CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# SARS & MERS CFR https://www.worldometers.info/coronavirus/coronavirus-death-rate/
# smallpox CFR https://en.wikipedia.org/wiki/Smallpox
# measles CFR https://www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf
# Mumps CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# Rubella CFR (infants & in utero) https://www.cdc.gov/rubella/about/in-the-us.html
toRender <- data.frame(
  name = c('Covid-19', 'Ebola', 'Influenza', 'MERS', 'SARS', 'Mumps', 'Rubella', 'Smallpox', 'Measles'),
  fileName = c('Covid-19.gif', 'Ebola.gif', 'Influenza.gif', 'MERS.gif', 'SARS.gif', 'Mumps.gif', 'Rubella.gif', 'Smallpox.gif', 'Measles.gif'),
  CFR = c(0.034, 0.828, 0.001, 0.344, 0.096, 0.01, 0.001, 0.30, 0.002),
  R0 = c(2.6, 2, 1.3, 0.55, 3.5, 5.5, 6, 6, 15)
);

apply(
  toRender[1,],
  1,
  function(d) pandemic(
    size=30,
    steps=15,
    file=(d[['fileName']]),
    name=(d[['name']]),
    R0=(as.numeric(d[['R0']])),
    CFR=(as.numeric(d[['CFR']]))
  )
)


