library('plot.matrix');
library('animation');
library('spatstat');
library('sp');
library('ggplot2');

# https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2016.0659
# ... assumes "susceptible depletion" kicking in
getRt <- function(N, stillSusceptible, R0) ((N - stillSusceptible)/N) * R0;
getEffectiveR <- function(df, R0) round(getRt(nrow(df), nrow(subset(df, status != 0)), R0), 2);

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
  df$distance <- spDistsN1(as.matrix(df[c("x", "y")]), c(ceiling(size/2), ceiling(size/2)));
  df <- df[order(df$distance),];
  cfrString = paste(CFR*100, '%', sep='');
  saveGIF({
    for (i in 1:steps) {
      print(ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
        geom_point(color=getColors(df$status), size = 5) +
        labs(
          subtitle=paste(name),
          caption=paste(
            'R0=',
            R0,
            '   R=',
            getEffectiveR(df, R0),
            '   CFR=',
            cfrString,
            sep=''
          )
        ) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.caption = element_text(colour = 'gray10', family='mono', hjust = 0.5, size=20, vjust = 0),
          plot.subtitle = element_text(colour = 'gray20', family='mono', hjust = 0.5, size=30),
          plot.margin=unit(c(1.2,1,1,1.2), "cm")
        )
      );
      copy <- df;
      for (x1 in 1:size) {
        for (y1 in 1:size) {
          if (df[df$x == x1 & df$y == y1,]$status == 1) {
            thisR0 <- getStochRound(R0);
            thisR <- getStochRound(getEffectiveR(df, thisR0));
            if (thisR0 > 0) {
              effectiveNaughtRatio = thisR/thisR0;
              copy[copy$status == 0,][1:thisR0,]$status = sample(1:2, thisR0, replace = TRUE, prob=c(effectiveNaughtRatio, (1 - effectiveNaughtRatio)));
              copy[copy$x == x1 & copy$y == y1,]$status = sample(3:4, 1, prob=c(1-CFR, CFR));
            }
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
  name = c('MERS', 'Influenza', 'Covid-19', 'Ebola', 'SARS', 'Mumps', 'Rubella', 'Smallpox', 'Measles'),
  fileName = c('MERS.gif', 'Influenza.gif', 'Covid-19.gif', 'Ebola.gif', 'SARS.gif', 'Mumps.gif', 'Rubella.gif', 'Smallpox.gif', 'Measles.gif'),
  CFR = c(0.344, 0.001, 0.034, 0.828, 0.096, 0.01, 0.001, 0.30, 0.002),
  R0 = c(0.8, 1.3, 2.6, 2, 3.5, 5.5, 6, 6, 15)
);

apply(
  toRender,
  1,
  function(d) pandemic(
    size=31,
    steps=15,
    file=(d[['fileName']]),
    name=(d[['name']]),
    R0=(as.numeric(d[['R0']])),
    CFR=(as.numeric(d[['CFR']]))
  )
)