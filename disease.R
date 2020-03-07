library('plot.matrix');
library('animation');
library('spatstat');
library('sp');
library('ggplot2');
library('rootSolve');

S0 <- 0.9;
initialInfected = 9;

# https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2016.0659
# assumes "susceptible depletion" kicking in
getR <- function(N, stillSusceptible, R0) (stillSusceptible/N) * R0;
getEffectiveR <- function(df, R0) round(getR(nrow(df), nrow(subset(df, status == 0)) * S0, R0), 2);

# https://mathematicsinindustry.springeropen.com/track/pdf/10.1186/s13362-019-0058-7
# "final size relation" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/
# s∞ = S(0)e^[–R0(1–s∞)]
getFinalUninfectedSusceptible <- function(R0) {
  return(multiroot(
    f = function(R0, Sinf) return(Sinf - (S0 * exp(-R0 * (1 - Sinf)))),
    start = 0,
    positive = TRUE,
    R0 = R0
  )$root);
}

getEffective0Ratio <- function(R0, R) if (R0 > 0) R/R0 else 0;

getPermUninfected <- function(R0, R, N) {
  return((getEffective0Ratio(R0, R) * (getFinalUninfectedSusceptible(R0) + (1 - S0))));
}

getInfectionProbabilities <- function(R0t, Rt, R0, R, N) {
  return(c(
    getEffective0Ratio(R0t, Rt), # infected
    getPermUninfected(R0, R, N) * (1 - getEffective0Ratio(R0t, Rt)) # perm uninfected
  ));
}

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

getColors <- function(statuses) lapply(statuses, function(s) getColor(s));

getDotSize <- function(size) 600 / (size * 4);

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
  df[1:initialInfected,]$status <- 1;
  cfrString = paste(CFR*100, '%', sep='');
  finalR <- getEffectiveR(df, R0);
  saveGIF({
    for (i in 1:steps) {
      R <- getEffectiveR(df, R0);
      print(ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
        geom_point(color=getColors(df$status), size = getDotSize(size)) +
        labs(
          title=name,
          subtitle=paste('R0=', R0, ', R=', R, ', CFR=', cfrString, sep='')
        ) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(colour = 'gray10', family='mono', hjust = 0.5, size=35, vjust = 0),
          plot.subtitle = element_text(colour = 'gray10', family='mono', hjust = 0.5, size=20, vjust = 0),
          plot.margin=unit(c(1.1, 0.25, 1, 0.25), "cm")
        )
      );
      copy <- df;
      apply(df, 1, function(d) {
        if (d[['status']] == 1) {
          R0t <- getStochRound(R0);
          Rt <- getStochRound(getEffectiveR(df, R0t));
          prob = getInfectionProbabilities(R0t, Rt, R0, R, size^2);
          # copy[is.na(copy)] <- 0;
          maxLength = min(R0t, length(copy[copy$status == 0,][,1]));

          if (max(prob) > 0 && maxLength > 0) {
            copy[copy$status == 0,][1:maxLength,]$status <<- sample(
              1:2,
              maxLength,
              replace = TRUE,
              prob
            );
          }
          copy[copy$x == d[['x']] & copy$y == d[['y']],]$status <<- sample(3:4, 1, prob=c(1-CFR, CFR));
        }
      });
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
  toRender[5,],
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