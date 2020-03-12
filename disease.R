library('animation');
library('spatstat');
library('sp');
library('ggplot2');
library('rootSolve');
library('gridExtra');
library('dplyr');
library('reshape2')

relativeHeight <- 0.85;
S0 <- 0.99;
initialInfected = 1;

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
getPermUninfected <- function(R0, R, N) getEffective0Ratio(R0, R) * (getFinalUninfectedSusceptible(R0) + (1 - S0));

getInfectionProbabilities <- function(R0t, Rt, R0, R, N) {
  c(
    getEffective0Ratio(R0t, Rt), # infected
    getPermUninfected(R0, R, N) * (1 - getEffective0Ratio(R0t, Rt)) # perm uninfected
  )
}

getOrder = function(status) {
  if (status == 1 | round(status, 0) == 3) 3
  else if (status == 0 | status == 2) 1
  else if (floor(status) == 3) 2
  else if (status == 4) 4
  else 1
}

getColor = function(status) {
  return(
    if (status == 0) 'ghostwhite' # 0 - unexposed
    else if (status == 1 ) 'firebrick1' # 1 - infected
    else if (round(status, 0) == 3) 'firebrick' # recovering
    else if (status == 2) 'ghostwhite' # 2 - non-infected
    else if (floor(status) == 3) 'azure3' # 3 - undead
    else if (status == 4) 'gray10' # 4 - dead
    else 'purple'
  );
}

getColors <- function(statuses) lapply(statuses, function(s) getColor(s));
getDotSize <- function(size, height) (height * relativeHeight) / (size * 4);

# http://coleoguy.blogspot.com/2016/04/stochasticprobabilistic-rounding.html
getStochRound <- function(x) {
  q = abs(x - trunc(x));
  adj = sample(0:1, size = 1, prob = c(1 - q, q));
  if(x < 0) adj = adj * -1;
  trunc(x) + adj;
}

getDots <- function(df, size, generation, name, R0, CFR, height) {
  ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
    geom_point(color=getColors(df$status), size = getDotSize(size, height)) +
    labs(
      title=name,
      subtitle=paste('Generation ', generation, ', R0=', R0, ', CFR=', CFR * 100, '%', sep='')
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(colour = 'gray10', family='Anonymous Pro', hjust = 0.5, size=35, vjust = 0),
      plot.subtitle = element_text(colour = 'gray10', family='Anonymous Pro', hjust = 0.5, size=25, vjust = 0),
      plot.margin=unit(c(1.1, 0.25, 0.5, 0.25), "cm")
    );
}

getHisto <- function(df) {
  df <- df[order(df$generation),][rev(order(apply(df, 1, function(e) getOrder(e[['status']])))),];
  ggplot(df, aes(x=generation, y=1, fill=getColors(df$status))) +
  geom_bar(stat = 'identity') +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none',
    panel.background = element_blank()
  )
}

epidemic <- function(size, generations, filename, name, R0, CFR, width = 650, interval = 2) {
  height <- width;
  df <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size, status = 0))));
  df[df$x == ceiling(size / 2) & df$y == ceiling(size / 2),]$status = 1;
  df$distance <- spDistsN1(as.matrix(df[c("x", "y")]), c(ceiling(size / 2), ceiling(size / 2)));
  df <- df[order(df$distance),];
  df[1:initialInfected,]$status <- 1;
  finalR <- getEffectiveR(df, R0);
  df_histogram <- cbind(df, generation=0);

  saveGIF({
    for (i in 1:generations) {
      R <- getEffectiveR(df, R0);
      df_histogram <- rbind(df_histogram, cbind(df, generation=i));
      print(grid.arrange(
        getDots(df, size, i, name, R0, CFR, height),
        getHisto(df_histogram),
        heights=c(height * relativeHeight, height * (1 - relativeHeight))
      ));
      
      copy <- df;
      apply(df, 1, function(d) {
        if (d[['status']] == 1) {
          R0t <- getStochRound(R0);
          Rt <- getStochRound(getEffectiveR(df, R0t));
          prob = getInfectionProbabilities(R0t, Rt, R0, R, size^2);
          maxLength = min(R0t, length(copy[copy$status == 0,][,1]));

          if (max(prob) > 0 && maxLength > 0) {
            copy[copy$status == 0,][1:maxLength,]$status <<- sample(
              1:2,
              maxLength,
              replace = TRUE,
              prob
            );
          }
          copy[copy$x == d[['x']] & copy$y == d[['y']],]$status <<- sample(3:4, 1, prob=c(1 - CFR, CFR));
        }
      });
      copy$status <- apply(copy, 1, function(e) ifelse(between(e[['status']], 3.0, 3.9), e[['status']] + 0.1, e[['status']]));
      df <- copy;
    }
  }, movie.name=filename, interval = interval, ani.width = width, ani.height = height);
}

epidemicFromβγ <- function(size, generations, filename, name, β, γ, CFR, width) {
  # β == transmission rate per infectious individual (e.g. 1.94 https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319)
  # γ == recovery rate; 1/γ == infectious period (e.g. 1/γ = 1.61 days https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319)
  R0 <- round(β / γ, 2);
  epidemic(size, generations, filename, name, R0, CFR, width);
}

epidemicFromRates <- function(size, generations, filename, name, τ, c, γ, CFR, width) {
  # τ = infection per contact
  # c = contact rate
  # β = contact rate * risk of infection
  β <- c * τ;
  R0 <- round(β / γ, 2);
  epidemic(size, generations, filename, name, R0, CFR, width);
}

# R0s from https://en.wikipedia.org/wiki/Basic_reproduction_number
# R0 SARS (1.77-1.85) https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1
# seasonal flu R0 https://www.ncbi.nlm.nih.gov/pubmed/19545404
# Covid-19 R0 1.4-2.5 https://www.who.int/news-room/detail/23-01-2020-statement-on-the-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov)
# Covid-19 R0 high (mean 3.28, median 2.79) https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319
# CFR Covid-19 highest (3.4%) https://www.who.int/dg/speeches/detail/who-director-general-s-opening-remarks-at-the-media-briefing-on-covid-19---3-march-2020
# CFR Covid-19 higher (2.3%) https://ourworldindata.org/coronavirus
# IFR Covid-19 (1.6%) https://www.medrxiv.org/content/10.1101/2020.03.04.20031104v1.full.pdf
# IFR Covid-19 low (0.3%–1%) https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200219-sitrep-30-covid-19.pdf?sfvrsn=3346b04f_2
# Ebola, SARS and MERS CFR https://ourworldindata.org/coronavirus
# Flu CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# smallpox CFR https://en.wikipedia.org/wiki/Smallpox
# measles CFR https://www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf
# Mumps CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# Rubella CFR (infants & in utero) https://www.cdc.gov/rubella/about/in-the-us.html
# renderFromR0 <- data.frame(
#   name = c('MERS', 'Influenza', 'Covid-19', 'Ebola', 'SARS', 'Mumps', 'Rubella', 'Smallpox', 'Measles'),
#   fileName = c('MERS.gif', 'Influenza.gif', 'Covid-19.gif', 'Ebola.gif', 'SARS.gif', 'Mumps.gif', 'Rubella.gif', 'Smallpox.gif', 'Measles.gif'),
#   CFR = c(0.34, 0.001, 0.023, 0.50, 0.10, 0.01, 0.001, 0.30, 0.002),
#   R0 = c(0.8, 1.3, 2.5, 2, 1.85, 5.5, 6, 6, 15)
# );

renderFromR0 <- data.frame(
   name = c('Covid-19 - Low', 'Covid-19 - High'),
   fileName = c('Covid-19-low.gif', 'Covid-19-high.gif'),
   CFR = c(0.01, 0.034),
   R0 = c(1.4, 3.28)
);

apply(
  renderFromR0,
  1,
  function(d) epidemic(
    size=51,
    generations=35,
    interval=1,
    width=650,
    file=(d[['fileName']]),
    name=(d[['name']]),
    R0=(as.numeric(d[['R0']])),
    CFR=(as.numeric(d[['CFR']]))
  )
)

renderFromβγ <- data.frame(
  name = c('Covid-19', 'Covid-19 Hypothetical - ¾ β', 'Covid-19 Hypothetical - ½ β', 'Covid-19 Hypothetical...'),
  fileName = c('Covid-19-1.gif', 'Covid-19-2.gif', 'Covid-19-3.gif', 'Covid-19-4.gif'),
  CFR = c(0.023, 0.023, 0.023, 0.023),
  β = c(1.6, 1.2, 0.8, (1.6 * 0.16)),
  γ = c(1/1.61, 1/1.61, 1/1.61, 1/1.61)
);

# apply(
#   renderFromβγ[4,],
#   1,
#   function(d) epidemicFromβγ(
#     size=101,
#     generations=30,
#     width=800,
#     file=(d[['fileName']]),
#     name=(d[['name']]),
#     CFR=(as.numeric(d[['CFR']])),
#     β=(as.numeric(d[['β']])),
#     γ=(as.numeric(d[['γ']]))
#   )
# )