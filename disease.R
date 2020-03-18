library('animation');
library('spatstat');
library('sp');
library('ggplot2');
library('rootSolve');
library('gridExtra');
library('dplyr');
library('reshape2')

relativeHeight <- 0.85;
S0 <- 1.0;
initialInfected <- 1;
standardInterval <- 2;

# https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2016.0659
# assumes "susceptible depletion" kicking in
getR <- function(N, stillSusceptible, R0) (stillSusceptible / N) * R0;
getEffectiveR <- function(df, R0) round(getR(nrow(df), nrow(subset(df, status == 0)) * S0, R0), 2);

# https://mathematicsinindustry.springeropen.com/track/pdf/10.1186/s13362-019-0058-7
# "final size relation" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/
# s∞ = S(0)e^[–R0(1–s∞)]
getFinalUninfectedSusceptible <- function(R0) {
  multiroot(
    f = function(R0, Sinf) return(Sinf - (S0 * exp(-R0 * (1 - Sinf)))),
    start = 0,
    positive = TRUE,
    R0 = R0
  )$root
}

getEffective0Ratio <- function(R0, R) if (R0 > 0) R / R0 else 0;

getPermUninfected <- function(R0) getFinalUninfectedSusceptible(R0) + (1 - S0);

getInfectionProbabilities <- function(R0t, Rt, R0) {
  c(
    (1 - getPermUninfected(R0)) * (1 - getEffective0Ratio(R0t, Rt)), # remain susceptible
    getEffective0Ratio(R0t, Rt), # infected
    getPermUninfected(R0) * (1 - getEffective0Ratio(R0t, Rt)) # uninfected due to final size relation
  )
}

getNextStatus <- function(status, serialInterval, recoveryInterval, CFR) {
  if (between(status, 3.0, 3.9)) {
    status = status + 0.01;
    if (status > 3 + (getStochRound(recoveryInterval / standardInterval) * 0.01)) {
      return(sample(5:4, 1, prob=c(1 - CFR, CFR)));
    }
  }

  return(status);
}

getOrder = function(status) {
  if (status == 0 | status == 2) 1
  else if (status == 5) 2
  else if (status == 1) 3
  else if (floor(status) == 3) 4
  else 5
}

getColor = function(status) {
  return(
    if (status == 0) 'ghostwhite' # 0 - unexposed
    else if (status == 1) 'lightsalmon' # 1 - infected
    else if (status == 2) 'ghostwhite' # 2 - non-infected
    else if (floor(status) == 3) 'salmon3' # recovering
    else if (status == 4) 'gray10' # 4 - dead
    else if (status == 5) 'azure3' # 5 - recovered
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

getDots <- function(df, size, name, R0, CFR, days, serialInterval, height) {
  ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
    geom_point(color=getColors(df$status), size = getDotSize(size, height)) +
    labs(
      title=name,
      subtitle=paste('Day ', days, ', R0=', R0, ', CFR=', CFR * 100, '%', ', SI=', serialInterval, ' days', sep='')
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
  df <- df[order(df$iteration),][rev(order(apply(df, 1, function(e) getOrder(e[['status']])))),];
  ggplot(df, aes(x=iteration + 1, y=1, fill=getColors(df$status))) +
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

getDays <- function(iteration, interval) round(iteration * interval, 0);

epidemic <- function(size, iterations, filename, name, R0, CFR, serialInterval = 5, recoveryInterval = 14, width = 650, interval = 2) {
  height <- width;
  df <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size, status = 0))));
  df$distance <- spDistsN1(as.matrix(df[c("x", "y")]), c(ceiling(size / 2), ceiling(size / 2)));
  df <- df[order(df$distance),];
  df[1:initialInfected,]$status <- 1;
  percentGeneration <- standardInterval / serialInterval;

  saveGIF({
    for (i in 0:iterations) {
      dfHistogram <- if (exists("dfHistogram")) rbind(dfHistogram, cbind(df, iteration=i)) else cbind(df, iteration=i);
      print(grid.arrange(
        getDots(df, size, name, R0, CFR, getDays(i, standardInterval), serialInterval, height),
        getHisto(dfHistogram),
        heights=c(height * relativeHeight, height * (1 - relativeHeight))
      ));
      
      copy <- df;
      apply(df, 1, function(d) {
        if (d[['status']] == 1 && getStochRound(percentGeneration) == 1) {
          R0t <- getStochRound(R0); # rounded R0 per iteration
          Rt <- getStochRound(getEffectiveR(df, R0t)); # rounded effective R per iteration
          prob = getInfectionProbabilities(R0t, Rt, R0);
          maxLength = min(R0t, length(copy[copy$status == 0,][,1]));

          if (max(prob) > 0 && maxLength > 0) {
            copy[copy$status == 0,][1:maxLength,]$status <<- sample(
              0:2,
              maxLength,
              replace = TRUE,
              prob
            );
          }
          copy[copy$x == d[['x']] & copy$y == d[['y']],]$status <<- 3;
        }
      });
      copy$status <- apply(copy, 1, function(e) getNextStatus(e[['status']], serialInterval, recoveryInterval, CFR));
      df <- copy;
    }
  }, movie.name=filename, interval = interval, ani.width = width, ani.height = height);
}

# R0s from https://en.wikipedia.org/wiki/Basic_reproduction_number
# R0 SARS (1.77-1.85) https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1
# seasonal flu R0 https://www.ncbi.nlm.nih.gov/pubmed/19545404
# Covid-19 R0 1.4-2.5 https://www.who.int/news-room/detail/23-01-2020-statement-on-the-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov)
# Covid-19 R0 high (mean 3.28, median 2.79) https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319
# CFR Covid-19 highest (3.4%) https://www.who.int/dg/speeches/detail/who-director-general-s-opening-remarks-at-the-media-briefing-on-covid-19---3-march-2020
# CFR Covid-19 higher (2.3%) https://ourworldindata.org/coronavirus
# IFR Covid-19 (1.6%) https://www.medrxiv.org/content/10.1101/2020.03.04.20031104v1.full.pdf
# IFR Covid-19 (0.94%) https://institutefordiseasemodeling.github.io/nCoV-public/analyses/first_adjusted_mortality_estimates_and_risk_assessment/2019-nCoV-preliminary_age_and_time_adjusted_mortality_rates_and_pandemic_risk_assessment.html
# IFR Covid-19 low (0.3%–1%) https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200219-sitrep-30-covid-19.pdf?sfvrsn=3346b04f_2
# Ebola, SARS and MERS CFR https://ourworldindata.org/coronavirus
# Flu CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# smallpox CFR https://en.wikipedia.org/wiki/Smallpox
# measles CFR https://www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf
# Mumps CFR https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates
# Rubella CFR (infants & in utero) https://www.cdc.gov/rubella/about/in-the-us.html
# MERS Serial (13 days) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930778/
# flu serial (3 days) https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200306-sitrep-46-covid-19.pdf?sfvrsn=96b04adf_2
# Covid-19 serial (5-6 days) https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200306-sitrep-46-covid-19.pdf?sfvrsn=96b04adf_2
# Ebola serial (15.3 days) https://www.sciencedirect.com/science/article/pii/S1755436515000341
# SARS serial interval (8.4 days) https://dash.harvard.edu/bitstream/handle/1/25620506/Transmission%20dynamics%20and%20control%20of%20severe%20acute%20respiratory%20syndrome.pdf?sequence=1
# Mumps serial interval (20 days) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223546/
# Rubella serial interval (18 days) https://academic.oup.com/aje/article/180/9/865/2739204
# Smallpox serial interval (18 days) https://academic.oup.com/aje/article/180/9/865/2739204
# Measles serial interval (12 days) https://academic.oup.com/aje/article/180/9/865/2739204
# renderFromR0 <- data.frame(
#   name = c('MERS', 'Influenza', 'Covid-19 (Unmitigated)', 'Ebola', 'SARS', 'Mumps', 'Rubella', 'Smallpox', 'Measles'),
#   fileName = c('MERS.gif', 'Influenza.gif', 'Covid-19.gif', 'Ebola.gif', 'SARS.gif', 'Mumps.gif', 'Rubella.gif', 'Smallpox.gif', 'Measles.gif'),
#   serialInterval = c(13, 3, 5, 15, 8, 20, 18, 18, 12),
#   recoveryInterval = c(14, 14, 14, 14, 14, 14, 14, 14, 14),
#   CFR = c(0.34, 0.001, 0.023, 0.50, 0.10, 0.01, 0.001, 0.30, 0.002),
#   R0 = c(0.8, 1.3, 2.5, 2, 1.85, 5.5, 6, 6, 15)
# );

renderFromR0 <- data.frame(
   name = c('Covid-19 (Low Estimates)', 'Covid-19 (High Estimates)', 'Covid-19 (Unmitigated)', 'Covid-19 (R = ¾ R0)', 'Covid-19 (R = 62.5% R0)', 'Covid-19 (R = ½ R0)'),
   fileName = c('Covid-19-low.gif', 'Covid-19-high.gif', 'Covid-19-mid.gif', 'Covid-19-75.gif', 'Covid-19-625.gif', 'Covid-19-half.gif'),
   serialInterval = c(5, 5, 5, 5, 5, 5),
   recoveryInterval = c(14, 14, 14, 14, 14, 14),
   CFR = c(0.003, 0.034, 0.023, 0.01, 0.01, 0.005),
   R0 = c(1.4, 3.28, 2.5, 1.875, 1.56, 1.25)
);

apply(
  renderFromR0[c(6),],
  1,
  function(d) epidemic(
    size=51,
    iterations=50, # ~ generations, if serialInterval == standardInterval
    interval=0.5,
    width=650,
    file=(d[['fileName']]),
    name=(d[['name']]),
    R0=(as.numeric(d[['R0']])),
    CFR=(as.numeric(d[['CFR']])),
    serialInterval=(as.numeric(d[['serialInterval']])),
    recoveryInterval=(as.numeric(d[['recoveryInterval']]))
  )
)