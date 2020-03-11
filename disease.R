library('animation');
library('spatstat');
library('sp');
library('ggplot2');
library('rootSolve');

S0 <- 0.95;
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
getDotSize <- function(size, width) width / (size * 4);

# http://coleoguy.blogspot.com/2016/04/stochasticprobabilistic-rounding.html
getStochRound <- function(x) {
  q = abs(x - trunc(x));
  adj = sample(0:1, size = 1, prob = c(1 - q, q));
  if(x < 0) adj = adj * -1;
  trunc(x) + adj;
}

epidemic <- function(size, generations, filename, name, R0, CFR, width = 600) {
  df <- data.frame(expand.grid(c(list(x = 1:size, y = 1:size, status = 0))));
  df[df$x == ceiling(size/2) & df$y == ceiling(size/2),]$status = 1;
  df$distance <- spDistsN1(as.matrix(df[c("x", "y")]), c(ceiling(size/2), ceiling(size/2)));
  df <- df[order(df$distance),];
  df[1:initialInfected,]$status <- 1;
  finalR <- getEffectiveR(df, R0);
  saveGIF({
    for (i in 1:generations) {
      R <- getEffectiveR(df, R0);
      print(ggplot(df, aes(x=x, y=y, fill=getColors(df$status))) + 
        geom_point(color=getColors(df$status), size = getDotSize(size, width)) +
        labs(
          title=name,
          subtitle=paste('Generation ', i, ', R0=', R0, ', CFR=', CFR*100, '%', sep='')
        ) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(colour = 'gray10', family='Anonymous Pro', hjust = 0.5, size=35, vjust = 0),
          plot.subtitle = element_text(colour = 'gray10', family='Anonymous Pro', hjust = 0.5, size=25, vjust = 0),
          plot.margin=unit(c(1.1, 0.25, 1, 0.25), "cm")
        )
      );
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
          copy[copy$x == d[['x']] & copy$y == d[['y']],]$status <<- sample(3:4, 1, prob=c(1-CFR, CFR));
        }
      });
      df <- copy;
    }
  }, movie.name=filename, interval = 2, ani.width = width, ani.height = width);
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
renderFromR0 <- data.frame(
  name = c('MERS', 'Influenza', 'Covid-19', 'Ebola', 'SARS', 'Mumps', 'Rubella', 'Smallpox', 'Measles'),
  fileName = c('MERS.gif', 'Influenza.gif', 'Covid-19.gif', 'Ebola.gif', 'SARS.gif', 'Mumps.gif', 'Rubella.gif', 'Smallpox.gif', 'Measles.gif'),
  CFR = c(0.34, 0.001, 0.023, 0.50, 0.10, 0.01, 0.001, 0.30, 0.002),
  R0 = c(0.8, 1.3, 2.5, 2, 1.85, 5.5, 6, 6, 15)
);

# renderFromR0 <- data.frame(
#    name = c('Covid-19 - Low', 'Covid-19 - Medium', 'Covid-19 - High'),
#    fileName = c('Covid-19-low.gif', 'Covid-19-medium.gif', 'Covid-19-high.gif'),
#    CFR = c(0.01, 0.023, 0.034),
#    R0 = c(1.4, 2.5, 3.28)
# );

apply(
  renderFromR0[2,],
  1,
  function(d) epidemic(
    size=31,
    generations=15,
    width=600,
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