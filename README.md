# diseasemodel
A simple disease spread model based on R0, CFR and the final size relation

Assumes [susceptible depletion](https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2016.0659) kicks in using the final size relation s∞ = S(0)e^[–R0(1–s∞)] [[1](https://mathematicsinindustry.springeropen.com/track/pdf/10.1186/s13362-019-0058-7) [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/)] 

## R0s
* [Various](https://en.wikipedia.org/wiki/Basic_reproduction_number)
* [SARS (1.77-1.85)](https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1)
* [Seasonal Flu](https://www.ncbi.nlm.nih.gov/pubmed/19545404)
* [Covid-19 (1.4-2.5)](https://www.who.int/news-room/detail/23-01-2020-statement-on-the-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov))
* [Covid-19 high (mean 3.28, median 2.79)](https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319)

## Fatality Rates
* [Covid-19 high 3.4%](https://www.who.int/dg/speeches/detail/who-director-general-s-opening-remarks-at-the-media-briefing-on-covid-19---3-march-2020)
* [Covid-19 CFR high (2.3%)](https://ourworldindata.org/coronavirus)
* [Covid-19 IFR (1.6%)](https://www.medrxiv.org/content/10.1101/2020.03.04.20031104v1.full.pdf)
* [Covid-19 IFR (0.94%)](https://institutefordiseasemodeling.github.io/nCoV-public/analyses/first_adjusted_mortality_estimates_and_risk_assessment/2019-nCoV-preliminary_age_and_time_adjusted_mortality_rates_and_pandemic_risk_assessment.html)
* [Covid-19 IFR low (0.3%–1%)](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200219-sitrep-30-covid-19.pdf?sfvrsn=3346b04f_2)
* [Ebola, SARS and MERS CFR](https://ourworldindata.org/coronavirus)
* [Flu CFR](https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates)
* [Smallpox CFR](https://en.wikipedia.org/wiki/Smallpox)
* [Measles CFR](https://www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf)
* [Mumps CFR](https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates)
* [Rubella CFR (infants & in utero)](https://www.cdc.gov/rubella/about/in-the-us.html)
