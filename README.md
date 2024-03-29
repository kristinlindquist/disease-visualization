# A simple disease model
A simple disease spread visualization based on R0, CFR, serial interval and the final size relation. [Associated Medium post here](https://medium.com/@k_a_lindquist/reproduction-numbers-fatality-rates-6660a1e9d933).

Assumes a completely isolated, small community in which [susceptible depletion](https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2016.0659) kicks in, using the final size relation s∞ = S(0)e^[–R0(1–s∞)] [[1](https://mathematicsinindustry.springeropen.com/track/pdf/10.1186/s13362-019-0058-7) [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/)]. Other assumptions include homogeneous contact & transmission rates (whereas better models consider heterogeneous contact & transmission rates, e.g. [[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4808916/), [2](https://www.uvm.edu/pdodds/research/papers/years/2005/watts2005a.pdf), [3](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0120701)]), that [serial interval == generation interval](https://nccid.ca/publications/glossary-terms-infectious-disease-modelling-proposal-consistent-language/) and that no behavior change occurs over the course of the epidemic. If you see any errors or have questions/concerns/comments, please [hit me up](mailto:kristin@rationally.io).

## Example (Smallpox)
![Smallpox](https://github.com/kristinlindquist/disease-visualization/assets/9382486/ce6440bf-4d83-42d7-a1bd-b28972bd2a25)

## Sources for parameters
### R0
* [Various](https://en.wikipedia.org/wiki/Basic_reproduction_number)
* [SARS (1.77-1.85)](https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1)
* [Seasonal Influenza (1.3)](https://www.ncbi.nlm.nih.gov/pubmed/19545404)
* [Covid-19 (1.4-2.5)](https://www.who.int/news-room/detail/23-01-2020-statement-on-the-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov))
* [Covid-19 high (mean 3.28, median 2.79)](https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319)
* [Covid-19 (1.32 with a corresponding SI of 3.96)](https://www.jwatch.org/na51171/2020/03/27/serial-interval-covid-19)

### Fatality Rate
* [Covid-19 high 3.4%](https://www.who.int/dg/speeches/detail/who-director-general-s-opening-remarks-at-the-media-briefing-on-covid-19---3-march-2020)
* [Covid-19 CFR high (2.3%)](https://ourworldindata.org/coronavirus)
* [Covid-19 IFR (1.6%)](https://www.medrxiv.org/content/10.1101/2020.03.04.20031104v1.full.pdf)
* [Covid-19 IFR (0.94%)](https://institutefordiseasemodeling.github.io/nCoV-public/analyses/first_adjusted_mortality_estimates_and_risk_assessment/2019-nCoV-preliminary_age_and_time_adjusted_mortality_rates_and_pandemic_risk_assessment.html)
* [Covid-19 IFR low (0.3%–1%)](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200219-sitrep-30-covid-19.pdf?sfvrsn=3346b04f_2)
* [Ebola, SARS and MERS CFR (50%, 10%, 34%)](https://ourworldindata.org/coronavirus)
* [Seasonal Influenza CFR (0.1%)](https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates)
* [Smallpox CFR (30%)](https://en.wikipedia.org/wiki/Smallpox)
* [Measles CFR (0.2%)](https://www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf)
* [Mumps CFR (1%)](https://en.wikipedia.org/wiki/List_of_human_disease_case_fatality_rates)
* [Rubella CFR (0.1% - fatalities almost entirely infants & in utero)](https://www.cdc.gov/rubella/about/in-the-us.html)

### Serial Interval
* [MERS (13 days)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930778/)
* [Seasonal Influenza (3 days)](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200306-sitrep-46-covid-19.pdf?sfvrsn=96b04adf_2)
* [Covid-19 (5-6 days)](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200306-sitrep-46-covid-19.pdf?sfvrsn=96b04adf_2)
* [Covid-19 (3.96 with correspond R0 of 1.32)](https://www.jwatch.org/na51171/2020/03/27/serial-interval-covid-19)
* [Ebola (15.3 days)](https://www.sciencedirect.com/science/article/pii/S1755436515000341)
* [SARS (8.4 days)](https://dash.harvard.edu/bitstream/handle/1/25620506/Transmission%20dynamics%20and%20control%20of%20severe%20acute%20respiratory%20syndrome.pdf?sequence=1)
* [Mumps (20 days)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223546/)
* [Rubella (18 days)](https://academic.oup.com/aje/article/180/9/865/2739204)
* [Smallpox (18 days)](https://academic.oup.com/aje/article/180/9/865/2739204)
* [Measles (12 days)](https://academic.oup.com/aje/article/180/9/865/2739204)
