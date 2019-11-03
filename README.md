
# [Backtesting Global Growth-at-Risk Replication Files](https://ssrn.com/abstract=3461214)

This repository contains the replication files for the paper <i>Backtesting Global Growth-at-Risk</i>
by Christian Brownlees and Andre B.M. Souza which is available on SSRN at the address
[https://ssrn.com/abstract=3461214](https://ssrn.com/abstract=3461214)

## Authors 
 [Christian Brownlees](http://www.econ.upf.edu/~cbrownlees/) and [Andre B.M. Souza](http://www.andrebmsouza.com)

## Software Requirements

[MATLAB](https://www.mathworks.com/) The code has been tested with the MATLAB releases R2017a and R2019a

## Instructions

To replicate the out-of-sample results run the script <tt>gar_replication.m</tt>.
The script will create Tables 4 to 6 of the paper. The tables will be stored as individual CSV files in the directory <tt>tables</tt>.

## Data

***Important Disclaimer:*** The data used in this study was downloaded from the following sources in June 2019.

 - [GDP](https://stats.oecd.org/sdmx-json/data/DP_LIVE/.QGDP.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en) from the OECD Database
 - [NFCI](https://www.imf.org/~/media/Files/Publications/GFSR/2017/October/chapter-3/csv-data/data-appendix.ashx?la=eni) from the IMF.
 - [ST_INTEREST](https://stats.oecd.org/sdmx-json/data/DP_LIVE/.STINT.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en) from the OECD database
 - [LT_INTEREST](https://stats.oecd.org/sdmx-json/data/DP_LIVE/.LTINT.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en) from the OECD database
 - [BAA_CREDIT](https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=BAA10Y&scale=left&cosd=1986-01-02&coed=2019-10-15&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-10-17&revision_date=2019-10-17&nd=1986-01-02) from the St. Louis Fed.
 - [AAA_CREDIT](https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=AAA10Y&scale=left&cosd=1983-01-03&coed=2019-10-15&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-10-17&revision_date=2019-10-17&nd=1983-01-03) from the St. Louis Fed.
 - [CREDIT_STATISTICS](https://www.bis.org/statistics/c_gaps/c_gaps.xlsx) from the BIS database
 - [PROPERTY PRICES](https://www.bis.org/statistics/pp/pp_detailed.xlsx)  from the BIS database
 - [WUI](https://www.policyuncertainty.com/media/WUI_Data.xlsx) from the Policy Uncertainty website
 - [GPR](https://www2.bc.edu/matteo-iacoviello/gpr_files/gpr_web_latest.xlsx) from the Policy Uncertainty website
 - EPU for several countries, all of which can be found in the Policy Uncertainty website: 
   * [Australia](https://www.policyuncertainty.com/media/Australia_Policy_Uncertainty_Data.xlsx)
   * [Canada](https://www.policyuncertainty.com/media/Canada_Policy_Uncertainty_Data.xlsx)
   * [Europe](https://www.policyuncertainty.com/media/Europe_Policy_Uncertainty_Data.xlsx)
   * [Mexico](https://www.policyuncertainty.com/media/Mexico_Policy_Uncertainty_Data.xlsx)
   * [South Korea](https://www.policyuncertainty.com/media/Korea_Policy_Uncertainty_Data.xlsx)
   * [United Kingdom](https://www.policyuncertainty.com/media/UK_Policy_Uncertainty_Data.xlsx)
   * [United States](https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx)
   * [Japan](https://www.policyuncertainty.com/media/Japan_Policy_Uncertainty_Data.xlsx)
   * [Greece](https://www.policyuncertainty.com/media/HKKS_Greece_Policy_Uncertainty_Data.xlsx)
   * [Ireland](https://www.policyuncertainty.com/media/Ireland_Policy_Uncertainty_Data.xlsx)
   * [Netherlands](https://www.policyuncertainty.com/media/Netherlands_Policy_Uncertainty_Data.xlsx)
   * [Sweden](https://www.policyuncertainty.com/media/Sweden_Policy_Uncertainty_Data.xlsx)

## Additional Resources
### [Vulnerable Growth Replication Files (Adrian et al, 2019)](https://www.aeaweb.org/articles?id=10.1257/aer.20161923)
 - rq.m: Function to compute quantile regression. Source: Vulnerable Growth Replication Files (Adrian et al, 2019)

### [MFE Toolbox](https://github.com/bashtage/mfe-toolbox)
 - covnw.m: Function to estimate HAC covariance matrices.
 - olsnw.m: Function to perform inference on ols parameters with HAC.
 - normloglik.m: Log Likelihood for the standard normal. 

### The following files are modified versions of the above. We use them to perform DQ tests for h>1 step ahead forecasts.
 - covnw_dq.m: Function to estimate HAC covariance matrices for binary hit variables.
 - olsnw_dq.m: Function to perform robust inference on binary hit variables.

### [Belloni, Chernozhukov (2011) code files](https://faculty.fuqua.duke.edu/~abn5/belloni-software.html)
 - L1QR.m: L1 penalized quantile regression. Depends on SDPT3-4.0. See below.

### Packages:
 - [SDPT3-4.0](https://github.com/sqlp/sdpt3): 
 

