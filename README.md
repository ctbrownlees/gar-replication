
# [Backtesting Global Growth-at-Risk Replication Files](https://ssrn.com/abstract=3461214)

#### Authors 
 [Christian Brownlees](http://www.econ.upf.edu/~cbrownlees/) and [Andre B.M. Souza](www.andrebmsouza.com)

#### Software Requirements

[MATLAB](https://www.mathworks.com/) The code has been tested with the MATLAB releases R2017a and R2019a

#### Instructions

To replicate the out-of-sample results of Brownlees and Souza (2019) run the script <tt>gar_replication.m</tt>.
The script will create Tables 4 to 6. The tables will be stored as individual CSV files in the directory <tt>tables</tt>.

#### Data

The data used in this study was downloaded from the following sources in June 2019.

 - GDP: 'https://stats.oecd.org/sdmx-json/data/DP_LIVE/.QGDP.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en'
 - NFCI: 'https://www.imf.org/~/media/Files/Publications/GFSR/2017/October/chapter-3/csv-data/data-appendix.ashx?la=eni'
 - ST_INTEREST: 'https://stats.oecd.org/sdmx-json/data/DP_LIVE/.STINT.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en'
 - LT_INTEREST: 'https://stats.oecd.org/sdmx-json/data/DP_LIVE/.LTINT.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en'
 - BAA_CREDIT:  'https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=BAA10Y&scale=left&cosd=1986-01-02&coed=2019-10-15&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-10-17&revision_date=2019-10-17&nd=1986-01-02'
 - AAA_CREDIT:  'https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=AAA10Y&scale=left&cosd=1983-01-03&coed=2019-10-15&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-10-17&revision_date=2019-10-17&nd=1983-01-03'
 - CREDIT_STATISTICS: 'https://www.bis.org/statistics/c_gaps/c_gaps.xlsx'
 - PROPERTY PRICES: 'https://www.bis.org/statistics/pp/pp_detailed.xlsx'
 - WUI: 'https://www.policyuncertainty.com/media/WUI_Data.xlsx'
 - GPR: 'https://www2.bc.edu/matteo-iacoviello/gpr_files/gpr_web_latest.xlsx'
 - EPU: 
   * 'https://www.policyuncertainty.com/media/Australia_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Canada_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Europe_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Europe_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Europe_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Mexico_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Korea_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/UK_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Japan_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/HKKS_Greece_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Ireland_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Netherlands_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Canada_Policy_Uncertainty_Data.xlsx'
   * 'https://www.policyuncertainty.com/media/Sweden_Policy_Uncertainty_Data.xlsx'

