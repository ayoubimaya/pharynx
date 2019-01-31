<!-- Readme_HTML for pharynx data -->

<h1>Background</h1>
<p> It is well known that tobacco and alcohol can have adverse effects on human health. 
	About 90% of people with <a href="https://en.wikipedia.org/wiki/Oral_cancer">oral cancer</a> used tobacco products, and 75-80% of those patients also consumed alcohol. 
	According to the <a href="https://www.cdc.gov/">CDC</a>, the survival rate is approximately 50%. 
	</p>

<h1>About Pharynx Analysis (R)</h1>
<p>The purpose of the study done in R is to assess what variables affect the survival rate of carcinoma of the oropharynx. 
<a href="https://en.wikipedia.org/wiki/Oral_cancer">Oral cancer</a> is hard to discover and diagnose, meaning that late stage discovery is more common. 
Due to late stage discovery, this cancer has often metastasized to another location thus reducing the time duration until death. 
A survival analysis on <a href="https://en.wikipedia.org/wiki/Oral_cancer">oral cancer</a> helps to answer what characteristics increase or decrease the probability of survival, what the proportion of the population that will survive past a given time, and what will be the rate at which those who survive die or fail.</p>

<h2>Methods (R)</h2>
<p>I chose to analyze the data using <a href="https://en.wikipedia.org/wiki/Proportional_hazards_model">Cox Proportional Hazards model</a> due to its advantage of fitting survival models without specifying the underlying hazard function. 
This is helpful, because we don't necessarily always know what the underlying hazard function is. 
The predictors I used to look at were condition (COND), tumor stage (T_STAGE), and node stage (N_STAGE). 
When running a <a href="https://en.wikipedia.org/wiki/Proportional_hazards_model">Cox Proportional Hazards model</a>, I noticed that node stage violated the proportional hazards assumption. 
I then decided to create an initial regression model where time (TIME) was my independent variable and condition and tumor stage as my dependent variables in order to look at an outlier test. 
I found that one observation highly influenced the model and decided to remove that outlier. 
Once removing that outlier I decided to start over with another <a href="https://en.wikipedia.org/wiki/Proportional_hazards_model">Cox Proportional Hazards model</a>, however node stage was still an insignificant predictor. 
I also tried to introduce other predictors as well as interactions, but it appears that tumor stage and condition were the only significant predictors. </p>

<h1>About Pharynx Analysis (SAS)</h1>
<p>The purpose of the study was to predict patient outcome of those diagnosed with carcinoma of the oropharynx. 
	I hypothesize that treatment type will have a significant effect on survival. </p>

<h2>Methods (SAS)</h2>
    <p>I collected data starting in January 1968 and continued until December 1972. 
	To be included in the study, I looked at patient eligibility criteria, and screened for extreme cases. 
	The study consists of 195 patients diagnosed with <a href="https://en.wikipedia.org/wiki/Oral_cancer">oral cancer</a> in varying stages of the cancer.
	The data was collected at 6 large cancer research institutions across the United States. 
	The patients were randomly assigned to either the standard treatment, radiation therapy, or the test treatment, radiation coupled with chemotherapy. 
	Approximately 27% of participants were censored due to the patient moving or transferring to another institution that did not participate in the study.
    I chose to analyze the data using <a href="https://en.wikipedia.org/wiki/Proportional_hazards_model">Cox Proportional Hazards model</a> as the data suggested. 
	To build our model I used a backward selection method to analyze the most relevant variables. 
	I found that condition, tumor stage and node stage were the most relevant. 
	I found that node stage violates the proportional hazard assumption. 
	After additional analysis, I found that we had an observation that heavily influenced the model. 
	I decided to remove the observation, as it may have been entered wrong, because the values did not make sense intuitively.
	<br>Analysis was completed using SAS 9.4.  </p>

	
<h1>Where to find the data</h1>
<p>Data for Carcinoma of the Oropharynx can be found <a href="http://www.umass.edu/statdata/statdata/data/">here</a>.
</p>
<p></p>
