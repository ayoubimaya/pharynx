proc format;
invalue bmtnum 'NO NODE' = 0  'SINGLE NODE <3CM' = 1  
'NODE 3+CM' = 2 'MULTIPLE NODES'=3;
value   bmtfmt 0 = 'NO NODE'  1 = 'SINGLE NODE <3CM'  
2 = 'NODE 3+CM'  3 = 'MULTIPLE NODES';

invalue bmtnum_ 'TUMOR<2CM' = 1  '2CM<TUMOR<4CM' = 2  
'TUMOR 4+CM' = 3 'MASSIVE TUMOR'=4;
value   bmtfmt_ 1 = 'TUMOR<2CM'  2 = '2CM<TUMOR<4CM'  
3 = 'TUMOR 4+CM'  4 = 'MASSIVE TUMOR';
run;


proc lifetest data=cancer NELSON PLOTS=(S loglogs logsurv);
    time time*status(0);
    strata tx / test=all;
    title "Survival Functions for Cancer Data";
run;

*overlay plot of survival function with confidence bands;
ods select survivalplot;
proc lifetest data=cancer PLOTS=survival(cb=hw test);
    time time*status(0);
    strata tx / test=logrank;
    title "Survival Functions with confidence for Cancer Data";
run;
title;


proc phreg data=cancer simple; *simple displays the simple descriptive statistics for each predictor variable in the model statement;
    model time*status(0)=SEX TX GRADE AGE COND SITE T_STAGE N_STAGE/ ties=exact rl=pl type3(lr); 
	*Using the exact method for handling ties, RL gives CIs for hazard ratios of main effects not involved in interactions, 
		Type3(lr) gives the Likelihood Ratio Type 3 test for each effect that is specified in the MODEL statement (default is Wald);
    title "Cox Proportional Hazards Model of cancer Data";
run;
title;

*Use assess to investigate PH assumptions;
ods graphics on;
ods select cumulativeresiduals scoreprocess;
proc phreg data=cancer;
    model time*status(0)=SEX TX GRADE AGE COND SITE T_STAGE N_STAGE / ties=exact;
    assess  ph / resample seed=90210;
	title 'Assessing PH Assumptions for Main Effects';
run;

*Interactions;
*PL Ratio Test comparing main efffects model and the interactions model;
ods graphics off;
ods output globaltests(match_all persist=proc)=lrtest;
ods select globaltests;
proc phreg data=cancer;
    model time*status(0)=COND  T_STAGE N_STAGE / ties=exact; *main effects only;
	Title 'Main Effects Model';
run;

ods select globaltests;
proc phreg data=cancer;
    model time*status(0)=COND|T_STAGE|N_STAGE @2 / ties=exact; *2 way interactions;
	Title 'Main Effects + Interactions Model';
run;

ods output close;

Title;
proc print data=lrtest;
	Title 'Main Effects Model';
run;
proc print data=lrtest1;
	Title 'Main Effects + Interactions Model';
run;

data testmod;
    merge lrtest(rename=(chisq=chisq1 df=df1))
          lrtest1(rename=(chisq=chisq2 df=df2));
    if test = 'Likelihood Ratio';
    df=df2-df1;
    chisq=chisq2-chisq1;
    pvalue=1-probchi(chisq,df);
	Title 'LR Test Main vs. Interaction';
run;

proc print data=testmod noobs;
    var chisq df pvalue;
    format pvalue pvalue.;
    title 'LRT for Interactions';
run;



*Goodness of Fit;
ods graphics off;
proc phreg data=cancer noprint;
    model TIME*STATUS(0) = COND  T_STAGE N_STAGE / ties=exact;
    output out=fit xbeta=score;
	Title 'Goodness of Fit';
run;

proc rank data=fit groups=10 out=ranks;
    var score;
    ranks bin;
run;

ods select type3;
proc phreg data=ranks;
    model TIME*STATUS(0) =COND T_STAGE N_STAGE / ties=exact type3(score);
    title 'Overall Goodness-of-Fit Test';
run;

title;


*Let's now investigate looking at the residuals;
ods graphics off;
proc phreg data=cancer noprint;
    model time*status(0)=COND T_STAGE N_STAGE / ties=exact;
    output out=resid resdev=deviance 
           ressco=scorecond scoret scoren
           ld=likedisp lmax=maxdisp survival=surv 
           dfbeta=dfcond dft dfn;
run;


data residuals;
    set resid;
    tx_status=compress(tx||status);
run;

proc sgplot data=residuals;
    scatter y=scorecond x=cond / group=tx_status;
    yaxis label='Score Residual';
    title 'Score Residuals for condition';
run;


proc sgplot data=residuals;
    scatter y=scoren x=N_STAGE / group=tx_status;
    yaxis label='Score Residual';
    title 'Score Residuals for n_sstage';
run;


proc sgplot data=residuals;
    scatter y=scoret x=t_STAGE / group=tx_status;
    yaxis label='Score Residual';
    title 'Score Residuals for T_STAGE';
run;

proc sgplot data=residuals;
    scatter y=maxdisp x=surv / group=tx_status;
    yaxis label='L-Max';
   * format tx_status censor.;
    title 'MAX Statistic by Survival Probabilities';
run;

proc sgplot data=residuals;
    scatter y=likedisp x=surv / group=tx_status;
    yaxis label='Likelihood Displacement';
   * format tx_status $censor.;
    title 'LD Statistic by Survival Probabilities';
run;






/*NOW FOR STATIFIED*/

*Stratified Cox Model;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE TX*COND TX*T_STAGE TX*N_STAGE/ ties=exact rl=pl;
    strata TX;
    title 'Stratified Cox Model for CANCER Data';
run;

*Stratified Cox Model dropping TX*N_STAGE;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE TX*COND TX*T_STAGE / ties=exact rl=pl;
    strata TX;
    title '2nd Stratified Cox Model for CANCER Data';
run;

*Stratified Cox Model dropping TX*T_STAGE ;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE TX*COND / ties=exact rl=pl;
    strata TX;
    title '3RD Stratified Cox Model for CANCER Data';
run;

*Stratified Cox Model dropping clinic*prison;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE/ ties=exact rl=pl;
    strata TX;
    title 'Stratified Cox Model for CANCER Data (Main Effects Only)';
run;

*output adjusted survival function data;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE/ ties=exact rl=pl;
    strata TX;
	baseline out=Adj survival=s lower=lcl_s upper=ucl_s; 
    title 'Stratified Cox Model';
run;
*Plotting the data;
proc sgplot data=Adj; 
    band x=TIME upper=ucl_s lower=lcl_s / group=TX;
    step x=TIME y=s /group=TX;
    title 'Stratified Cox Model';
run;


/*Code for creating overlaid survival plots beyond SAS 9.3.*/
ods graphics;
proc phreg data=CANCER plots(overlay cl)=s;
      model TIME*status(0)=COND N_STAGE T_STAGE
        / ties=exact rl=pl;
    strata TX;
    baseline / rowid=TX;
    title 'Stratified Cox Model for CANCER Data';
run;
title;

*Fitting CREATED Time Dependent Variables;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE TIME*N_STAGE / ties=exact;
	Title 'Time Interaction with N_STAGE';
run;

*We see that the PH assumption for clinic is still violated.  Let's estimate the HR for clinic at 1, 2, 3, 4 AND 5 YEARS;
proc phreg data=CANCER;
	model TIME*status(0)=COND N_STAGE T_STAGE
                         TIME*N_STAGE / ties=exact;
    hazardratio N_STAGE / at (TIME=365 730 1095 1460 1825);
    title 'Model N_STAGE Using Time Interaction';
run;


*Piecewise Cox Model;
ods graphics off;
proc phreg data=CANCER;
    model TIME*status(0)=COND N_STAGE T_STAGE
		N_INT1 N_INT2 N_INT3 N_INT4 N_INT5 / ties=exact rl=pl;
    N1=(N_STAGE=1);
    N_INT1=N1*(TIME lt 365);
    N_INT2=N1*(365 le TIME lt 730);
    N_INT3=N1*(730 le TIME lt 1095);
	N_INT4=N1*(1095 le TIME lt 1460);
    N_INT5=N1*(TIME ge 1460);
    PH:test N_INT1=N_INT2=N_INT3=N_INT4=N_INT5;
    title 'Piecewise Cox Model';
run;



/*Inspecting dfbetas to assess influence of observations 
on individual regression coefficients*/
proc phreg data = cancer;
model TIME*status(0)=COND N_STAGE T_STAGE;
strata tx;
output out = dfbeta dfbeta=dfcond dfn dft ;
run;

proc sgplot data = dfbeta;
scatter x = cond y=dfcond / markerchar=case;
run;
proc sgplot data = dfbeta;
scatter x = N_STAGE y=dfn / markerchar=case;
run;
proc sgplot data = dfbeta;
scatter x = T_STAGE y=dft / markerchar=case;
run;


proc print data = CANCER(where=(CASE=141 or CASE=159));
var TX TIME COND N_STAGE T_STAGE;
run;


proc phreg data = CANCER(where=(CASE^=141 and CASE^=159));
model TIME*status(0)=COND N_STAGE T_STAGE;
strata tx;
output out = dfbeta dfbeta=dfcond dfn dft;
run;


/*liklihood displacement*/
proc phreg data = cancer;
model TIME*status(0)=COND N_STAGE T_STAGE;
strata tx;
output out=ld ld=ld;
run;

proc sgplot data=ld;
scatter x=time y=ld / markerchar=case;
run;


/*case 141 has a lot of influence on the model. */

DATA CANCER2;
SET CANCER;
IF CASE=141 THEN DELETE;
RUN;

proc phreg data=cancer2;
    model time*status(0)=SEX GRADE AGE COND SITE T_STAGE N_STAGE / ties=exact SELECTION=BACKWARD;
	STRATA TX;
run;

ods graphics on;
ods select cumulativeresiduals scoreprocess;
proc phreg data=cancer2;
    model time*status(0)=COND  T_STAGE  / ties=exact;
	STRATA TX;
    assess  ph / resample seed=90210;
	title 'Assessing PH Assumptions for Main Effects';
run;
