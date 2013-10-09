/**********************************************************************************************************************
* Author: Alex Keil
* Program: arsen_gest_01.sas
* Date: 7/16/2012
* Project: Copper smelter cohort - arsenic exposed workers
* Tasks: Apply g-estimation to cohort
* Data in: smdat01.sas7bdat, smdat02.sas7bdat
* Data out:
* Description: 
**********************************************************************************************************************/
*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';

LIBNAME dat "E:/School/Outside Projects/CopperArsenic/MSERRdevel/datlib/";

DATA a;
	MERGE dat.lubin_exp dat.lubin_ind;
	BY smid;
	IF age_exit > age >= age_entry;
	in = age;
	out = min(age+1, age+(age_exit-age));
	IF last.smid THEN d = respcan;
	logpy = log(out-in);

	*calendar year;
RUN;

PROC PRINT DATA = a (OBS=10);
	TITLE "Merged data - Lubin's variables";
PROC CONTENTS;
RUN;

*POISSON MODEL;
PROC GENMOD DATA = a;
	TITLE "Crude poisson rate model - respiratory cancer";
	MODEL respcan = d_lubin / LINK=log D=p OFFSET=logpy;
	ESTIMATE "Rate ratio" INTERCEPT 0 D_lubin 1;
RUN;

PROC GENMOD DATA = a;
	TITLE "Crude poisson rate model - respiratory cancer";
	MODEL respcan = d_lubin / LINK=log D=p OFFSET=logpy;
	ESTIMATE "Rate ratio" INTERCEPT 0 D_lubin 1;
RUN;

%MACRO gestit(minpsi=-1, maxpsi=1, step=0.1, eof=1825, idvar=id, data=d, 
			estim=deltanaught, out=gestout, plot=1, gpath=&gpath, opts= NONOTES NOMPRINT NOSYMBOLGEN;
);
OPTIONS &OPTS;
ODS LISTING CLOSE;
%LET iter = 1;
%DO psi_step = %SYSEVALF(&MINPSI*1000) %TO %SYSEVALF(&MAXPSI*1000) %BY %SYSEVALF(&STEP*1000);
%LET psi = %SYSEVALF(&PSI_STEP/1000);
DATA _gestdata;
	SET &data;
	BY &IDVAR;
	
	RETAIN cumdaysused;
	IF first.&idvar THEN cumdaysused=0;
	daysused = exp(&PSI*d_lubin);
	cumdaysused = cumdaysused + daysused; 
	
	%IF %EVAL(&PSI_step=0) %THEN cnaught = 1*&EOF;
	%ELSE %IF %EVAL(&PSI_step>0) %THEN cnaught = exp(0*&PSI)*&EOF;
	%ELSE %IF %EVAL(&PSI_step<0) %THEN cnaught = exp(1*&PSI)*&EOF;
;
	IF last.&idvar THEN DO;
		deltanaught = (cumdaysused < cnaught)*d;
		xnaught = MIN(cumdaysused, cnaught);
		tnaught = cumdaysused;
		xdelta = xnaught*deltanaught;
	END;

PROC MEANS DATA = _gestdata MAX NOPRINT;
	CLASS ID;
	VAR deltanaught xnaught tnaught xdelta;
	OUTPUT OUT= endvars MAX = deltanaught xnaught tnaught xdelta;
DATA _gestdata;
	MERGE _gestdata(DROP=deltanaught xnaught tnaught xdelta) ENDVARS (KEEP = ID deltanaught xnaught tnaught xdelta);
	BY &idvar;
RUN;
*todo - model for continuous exposure!;
PROC PHREG DATA = _gestdata COVS ;
	ID &idvar;
	MODEL (yesterday day)*d_lubin(0) = &ESTIM /*put covariates here*/;
	ODS OUTPUT parameterestimates = PARMS;
RUN;

DATA parms(KEEP=psi estimate stder:);
	SET parms;
	IF parameter = "&ESTIM";
	psi = &PSI;
DATA sampsize;
	SET _gestdata;
	BY id day;
	IF last.id;
RUN;
DATA sampsize (KEEP=ndead);
	SET sampsize END=eof;
	ndead + deltanaught;
	IF eof THEN OUTPUT;
RUN;
DATA parms;
	SET parms;
	IF _N_=1 THEN set sampsize;
RUN;
%IF &ITER=1 %THEN %DO;
	DATA temp; SET parms;
%END;
%ELSE %DO;
	DATA temp; SET temp parms;
%END;
%PUT iteration &ITER of %SYSEVALF((&maxpsi-&minpsi)/&step);
%LET iter = %EVAL(&ITER+1)  ;
%END;
DATA temp;
	SET temp;
	Z = estimate/stderr;
	Zsquared = z**2;
PROC MEANS DATA = temp MIN NOPRINT;
	VAR zsquared;
	OUTPUT OUT=_chimin(DROP=_TYPE_ _FREQ_) MIN=_chimin;

DATA &OUT;
	SET temp;
	IF _N_=1 THEN SET _chimin;
	min=0;
	IF zsquared=_chimin THEN DO;
		min=1;
		CALL symput("chimin", psi);
	END;
ODS LISTING;
%IF &PLOT=1 %THEN %DO;
	ODS LISTING  GPATH = "&gpath" STYLE=MINIMAL;
	ODS GRAPHICS / IMAGENAME="Psi_function_&estim._SAS" RESET=index IMAGEFMT=PNG;
	ODS ESCAPECHAR = '^';
 
	TITLE "G-function for SNM - &ESTIM estimator";
	PROC SGPLOT DATA = &OUT;
		SERIES X=psi Y = Z;
		REFLINE 0 / AXIS=Y LINEATTRS=(COLOR="BLUE");
		REFLINE &chimin / AXIS=X LINEATTRS=(COLOR="BLUE") LABEL="^{unicode psi}^{unicode hat} = &chimin" ;
		XAXIS LABEL="^{unicode psi}^{unicode bar}";
		YAXIS LABEL = "Wald Z";
	PROC SGPLOT DATA = &OUT;
		SERIES X=psi Y = Zsquared;
		REFLINE 0 / AXIS=Y LINEATTRS=(COLOR="BLUE");
		REFLINE &chimin / AXIS=X LINEATTRS=(COLOR="BLUE") LABEL="^{unicode psi}^{unicode hat} = &chimin" ;
		XAXIS LABEL="^{unicode psi}^{unicode bar}";
		YAXIS LABEL = "Wald Z^{unicode '00B2'x}";
%END;
RUN;
OPTIONS NOTES;
%MEND;

LIBNAME dout  "&dpath";

%GESTIT(minpsi=-3, maxpsi=3, step=.01, eof=1825, idvar=id, data=d, 
	estim=xnaught, out=dout.gestoutx, plot=1, gpath=&gpath, opts=nonotes nomprint nosymbolgen);

RUN;QUIT;RUN;
