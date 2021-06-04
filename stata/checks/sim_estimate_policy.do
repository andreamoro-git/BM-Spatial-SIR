//

// set some filenames and file types
local filetype = "dens"
local nobehfile "../output/nc5-`filetype'-20-80-25pc.csv"
local filetype = "dens-beh_p"
local behfile "../output/nc5-`filetype'-20-80-25pc.csv"
// drop observations before peak or not
*local drop " if t>tm"
local drop "if _n<0"
*local drop "if t>=40"
set seed 2443


// start with no behavior
insheet using `nobehfile', comma clear

// generate some variables
quietly xtset naics t
set more off
gen treated = t>=20 & npi_date==20
replace outside = 1 if t==0
summ

sort naics t
gen pctactive = active/q_popsize
gen growthi = D.F1.active/active
tab density, gen(ddens)
tab dens

egen maxactive = max(active), by(naics)
gen tmax = t if maxactive==active
egen tm = max(tmax), by(naics)

gen check = tm==20
tab check npi_date, col

tab npi_date treated
*list naics t check npi_date treated if treated==1

* drop observations after the peak if requested
drop `drop'

* drop some observations for symmetry
foreach dens in 0.52 0.54 0.56 0.58 1.5 {
  drop if abs(density - `dens') <0.01
}

gen mysample = 1
* selectively drop half cities by density and treatment
foreach dens in 0.5 0.7 0.9 1.1 1.3 {
  replace mysample=0 if npi_date==20 & abs(density-`dens')<0.00001
}
foreach dens in 0.6 0.8 1 1.2 1.4 {
  replace mysample =0 if npi_date==999 & abs(density-`dens')<0.0001
}

xtreg growthi i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, Estimated) replace

predict prgr
summ prgr growthi
summ prgr growthi if mysample==1
sum prgr growthi if mysample==0 & t>20


xtreg growthi i.t treated, fe
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, True) append
xtreg outside i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, Estimated) append
xtreg outside i.t treated, fe
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, True) append
xtreg pctactive i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, "\% Active", Model, Baseline, Specification, Estimated) append
xtreg pctactive i.t treated, fe
regsave using results, addlabel(Outcome, "\% Active", Model, Baseline, Specification, True) append
//xtreg growthi i.t##mysample treated##mysample


/////////////////////////////////////////////////////////////////////
* now with behavior
insheet using `behfile', comma clear
quietly xtset naics t
gen treated = t>=20 & npi_date==20
replace outside = 1 if t==0

summ

sort naics t
gen pctactive = active/q_popsize
gen growthi = D.F1.active/active
tab density, gen(ddens)
tab dens

egen maxactive = max(active), by(naics)
gen tmax = t if maxactive==active
egen tm = max(tmax), by(naics)

drop `drop'

* drop some observations for symmetry
foreach dens in 0.52 0.54 0.56 0.58 1.5 {
  drop if abs(density - `dens') <0.01
}

gen mysample = 1

* selectively drop half cities by density and treatment
foreach dens in 0.5 0.7 0.9 1.1 1.3 {
  replace mysample=0 if npi_date==20 & abs(density-`dens')<0.00001
}
foreach dens in 0.6 0.8 1 1.2 1.4 {
  replace mysample =0 if npi_date==999 & abs(density-`dens')<0.0001
}


gen post = t>=20
xtreg growthi i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, Growth rate, Model, With behavior, Specification, Estimated) append
xtreg growthi i.t treated, fe
regsave using results, addlabel(Outcome, Growth rate, Model, With behavior, Specification, True) append
xtreg outside i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, Contacts, Model, With behavior, Specification, Estimated) append
xtreg outside i.t treated, fe
regsave using results, addlabel(Outcome, Contacts, Model, With behavior, Specification, True) append
xtreg pctactive i.t treated if mysample==1, fe
regsave using results, addlabel(Outcome, "\% Active", Model, With behavior, Specification, Estimated) append
xtreg pctactive i.t treated, fe
regsave using results, addlabel(Outcome, "\% Active", Model, With behavior, Specification, True) append
//xtreg growthi i.t##mysample treated##mysample
//ßxtreg outside i.t##mysample treated##mysample

use results, clear
keep if var=="treated"

foreach spec in "True" "Estimated"{
  gen _`spec' = coef if Spec == "`spec'"
  gen _`spec'sd = stderr if Spec == "`spec'"
  gen _`spec'N = N if Spec=="`spec'"
  egen __`spec' = max(_`spec'), by(Model Outcome)
  egen __`spec'sd = max(_`spec'sd), by(Model Outcome)
  egen `spec'N = max(_`spec'N), by(Model Outcome)
  gen `spec' = string(__`spec',"%5.3f") + " (" + string(__`spec'sd,"%5.3f") + ")"
}

sort Outcome Model
replace Outcome = "\multirow{2}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model Est* True*
*texsave using table.tex, replace nofix preamble("\usepackage{multirow}") align(llcc) width(.8\linewidth)

listtex Outcome Model True TrueN Estimated EstimatedN, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & N & Estimated & N \\ \midrule") foot("\bottomrule \\n \end{tabular}")

list Out Mod
