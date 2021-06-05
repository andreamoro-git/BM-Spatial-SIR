* generates table with estimates and "true" parameters when there are
* multiple treatment times

local basedir = "output/nc5-"
local filetype = "dens"
capture log close
log using `basedir'sim_estimate_policies_multitimes.log, text replace

insheet using `basedir'`filetype'-times-25pc.csv, comma clear
set seed 2443

// drop observations before peak or not
*local drop " if t>tm"
local drop "if _n<0"


quietly xtset naics t
set more off
*gen treated = t>=20 & npi_date==20
replace outside = 1 if t==0
summ

sort naics t
gen growthi = D.F1.active/active
tab density, gen(ddens)

egen maxactive = max(active), by(naics)
gen tmax = t if maxactive==active
egen tm = max(tmax), by(naics)
drop `drop'

tab density
tab density npi_date

*generate a "fake" sample with one half cities
generate random = runiform()
egen randcity = mean(random), by(naics)

*gen mysample = randcity>.5
*tab density npi_date if mysample==1

keep if density==1.5 | density==0.5 | density==1
keep if npi_date== 15 | npi_date==40 | npi_date==999

gen mysample = 0
replace mysample = 1 if npi_date==999 & density==1
replace mysample = 1 if npi_date==40 & density==0.5
replace mysample = 1 if npi_date==15 & density==1.5
gen treated = 0
replace treated = 1 if npi_date==40 & t>=40
replace treated = 1 if npi_date==15 & t>=15


xtreg growthi i.t treated##c.density if mysample==1, fe
scalar g0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Growth rate, Model, Baselineint, Specification, Estimated) replace
xtreg growthi i.t treated##c.density, fe
scalar gall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Growth rate, Model, Baselineint, Specification, True) replace
xtreg outside i.t treated##c.density if mysample==1, fe
scalar o0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Contacts, Model, Baselineint, Specification, Estimated) append
xtreg outside i.t treated##c.density, fe
scalar oall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Contacts, Model, Baselineint, Specification, True) append
xtreg active i.t treated##c.density if mysample==1, fe
scalar ac0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Infected, Model, Baselineint, Specification, Estimated) append
xtreg active i.t treated##c.density, fe
scalar acall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Infected, Model, Baselineint, Specification, True) append

preserve
use `basedir'resultsest, clear
keep if var=="1.treated" | var=="1.treated#c.density"
rename coef Estimated
gen Estimated2 = Estimated if var=="1.treated#c.density"
replace Estimated = . if var=="1.treated#c.density"
collapse (max)Estimate*, by(Outcome Model)
save `basedir'resultsest, replace
list
use `basedir'resultstrue, clear
keep if var=="1.treated" | var=="1.treated#c.density"
rename coef True
gen True2 = True if var=="1.treated#c.density"
replace True = . if var=="1.treated#c.density"
collapse (max)True*, by(Outcome Model)
save `basedir'resultstrue, replace
list
restore

// no interaction
xtreg growthi i.t treated if mysample==1, fe
scalar g0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Growth rate, Model, Baseline, Specification, Estimated) append
xtreg growthi i.t treated, fe
scalar gall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Growth rate, Model, Baseline, Specification, True) append
xtreg outside i.t treated if mysample==1, fe
scalar o0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Contacts, Model, Baseline, Specification, Estimated) append
xtreg outside i.t treated, fe
scalar oall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Contacts, Model, Baseline, Specification, True) append
xtreg active i.t treated if mysample==1, fe
scalar ac0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Infected, Model, Baseline, Specification, Estimated) append
xtreg active i.t treated, fe
scalar acall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Infected, Model, Baseline, Specification, True) append

use `basedir'resultsest, clear
keep if var=="treated" | Model=="Baselineint"
replace Estimated = coef if var=="treated"
save `basedir'resultsest, replace
sort Outcome
list Out Mod Est*

use `basedir'resultstrue, clear
keep if var=="treated" | Model=="Baselineint"
replace True = coef if var=="treated"
save `basedir'resultstrue, replace
merge 1:1 Outcome Model using `basedir'resultsest
sort Outcome
list Out Mod True* Est*
gen sortvar = 1 if Outcome=="Infected"
replace sortvar = 2 if Outcome=="Contacts"
replace sortvar = 3 if Outcome=="Growth rate"

drop _merge
sort Outcome Model
save `basedir'allresults, replace

********
* now with behavioral responses
********

local filetype = "dens-beh_p"

insheet using `basedir'`filetype'-times-25pc.csv, comma clear
set seed 2443
// drop observations before peak or not
*local drop " if t>tm"
local drop "if _n<0"


quietly xtset naics t
set more off
*gen treated = t>=20 & npi_date==20
replace outside = 1 if t==0
summ

sort naics t
gen growthi = D.F1.active/active
tab density, gen(ddens)

egen maxactive = max(active), by(naics)
gen tmax = t if maxactive==active
egen tm = max(tmax), by(naics)
drop `drop'

tab density
tab density npi_date

*generate a "fake" sample with one half cities
generate random = runiform()
egen randcity = mean(random), by(naics)

*gen mysample = randcity>.5
*tab density npi_date if mysample==1

keep if density==1.5 | density==0.5 | density==1
keep if npi_date== 15 | npi_date==40 | npi_date==999

gen mysample = 0
replace mysample = 1 if npi_date==999 & density==1
replace mysample = 1 if npi_date==40 & density==0.5
replace mysample = 1 if npi_date==15 & density==1.5
gen treated = 0
replace treated = 1 if npi_date==40 & t>=40
replace treated = 1 if npi_date==15 & t>=15


xtreg growthi i.t treated##c.density if mysample==1, fe
scalar g0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Growth rate, Model, Beh_Baselineint, Specification, Estimated) replace
xtreg growthi i.t treated##c.density, fe
scalar gall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Growth rate, Model, Beh_Baselineint, Specification, True) replace
xtreg outside i.t treated##c.density if mysample==1, fe
scalar o0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Contacts, Model, Beh_Baselineint, Specification, Estimated) append
xtreg outside i.t treated##c.density, fe
scalar oall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Contacts, Model, Beh_Baselineint, Specification, True) append
xtreg active i.t treated##c.density if mysample==1, fe
scalar ac0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Infected, Model, Beh_Baselineint, Specification, Estimated) append
xtreg active i.t treated##c.density, fe
scalar acall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Infected, Model, Beh_Baselineint, Specification, True) append

preserve
use `basedir'resultsest, clear
keep if var=="1.treated" | var=="1.treated#c.density"
rename coef Estimated
gen Estimated2 = Estimated if var=="1.treated#c.density"
replace Estimated = . if var=="1.treated#c.density"
collapse (max)Estimate*, by(Outcome Model)
save `basedir'resultsest, replace
list
use `basedir'resultstrue, clear
keep if var=="1.treated" | var=="1.treated#c.density"
rename coef True
gen True2 = True if var=="1.treated#c.density"
replace True = . if var=="1.treated#c.density"
collapse (max)True*, by(Outcome Model)
save `basedir'resultstrue, replace
list
restore


// no interaction
xtreg growthi i.t treated if mysample==1, fe
scalar g0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Growth rate, Model, Beh_Baseline, Specification, Estimated) append
xtreg growthi i.t treated, fe
scalar gall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Growth rate, Model, Beh_Baseline, Specification, True) append
xtreg outside i.t treated if mysample==1, fe
scalar o0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Contacts, Model, Beh_Baseline, Specification, Estimated) append
xtreg outside i.t treated, fe
scalar oall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Contacts, Model, Beh_Baseline, Specification, True) append
xtreg active i.t treated if mysample==1, fe
scalar ac0 = e(ll)
regsave using `basedir'resultsest, addlabel(Outcome, Infected, Model, Beh_Baseline, Specification, Estimated) append
xtreg active i.t treated, fe
scalar acall = e(ll)
regsave using `basedir'resultstrue, addlabel(Outcome, Infected, Model, Beh_Baseline, Specification, True) append

use `basedir'resultsest, clear
keep if var=="treated" | Model=="Beh_Baselineint"
replace Estimated = coef if var=="treated"
save `basedir'resultsest, replace
sort Outcome
list Out Mod Est*

use `basedir'resultstrue, clear
keep if var=="treated" | Model=="Beh_Baselineint"
replace True = coef if var=="treated"
save `basedir'resultstrue, replace
merge 1:1 Outcome Model using `basedir'resultsest
sort Outcome
list Out Mod True* Est*
gen sortvar = 1 if Outcome=="Infected"
replace sortvar = 2 if Outcome=="Contacts"
replace sortvar = 3 if Outcome=="Growth rate"

sort Outcome Model
drop _merge

merge Outcome Model using `basedir'allresults
gen modsort= Model=="Beh_Baseline" | Model=="Beh_Baselineint"

replace Outcome = "\multirow{2}{*}{"+Outcome+"}"
replace Outcome = " " if Model=="Beh_Baselineint" | Model=="Baselineint"
local nob = _N+2
set obs `nob'
replace modsort = -1 if _n==_N-1
replace Outcome = " \\ \multicolumn{5}{c}{Without behavioral responses} \\ \midrule " if _n==_N-1
replace modsort = 0.5 if _n==_N
replace Outcome = " \\ \multicolumn{5}{c}{With behavioral responses} \\ \midrule " if _n==_N
sort  modsort sortvar Model

list modsort Model Out True Estimated True2 Estimated2

format True* Estim* %10.3f
list True Estimated True2 Estimated2
listtex Out True Estimated True2 Estimated2 using ../output/estimates-multitimes.tex, replace type rstyle(tabular) head(\begin{tabular}{lcccc} \toprule & \multicolumn{2}{c}{Treated} & \multicolumn{2}{c}{Treated\#Density}  \\ Outcome  & True & Estimated & True & Estimated \\ \midrule \\ ) foot("\bottomrule \end{tabular}")

log close
