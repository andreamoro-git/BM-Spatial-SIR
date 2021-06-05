//

capture log close
log using sim_estimate_policies_multitimes.log, text replace
local filetype = "dens-beh_p"
*local filetype = "dens"

insheet using ../output/nc5-`filetype'-times-25pc.csv, comma clear
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
regsave using results, addlabel(Outcome, Growth rate, Model, Baselineint, Specification, Estimated) replace
xtreg growthi i.t treated##c.density, fe
scalar gall = e(ll)
regsave using results, addlabel(Outcome, Growth rate, Model, Baselineint, Specification, True) append
xtreg outside i.t treated##c.density if mysample==1, fe
scalar o0 = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baselineint, Specification, Estimated) append
xtreg outside i.t treated##c.density, fe
scalar oall = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baselineint, Specification, True) append
xtreg active i.t treated##c.density if mysample==1, fe
scalar ac0 = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baselineint, Specification, Estimated) append
xtreg active i.t treated##c.density, fe
scalar acall = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baselineint, Specification, True) append

// no interaction
xtreg growthi i.t treated if mysample==1, fe
scalar g0 = e(ll)
regsave using resultsnoint, addlabel(Outcome, Growth rate, Model, Baseline, Specification, Estimated) replace
xtreg growthi i.t treated, fe
scalar gall = e(ll)
regsave using resultsnoint, addlabel(Outcome, Growth rate, Model, Baseline, Specification, True) append
xtreg outside i.t treated if mysample==1, fe
scalar o0 = e(ll)
regsave using resultsnoint, addlabel(Outcome, Contacts, Model, Baseline, Specification, Estimated) append
xtreg outside i.t treated, fe
scalar oall = e(ll)
regsave using resultsnoint, addlabel(Outcome, Contacts, Model, Baseline, Specification, True) append
xtreg active i.t treated if mysample==1, fe
scalar ac0 = e(ll)
regsave using resultsnoint, addlabel(Outcome, Active, Model, Baseline, Specification, Estimated) append
xtreg active i.t treated, fe
scalar acall = e(ll)
regsave using resultsnoint, addlabel(Outcome, Active, Model, Baseline, Specification, True) append

use resultsnoint, clear
keep if var=="treated"

foreach spec in "True" "Estimated"{
  gen _`spec' = coef if Spec == "`spec'"
  gen _`spec'sd = stderr if Spec == "`spec'"
  gen _`spec'N = N if Spec=="`spec'"
  egen __`spec' = max(_`spec'), by(Model Outcome)
  egen __`spec'sd = max(_`spec'sd), by(Model Outcome)
  egen `spec'N = max(_`spec'N), by(Model Outcome)
  gen `spec' = string(__`spec',"%10.3f") + " (" + string(__`spec'sd,"%10.3f") + ")"
  gen `spec'_nosd = string(__`spec',"%10.3f")
}

sort Outcome Model
replace Outcome = "\multirow{2}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model True_nosd Estimated
listtex Outcome Model True_nosd Estimated, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & Estimated  \\ \midrule") foot("\bottomrule \end{tabular}")



use results, clear
keep if var=="1.treated"

foreach spec in "True" "Estimated"{
  gen _`spec' = coef if Spec == "`spec'"
  gen _`spec'sd = stderr if Spec == "`spec'"
  gen _`spec'N = N if Spec=="`spec'"
  egen __`spec' = max(_`spec'), by(Model Outcome)
  egen __`spec'sd = max(_`spec'sd), by(Model Outcome)
  egen `spec'N = max(_`spec'N), by(Model Outcome)
  gen `spec' = string(__`spec',"%10.3f") + " (" + string(__`spec'sd,"%10.3f") + ")"
  gen `spec'_nosd = string(__`spec',"%10.3f")
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model True_nosd Estimated
listtex Outcome Model True_nosd Estimated, type rstyle(tabular) head("\begin{tabular}{llcc}" "\midrule Outcome & Model & True & Estimated  \\ \midrule") foot("\bottomrule \end{tabular}")

use results, clear
keep if var=="1.treated#c.density"

foreach spec in "True" "Estimated"{
  gen _`spec' = coef if Spec == "`spec'"
  gen _`spec'sd = stderr if Spec == "`spec'"
  gen _`spec'N = N if Spec=="`spec'"
  egen __`spec' = max(_`spec'), by(Model Outcome)
  egen __`spec'sd = max(_`spec'sd), by(Model Outcome)
  egen `spec'N = max(_`spec'N), by(Model Outcome)
  gen `spec' = string(__`spec',"%10.3f") + " (" + string(__`spec'sd,"%10.3f") + ")"
  gen `spec'_nosd = string(__`spec',"%10.3f")
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model True_nosd Estimated
listtex Outcome Model True_nosd Estimated, type rstyle(tabular) head("\begin{tabular}{llcc}" "\midrule Outcome & Model & True &  Estimated \\ \midrule") foot("\bottomrule \end{tabular}")

log close
