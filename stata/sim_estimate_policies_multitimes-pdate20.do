//
capture log close
log using sim_estimate_policies_multitimes-pdate20.txt, text replace

*local filetype = "dens-beh_p"
local filetype = "dens"
*insheet using ../output/nc5-`filetype'-times-25pc.csv, comma clear
insheet using ../output/nc5-`filetype'-20-80-25pc.csv, comma clear
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

drop if density >.5 & density < .6

gen mysample = 0
replace mysample = 1 if npi_date==999 & density>=1
replace mysample = 1 if npi_date==20 & density<1
gen treated = 0
replace treated = 1 if npi_date==20 & t>=20


xtreg growthi i.t treated##c.density if mysample==1, fe
scalar g0 = e(ll)
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, Estimated) replace
xtreg growthi i.t treated##c.density, fe
scalar gall = e(ll)
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, True) append
xtreg outside i.t treated##c.density if mysample==1, fe
scalar o0 = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, Estimated) append
xtreg outside i.t treated##c.density, fe
scalar oall = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, True) append
xtreg active i.t treated##c.density if mysample==1, fe
scalar ac0 = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baseline, Specification, Estimated) append
xtreg active i.t treated##c.density, fe
scalar acall = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baseline, Specification, True) append


preserve
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
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model Est* True*
*texsave using table.tex, replace nofix preamble("\usepackage{multirow}") align(llcc) width(.8\linewidth)

listtex Outcome Model True TrueN Estimated EstimatedN, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & N & Estimated & N \\ \midrule") foot("\bottomrule \end{tabular}")

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
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model Est* True*
*texsave using table.tex, replace nofix preamble("\usepackage{multirow}") align(llcc) width(.8\linewidth)

listtex Outcome Model True TrueN Estimated EstimatedN, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & N & Estimated & N \\ \midrule") foot("\bottomrule \end{tabular}")

restore



* BEGIN attempts to do predictions TO REVIEW
preserve
tab density
gen temp = active if npi_date==999
egen actnotreat = mean(active), by(density t)
drop temp
xtreg active i.t treated##c.density, fe
predict predactive
gen realtreated = treated
gen notreated = npi_date==999
tab notreat npi_d
replace treated=1 if t>=20
predict myactive20
*keep if realtreated==1
collapse active myactive20 npi_date actnotreat , by(t realtreated)
keep if t>=20
twoway line active t if realtreated==1 || line myactive t if realtreated==1|| line actnotreat t if realtreated==1 || line active t if realtreated==0|| line myactive t if realtreated==0 & t<60
summ active myactive
gen diff = actnotreat - myactive20
summ diff if t>=20
sort t
keep if realtreated==1
keep active myactive t
save predictions, replace

restore


preserve
drop if npi_date==40

gen temp = active if npi_date==999
egen actnotreat = mean(active), by(density t)
drop temp
qui xtreg active i.t treated##c.density if mysample==1, fe
predict predactive
gen realtreated = treated
replace treated=1 if t>=20
predict myactive20
keep if realtreated==1
collapse active myactive20 npi_date actnotreat , by(t)
twoway line active t|| line myactive t  || line actnotreat t, legend(pos(1))
summ active myactive
keep active myactive t
rename active active_dbias
rename myactive myactive20_dbias
sort t
merge 1:1 t using predictions
drop _merge
save predictions, replace
restore

preserve
drop if npi_date==40
gen temp = outside if npi_date==999
egen actnotreat = mean(outside), by(density t)
drop temp
qui xtreg outside i.t treated##c.density if mysample==1, fe
predict predoutside
gen realtreated = treated
replace treated=1 if t>=20
predict myoutside20
keep if realtreated==1
collapse outside myoutside20 npi_date actnotreat , by(t)
twoway line outside t|| line myoutside t // || line actnotreat t, legend(pos(1))
summ outside myoutside
keep outside myoutside t
merge 1:1 t using predictions
drop _merge
save predictions, replace
restore


preserve
drop if npi_date==40
gen temp = outside if npi_date==999
egen actnotreat = mean(outside), by(density t)
drop temp
qui xtreg outside i.t treated##c.density if mysample==1, fe
predict predoutside
gen realtreated = treated
replace treated=1 if t>=20
predict myoutside20
keep if realtreated==1
collapse outside myoutside20 npi_date actnotreat , by(t)
twoway line outside t|| line myoutside t  || line actnotreat t, legend(pos(1))
summ outside myoutside
keep outside myoutside t
rename outside outside_dbias
rename myoutside myoutside20_dbias
merge 1:1 t using predictions
drop _merge
save predictions, replace

restore

preserve
drop if npi_date==40
drop if growthi==.
tab mysample density
qui xtreg growthi i.t treated##c.density, fe
predict predgrowthi
gen realtreated = treated
replace treated=1 if t>=20
predict mygrowthi20
keep if realtreated==1
collapse growthi mygrowthi20 npi_date, by(t)
twoway line growthi t|| line mygrowthi t , legend(pos(1))
summ growthi mygrowthi
keep growthi mygrowthi t
merge 1:1 t using predictions
drop _merge
save predictions, replace
restore

preserve
drop if npi_date==40
drop if growthi==.
qui xtreg growthi i.t treated##c.density if mysample==1, fe
predict predgrowthi
gen realtreated = treated
replace treated=1 if t>=20
predict mygrowthi20
keep if realtreated==1
collapse growthi mygrowthi20 npi_date, by(t)
twoway line growthi t|| line mygrowthi t
summ growthi mygrowthi
keep growthi mygrowthi t
rename growthi growthi_dbias
rename mygrowthi mygrowthi20_dbias
merge 1:1 t using predictions
drop _merge
save predictions, replace
restore


***********************
*Now with behavior

local filetype = "dens-beh_p"
insheet using ../output/nc5-`filetype'-20-80-25pc.csv, comma clear
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

drop if density >.5 & density < .6

gen mysample = 0
replace mysample = 1 if npi_date==999 & density>=1
replace mysample = 1 if npi_date==20 & density<1
gen treated = 0
replace treated = 1 if npi_date==20 & t>=20


xtreg growthi i.t treated##c.density if mysample==1, fe
scalar g0 = e(ll)
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, Estimated) replace
xtreg growthi i.t treated##c.density, fe
scalar gall = e(ll)
regsave using results, addlabel(Outcome, Growth rate, Model, Baseline, Specification, True) append
xtreg outside i.t treated##c.density if mysample==1, fe
scalar o0 = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, Estimated) append
xtreg outside i.t treated##c.density, fe
scalar oall = e(ll)
regsave using results, addlabel(Outcome, Contacts, Model, Baseline, Specification, True) append
xtreg active i.t treated##c.density if mysample==1, fe
scalar ac0 = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baseline, Specification, Estimated) append
xtreg active i.t treated##c.density, fe
scalar acall = e(ll)
regsave using results, addlabel(Outcome, Active, Model, Baseline, Specification, True) append


preserve
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
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model Est* True*
*texsave using table.tex, replace nofix preamble("\usepackage{multirow}") align(llcc) width(.8\linewidth)

listtex Outcome Model True TrueN Estimated EstimatedN, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & N & Estimated & N \\ \midrule") foot("\bottomrule \end{tabular}")

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
}

sort Outcome Model
replace Outcome = "\multirow{1}{*}{"+Outcome+"}"
replace Outcome = "" if Model=="With behavior"
drop if Spec == "Estimated"
drop _* r2 coef stderr var Spec N

list Outcome Model Est* True*
*texsave using table.tex, replace nofix preamble("\usepackage{multirow}") align(llcc) width(.8\linewidth)

listtex Outcome Model True TrueN Estimated EstimatedN, type rstyle(tabular) head("\begin{tabular}{llcccc}" "\midrule Outcome & Model & True & N & Estimated & N \\ \midrule") foot("\bottomrule \end{tabular}")

restore



* BEGIN attempts to do predictions TO REVIEW
preserve
tab outside
gen temp = active if npi_date==999
egen actnotreat = mean(active), by(density t)
drop temp
xtreg active i.t treated##c.density, fe
predict predactive
gen realtreated = treated
gen notreated = npi_date==999
tab notreat npi_d
replace treated=1 if t>=20
predict myactive20
*keep if realtreated==1
collapse active myactive20 npi_date actnotreat , by(t realtreated)
keep if t>=20
twoway line active t if realtreated==1 || line myactive t if realtreated==1|| line actnotreat t if realtreated==1 || line active t if realtreated==0|| line myactive t if realtreated==0 & t<60
summ active myactive
gen diff = actnotreat - myactive20
summ diff if t>=20
sort t
keep if realtreated==1
keep active myactive t
save predictions_beh, replace

restore

preserve
drop if npi_date==40

gen temp = active if npi_date==999
egen actnotreat = mean(active), by(density t)
drop temp
qui xtreg active i.t treated##c.density if mysample==1, fe
predict predactive
gen realtreated = treated
replace treated=1 if t>=20
predict myactive20
keep if realtreated==1
collapse active myactive20 npi_date actnotreat , by(t)
twoway line active t|| line myactive t  || line actnotreat t, legend(pos(1))
summ active myactive
keep active myactive t
rename active active_dbias
rename myactive myactive20_dbias
sort t
merge 1:1 t using predictions_beh
drop _merge
save predictions_beh, replace
restore

preserve
drop if npi_date==40
gen temp = outside if npi_date==999
egen outnotreat = mean(outside), by(density t)
drop temp
qui xtreg outside i.t treated##c.density, fe
predict predoutside
gen realtreated = treated
replace treated=1 if t>=20
predict myoutside20
keep if realtreated==1
collapse outside myoutside20 npi_date outnotreat , by(t)
twoway line outside t|| line myoutside t // || line actnotreat t, legend(pos(1))
summ outside myoutside
keep outside myoutside t
merge 1:1 t using predictions_beh
drop _merge
save predictions_beh, replace
restore

preserve
drop if npi_date==40
gen temp = outside if npi_date==999
egen actnotreat = mean(outside), by(density t)
drop temp
qui xtreg outside i.t treated##c.density if mysample==1, fe
predict predoutside
gen realtreated = treated
replace treated=1 if t>=20
predict myoutside20
keep if realtreated==1
collapse outside myoutside20 npi_date actnotreat , by(t)
twoway line outside t|| line myoutside t  || line actnotreat t, legend(pos(1))
summ outside myoutside
keep outside myoutside t
rename outside outside_dbias
rename myoutside myoutside20_dbias
merge 1:1 t using predictions_beh
drop _merge
save predictions_beh, replace

restore

preserve
drop if npi_date==40
drop if growthi==.
tab mysample density
qui xtreg growthi i.t treated##c.density, fe
predict predgrowthi
gen realtreated = treated
replace treated=1 if t>=20
predict mygrowthi20
keep if realtreated==1
collapse growthi mygrowthi20 npi_date, by(t)
twoway line growthi t|| line mygrowthi t , legend(pos(1))
summ growthi mygrowthi
keep growthi mygrowthi t
merge 1:1 t using predictions_beh
drop _merge
save predictions_beh, replace
restore

preserve
drop if npi_date==40
drop if growthi==.
qui xtreg growthi i.t treated##c.density if mysample==1, fe
predict predgrowthi
gen realtreated = treated
replace treated=1 if t>=20
predict mygrowthi20
keep if realtreated==1
collapse growthi mygrowthi20 npi_date, by(t)
twoway line growthi t|| line mygrowthi t
summ growthi mygrowthi
keep growthi mygrowthi t
rename growthi growthi_dbias
rename mygrowthi mygrowthi20_dbias
merge 1:1 t using predictions_beh
drop _merge
save predictions_beh, replace
restore

use predictions_beh, clear
outsheet using ../output/predictions_beh.csv, comma replace

use predictions, clear
outsheet using ../output/predictions.csv, comma replace

log close
