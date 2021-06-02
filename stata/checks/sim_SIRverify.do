insheet using ../output/nc5-SIR-policy-20-30pc.csv, comma clear

xtset naics t
summ
set more off
sort naics t

replace npi_date = . if npi_date == 999
gen treat = t>=npi_date & npi_date!=.
gen post = t>20
gen X = susceptible*active/25600
gen deltai = D.F1.active

regress deltai X active



, noconstant

gen logdelta = log(F1.active)-log(active)
gen growthi = deltai/active
gen X1 = susceptible/25600

* basic regression
reg growthi X1

* verify bias of using logs
regress logdelta logS

* basic regression with treatment controls

reg growthi X1 treat##post
