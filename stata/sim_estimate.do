log using sim_estimate.log, text replace
insheet using ../output/nc5-policy-20-80-30pc.csv, comma clear
gen spatial = 1
save alldata, replace
insheet using ../output/nc5-SIR-policy-20-30pc.csv, comma clear
gen spatial = 0
replace naics = naics+1000
append using alldata

xtset naics t
summ
set more off
sort naics t


*keep if t<80

replace npi_date = . if npi_date == 999
gen growth = D.active/L.active
gen growth2 = ln(active/L.active)
gen delta = D.active/L.active

summ gr*


gen treat = t>=npi_date & npi_date!=.
gen X = susceptible*active/25600
regress D.active L.X L.active if treat==0 & spatial==0, noconstant

reg growth X if naics==1

reg growth treat##naics

xtreg growth treat##naics i.t, fe
bysort treat: summ growthi
reg growthi treat##naics

summ

// try with SIR model to reproduce unit slope missed in sim_estimate_dens
insheet using ../output/nc5-SIR-dens.csv, comma clear
quietly xtset naics t
summ

sort naics t
drop if npi_date==20
egen maxactive = max(active), by(naics)
// gen tmax = t if maxactive==active
// egen tm = max(tmax), by(naics)
// drop if t > tm

* difference in infection between t+1 and t
gen deltai = D.F1.active
tab density, gen(densdum)
gen growthi = deltai/active
gen X1 = susceptible/25600

reg growth c.X1##densdum*, cluster(naics)

*matrix list e(b)
matrix b = e(b)

*figure out how many betas there are and correct
unab vars: densdum*
local nbetas =  `: word count `vars''
di `nbetas'

local beta`nbetas' = b[1,1]
local lastbeta = `nbetas'-1
forvalues i = 1(1)`lastbeta' {
  local beta`i' = b[1,1+`nbetas'*2 +`i'*2] + b[1,1]
}
di b[1,33]
di `beta11'
di `beta3'
gen beta = .

forvalues i = 1(1)`nbetas' {
  qui replace beta = `beta`i'' if densdum`i'==1
}

preserve
collapse beta, by(density)
list

gen normbeta = beta/.15
list
restore

reg growth c.X1##c.density, cluster(naics)

gen dX = density*X1
reg growth X1 dX
restore
log close
