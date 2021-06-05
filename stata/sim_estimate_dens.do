* generates data comparing SIR and Spatial-SIR estimate of beta
* with data generated fromm Spatial-SIR with different densities

local basedir = "output/nc5-"
local filetype = "dens"
capture log close
log using `basedir'sim_estimate_dens.log, text replace

*local filetype = "pop_cl"
*local filetype = "pop"

local filetype = "dens"
insheet using ../output/nc5-`filetype'-20-80-25pc.csv, comma clear
quietly xtset naics t
summ

sort naics t
*drop if npi_date==20
egen maxactive = max(active), by(naics)

* difference in infection between t+1 and t
gen deltai = D.F1.active
tab density, gen(densdum)
gen growthi = deltai/active
gen X1 = susceptible/q_popsize

reg growthi c.X1##densdum*, cluster(naics)
*reg growthi2 c.X1##densdum*, cluster(naics)

*matrix list e(b)
matrix b = e(b)
matrix se = e(V)
*figure out how many betas there are and correct

set more off
di sqrt(se[3,3])
di sqrt(se[5,5])

unab vars: densdum*
local nbetas =  `: word count `vars''
di `nbetas'

local beta`nbetas' = b[1,1]
local se`nbetas' = sqrt(se[1,1])
local lastbeta = `nbetas'-1
forvalues i = 1(1)`lastbeta' {
  local beta`i' = b[1,1+`nbetas'*2 +`i'*2] + `beta15'
  local se`i' = sqrt(se[1+`nbetas'*2 +`i'*2,1+`nbetas'*2 +`i'*2])
}
di `beta15'
di `beta1'
di `se15'
di `se1'
gen beta = .
gen std = .

forvalues i = 1(1)`nbetas' {
  qui replace beta = `beta`i'' if densdum`i'==1
  qui replace std = `se`i'' if densdum`i'==1
}
list beta std density in 1/5
preserve
collapse beta std, by(density)
list

outsheet density std beta using ../output/nc5-`filetype'betas.csv, comma replace
restore

// // now do only for observations before the peak
//
// drop beta std.
// gen tmax = t if maxactive==active
// egen tm = max(tmax), by(naics)
// drop if t > tm
//
// reg growthi c.X1##densdum*, cluster(naics)
// *reg growthi2 c.X1##densdum*, cluster(naics)
//
// *matrix list e(b)
// matrix b = e(b)
// matrix se = e(V)
// *figure out how many betas there are and correct
//
// set more off
// di sqrt(se[3,3])
// di sqrt(se[5,5])
//
// unab vars: densdum*
// local nbetas =  `: word count `vars''
// di `nbetas'
//
// local beta`nbetas' = b[1,1]
// local se`nbetas' = sqrt(se[1,1])
// local lastbeta = `nbetas'-1
// forvalues i = 1(1)`lastbeta' {
//   local beta`i' = b[1,1+`nbetas'*2 +`i'*2] + `beta15'
//   local se`i' = sqrt(se[1+`nbetas'*2 +`i'*2,1+`nbetas'*2 +`i'*2])
// }
// gen beta = .
// gen std = .
//
// forvalues i = 1(1)`nbetas' {
//   qui replace beta = `beta`i'' if densdum`i'==1
//   qui replace std = `se`i'' if densdum`i'==1
// }
// preserve
// collapse beta std, by(density)
// list
//
// outsheet density std beta using ../output/nc5-`filetype'betas-beforepeak.csv, comma replace
// restore

*** Now with behavioral

local filetype = "dens-beh_p"
insheet using ../output/nc5-`filetype'-20-80-25pc.csv, comma clear
quietly xtset naics t
summ

sort naics t
*drop if npi_date==20
egen maxactive = max(active), by(naics)

* difference in infection between t+1 and t
gen deltai = D.F1.active
tab density, gen(densdum)
gen growthi = deltai/active
gen X1 = susceptible/q_popsize

reg growthi c.X1##densdum*, cluster(naics)
*reg growthi2 c.X1##densdum*, cluster(naics)

*matrix list e(b)
matrix b = e(b)
matrix se = e(V)
*figure out how many betas there are and correct

set more off
di sqrt(se[3,3])
di sqrt(se[5,5])

unab vars: densdum*
local nbetas =  `: word count `vars''
di `nbetas'

local beta`nbetas' = b[1,1]
local se`nbetas' = sqrt(se[1,1])
local lastbeta = `nbetas'-1
forvalues i = 1(1)`lastbeta' {
  local beta`i' = b[1,1+`nbetas'*2 +`i'*2] + `beta15'
  local se`i' = sqrt(se[1+`nbetas'*2 +`i'*2,1+`nbetas'*2 +`i'*2])
}
di `beta15'
di `beta1'
di `se15'
di `se1'
gen beta = .
gen std = .

forvalues i = 1(1)`nbetas' {
  qui replace beta = `beta`i'' if densdum`i'==1
  qui replace std = `se`i'' if densdum`i'==1
}
list beta std density in 1/5
preserve
collapse beta std, by(density)
list

outsheet density std beta using ../output/nc5-`filetype'betas.csv, comma replace
restore
log close
