** this file imports data about Lombardy's infection from JHU's github covid-19 data repository

clear all

import delimited https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv
gen day = date(data,"YMD#")
format %td day

save "rawdata", replace
outsheet using rawdata.csv, comma replace


use rawdata
replace codice_regione=21 if denominazione_regione == "P.A. Bolzano"

xtset codice_regione day

rename codice_regione region
rename deceduti deaths
rename totale_positivi infect

*xtline deceduti

bys region: gen dgrowth = D.deaths/L.deaths
bys region: gen cgrowth = D.infect/L.infect

generate moveave3 = (F1.dgrowth + dgrowth + L1.dgrowth) / 3
generate moveave5 = (F1.dgrowth + F2.dgrowth + dgrowth + L1.dgrowth + L2.dgrowth) / 5
generate moveave7 = (F3.dgrowth + F1.dgrowth + F2.dgrowth + dgrowth + L1.dgrowth + L2.dgrowth + L3.dgrowth) / 7

outsheet moveave3 if region==3 & moveave3!=. using ../input/drlombardia.txt, replace noname
outsheet moveave5 if region==3 & moveave5!=. using ../input/drlombardia5.txt, replace noname
