python sim-paper.py
/Applications/Stata/StataSE.app/Contents/MacOS/StataSE -b do stata/sim_estimate_dens.do &
/Applications/Stata/StataSE.app/Contents/MacOS/StataSE -b do stata/sim_estimate_policies_multitimes-pdate20.do &
/Applications/Stata/StataSE.app/Contents/MacOS/StataSE -b do stata/sim_estimate_policies_multitimes.do &
python fig-paper.py
