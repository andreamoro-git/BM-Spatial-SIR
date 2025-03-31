stataexecutable="/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b"
python sim-paper.py
$stataexecutable do stata/sim_estimate_dens.do &
read -t 5 'running stata sim_estimate_dens.do'
$stataexecutable do stata/sim_estimate_policies_multitimes-pdate20.do &
read -t 5 'running stata sim_estimate_policies_multitimes_pdate20.do'
$stataexecutable do stata/sim_estimate_policies_multitimes.do &
read -t 5 'running stata sim_estimate_policies_multitimes.do'
python fig-paper.py
