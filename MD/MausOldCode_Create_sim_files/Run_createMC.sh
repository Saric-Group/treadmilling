#!/usr/bin/tcsh

set steps=100000
set Lx=25.0
set Ly=25.0
set MC_steps=100
# foreach Tau ( 10000.0 100000.0 1000000.0 10000000.0 )

foreach PBB ( 1.0 )
foreach Tau ( 675.0 )
foreach mu ( 4.0 )
foreach Nsim ( 3 )


python3 Create_config_MC.py ${Tau} ${mu} ${PBB} ${Nsim} ${Lx} ${Ly} ${steps} ${MC_steps}

./lmp_serial -in MC_One_Filament/Tau_${Tau}_mu_${mu}_PBB_${PBB}_Lx_${Lx}_Ly_${Ly}/Inputs_MC/in.MC_OneFilament_Tau_${Tau}_mu_${mu}_PBB_${PBB}_Lx_${Lx}_Ly_${Ly}_Nsim_${Nsim}

end 
end 
end 
end 