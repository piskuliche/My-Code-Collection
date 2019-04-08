topo readlammpsdata data.acncx
mol addfile traj.xyz type lammpstrj first 0 last -1 step 5 waitfor -1
pbc wrap -compound res -all
pbc box
mol selection type 1 2 3
mol material EdgyShiny
mol rep cpk 1.0 0.5 40 40
mol addrep 0
mol selection type 4 5 
mol material EdgyShiny
mol rep cpk 1.0 0.5 40 40
mol addrep 0
mol selection type 6 7 8
mol material EdgyShiny
mol rep vdw 2.0 40
mol addrep 0

color Display Background white
display depthcue off
display projection orthographic

