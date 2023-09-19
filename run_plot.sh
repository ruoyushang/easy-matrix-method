
sh clean.sh

#epoch='V5'
#epoch='V6'
#epoch='V5V6'
epoch='V4V5V6'
#onoff='OFF'
onoff='ON'

#energy_start=0
#energy_middle=4
#energy_end=8
#energy_start=1
#energy_middle=2
#energy_end=8

#energy_start=5
#energy_middle=8
#energy_end=11

#energy_start=0
#energy_middle=6
#energy_end=12

energy_start=2
energy_middle=6
energy_end=12
#energy_start=6
#energy_middle=9
#energy_end=12

#src='CrabNebula_elev_80_90'
#src='CrabNebula_elev_70_80'
#src='CrabNebula_elev_60_70'
#src='CrabNebula_elev_50_60'
#src='CrabNebula_elev_40_50'

#src='CrabNebula_1p0wobble'
#src='CrabNebula_1p5wobble'
#
#src='UrsaMajorII'  
#src='UrsaMinor' 
#src='RGB_J0710_p591' 
#src='1ES0229' 
#src='PKS1424' 
#src='PG1553' 
#src='3C273' 
#src='Segue1' 
#src='NGC1275' 
#src='H1426' 
#src='OJ287' 
#src='Draco' 
#src='BLLac' 
#src='3C264' 
#src='1ES0502' 
#src='M82'
#src='1ES0414'
#src='1ES1011'
#src='1ES0647'

#src='PSR_J1907_p0602'
#src='PSR_J1856_p0245'
#src='SS433'

#src='SNR_G189_p03' # IC 443
#src='PSR_J2021_p4026' # Gamma Cygni
#src='Geminga'
#src='PSR_J2021_p3651' # Dragonfly
src='PSR_J2032_p4127'
#src='SS433_0p5deg'
#src='SS433_2p0deg'
#src='Cas_A'

#src='InputList_PSR_J2032_p4127_baseline'
#src='InputList_PSR_J2032_p4127_binary'

#src='PSR_J0631_p1036'
#src='CTA1'
#src='Tycho'
#src='SNR_G150_p4'
#src='2HWC_J1953_p294'
#src='Cisne_HS_2013'

#src='SNR_G067p6_p0p9'
#src='CTB109'
#src='LHAASO_J1956_p2845'
#src='PSR_J0359_p5414'
#src='LHAASO_J0622_p3754'

#src='InputList_Mrk421'

#src='SgrA'

python3 ImposterAnalysis.py $src $epoch $onoff $energy_start $energy_middle $energy_end
