#! /bin/sh

walnut_ida=("2011SBb_R3_T9" "2014SBa_R5_T43" "2014SBa_R5_T90" "2014SBa_R7_T19" \
"201SBa_R5_T55" "2012SB_R11_T29" "2014SBa_R5_T47" "2014SBa_R5_T82" "2014SBa_R5_T91" \
"2014SBa_R7_T20" "2012SB_R12_T62" "2014SBa_R5_T48" "2014SBa_R5_T92" "2014SBa_R7_T21" \
"2012SB_R16_T66" "2014SBa_R5_T85" "2014SBa_R7_T23" "2014SBa_R5_T51" "2014SBa_R5_T86" \
"2014SBa_R6_T48" "2014SBa_R7_T24" "2014SBa_R5_T52" "2014SBa_R5_T87" "2014SBa_R7_T25" \
"2014SBa_R1_T33" "2014SBa_R5_T53" "2014SBa_R5_T89" "2014SBa_R6_T66" "2014SBa_R7_T26")

walnut_idb=("2014SBa_R6_T35" "2014SBa_R5_T9" "2014SBa_R1_T35" "2014SBa_R1_T7" \
"2012SB_R16_T64" "2014SBa_R6_T12" "2014SBa_R5_T19" "2014SBa_R1_T26" "2014SBa_R1_T23" \
"2012SB_R16_T58" "2014SBa_R7_T28" "2014SBa_R5_T84" "2014SBa_R5_T21" "2014SBa_R1_T27" \
"2012SB_R11_T69" "2014SBa_R6_T49" "2014SBa_R1_T6" "2012SB_R12_T57" "2014SBa_R6_T13" \
"2014SBa_R1_T20" "2014SBa_R1_T18" "2014SBa_R7_T18")

walnut_idc=("SelectF_R34_T27" "NewStuke_R7_T13" "NewStuke_R5_T9" "NewStuke_R3_T8" "GB_R11_T2" "2013SB_R22_T49" "2012SB_R16_T16" "2008SB_R5_T2" "SelectF_R34_T36" "SelectF_R30_T35" "NewStuke_R5_T7" "GB_R5_T9" "NewStuke_R7_T7" "NewStuke_R1_T6" "2014SBa_R6_T46" "2010SB_R1_T27" "2008SB_R8_T10" "SelectD_R1_T3" "NewStuke_R9_T15" "NewStuke_R9_T11" "GB_R9_T3" "GB_R1_T1" "2014SBa_R4_T19" "2013SB_R21_T42" "2012SB_R15_T56" "2012SB_R14_T63" "2010SB_R4_T39" "SelectF_R32_T21" "SelectF_R25_T25" "GB_R9_T2" "GB_R5_T12" "2012SB_R14_T34" "2011SBa_R22_T20" "SelectF_R34_T26" "SelectF_R32_T13" "SelectC_R1_T11" "NewStuke_R9_T18" "GB_R11_T5" "2012SB_R8_T68" "2010SB_R10_T26" "2010SB_R2_T24" "SelectF_R30_T13" "SelectF_R28_T25" "SelectD_R7_T1" "NewStuke_R3_T17" "GB_R7_T3" "2014SBa_R1_T36" "2014SBa_R1_T29" "2013SB_R23_T45" "2012SB_R17_T41" "GB_R3_T2" "2014SBa_R1_T30" "2013SB_R23_T31" "2011SBb_R4_T20" "2008SB_R4_T13" "2014SBa_R7_T17" "2014SBa_R5_T13" "2014SBa_R1_T32" "2012SB_R16_T67" "SelectF_R32_T17" "NewStuke_R9_T6" "NewStuke_R7_T9" "NewStuke_R1_T18" "NewStuke_R1_T17" "GB_R11_T3" "2014SBa_R6_T40" "2014SBa_R5_T20" "2014SBa_R1_T10" "NewStuke_R3_T14" "2014SBa_R4_T17" "2014SBa_R1_T21" "2012SB_R11_T22" "2012SB_R6_T6" "2011SBb_R4_T36" "2011SBb_R3_T68" "GB_R13_T3" "GB_R13_T2" "GB_R5_T3" "2014SBa_R6_T34" "2012SB_R15_T66" "2014SBa_R6_T36" "2011SBb_R3_T20" "GB_R5_T4" "2011SBb_R4_T26" "2011SBb_R3_T58" "2014SBa_R6_T45" "2014SBa_R6_T17" "2014SBa_R5_T17" "2014SBa_R5_T14" "2011SBb_R4_T65" "2014SBa_R1_T28" "2012SB_R11_T43" "2014SBa_R6_T65" "2014SBa_R6_T55" "2014SBa_R1_T15" "2014SBa_R6_T54" "2014SBa_R4_T16" "2014SBa_R1_T11" "2011SBb_R3_T60" "2014SBa_R6_T53" "2014SBa_R6_T14" "2014SBa_R5_T54" "2014SBa_R5_T16" "2014SBa_R1_T57" "2014SBa_R1_T5" "2011SBb_R4_T44" "2011SBb_R3_T56" "2014SBa_R6_T43" "2014SBa_R5_T15" "2014SBa_R1_T17" "2012SB_R12_T60" "2011SBb_R3_T65" "2014SBa_R6_T47")

walnut_idd=("SelectD_R9_T1" "2014SBa_R6_T64" "2014SBa_R6_T16" "2014SBa_R5_T83" "2014SBa_R5_T81" "2014SBa_R5_T49" "2014SBa_R1_T3" "2014SBa_R1_T25" "2014SBa_R16_T58")

walnuts=("${walnut_ida[@]}" "${walnut_idb[@]}" "${walnut_idc[@]}" "${walnut_idd[@]}")

for wid in ${walnuts[@]}
do
	#echo "python3 Code/walnut/01_density_normalization.py Stacks/Walnuts/ Results/walnut/clean/ "$wid
	#echo "python3 Code/walnut/02_watershed_segmentation.py Results/walnut/clean/ Results/walnut/watershed/ "$wid
	#echo "python3 Code/walnut/03_nut_alignment.py Results/walnut/clean/ Results/walnut/rotated/ "$wid
	#echo "python3 Code/walnut/04_interior_shell.py Results/walnut/ Results/walnut/watershed/ "$wid
#	echo "python3 Code/walnut/05_traditional_phenotyping.py Results/walnut/ Results/walnut/traditional/ "$wid
	#echo "python3 Code/walnut/06_ect_on_walnut.py Results/walnut/ Results/walnut/topology/ "$wid
	echo "python3 Code/walnut/07_ect_on_kernel.py Results/walnut/ Results/walnut/topology/ "$wid
	
done

