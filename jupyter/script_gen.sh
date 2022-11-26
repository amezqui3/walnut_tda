#! /bin/sh

walnut_id=("2011SBb_R3_T9" "2014SBa_R5_T43" "2014SBa_R5_T90" "2014SBa_R7_T19" \
"201SBa_R5_T55" "2012SB_R11_T29" "2014SBa_R5_T47" "2014SBa_R5_T82" "2014SBa_R5_T91" \
"2014SBa_R7_T20" "2012SB_R12_T62" "2014SBa_R5_T48" "2014SBa_R5_T92" "2014SBa_R7_T21" \
"2012SB_R16_T66" "2014SBa_R5_T85" "2014SBa_R7_T23" "2014SBa_R5_T51" "2014SBa_R5_T86" \
"2014SBa_R6_T48" "2014SBa_R7_T24" "2014SBa_R5_T52" "2014SBa_R5_T87" "2014SBa_R7_T25" \
"2014SBa_R1_T33" "2014SBa_R5_T53" "2014SBa_R5_T89" "2014SBa_R6_T66" "2014SBa_R7_T26")

for wid in ${walnut_id[@]}
do
	echo "python3 Code/walnut/03_nut_alignment.py Results/walnut/clean/ Results/walnut/rotated/ "$wid
done

