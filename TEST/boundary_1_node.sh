#!/bin/bash

ScriptDir=`pwd`

# QHA_steping_stone work:
cp QHA_steping_stone.py  ./Files_Outputs/Calcite_I
cp QHA_steping_stone.py  ./Files_Outputs/Calcite_II

cd ./Files_Outputs/Calcite_I
echo "
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You are now perforing the QHA_steping_stone work for Calcite I"
python QHA_steping_stone.py
cd $ScriptDir

cd ./Files_Outputs/Calcite_II
echo "
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You are now perforing the QHA_steping_stone work for Calcite II"
python QHA_steping_stone.py
cd $ScriptDir

rm -Rf EL_level
mkdir EL_level

rm -Rf EL_plus_E0_level
mkdir EL_plus_E0_level

rm -Rf G_PT
mkdir G_PT

rm -Rf All_P_H_real_intersection_at_Ts
mkdir All_P_H_real_intersection_at_Ts

rm -Rf ET
mkdir ET

rm -Rf TS
mkdir TS

rm -Rf H_PT
mkdir H_PT


mkdir -p EL_level
cp  InTerSect_EL_level.py ./EL_level

mkdir -p ./EL_plus_E0_level
cp InTerSect_EL_plus_E0_level.py ./EL_plus_E0_level


cd ./EL_level
python InTerSect_EL_level.py
cd $ScriptDir

cd ./EL_plus_E0_level
python InTerSect_EL_plus_E0_level.py
cd $ScriptDir


FOLDERS="
10.00K
30.10K
50.20K
70.30K
90.40K
110.51K
130.61K
150.71K
170.81K
190.91K
211.01K
231.11K
251.21K
271.31K
291.41K
311.52K
331.62K
351.72K
371.82K
391.92K
412.02K
432.12K
452.22K
472.32K
492.42K
512.53K
532.63K
552.73K
572.83K
592.93K
"
#613.03K
#633.13K
#653.23K
#673.33K
#693.43K
#713.54K
#733.64K
#753.74K
#773.84K
#793.94K
#814.04K
#834.14K
#854.24K
#874.34K
#894.44K
#914.55K
#934.65K
#954.75K
#974.85K
#994.95K
#1015.05K
#1035.15K
#1055.25K
#1075.35K
#1095.45K
#1115.56K
#1135.66K
#1155.76K
#1175.86K
#1195.96K
#1216.06K
#1236.16K
#1256.26K
#1276.36K
#1296.46K
#1316.57K
#1336.67K
#1356.77K
#1376.87K
#1396.97K
#1417.07K
#1437.17K
#1457.27K
#1477.37K
#1497.47K
#1517.58K
#1537.68K
#1557.78K
#1577.88K
#1597.98K
#1618.08K
#1638.18K
#1658.28K
#1678.38K
#1698.48K
#1718.59K
#1738.69K
#1758.79K
#1778.89K
#1798.99K
#1819.09K
#1839.19K
#1859.29K
#1879.39K
#1899.49K
#1919.60K
#1939.70K
#1959.80K
#1979.90K
#2000.00K
#"

for i in ${FOLDERS}; do

mkdir -p G_PT/${i}
cp  InTerSect_G_PT_level.py ./G_PT/${i}

cd ./G_PT/${i}
INPUT="filefolder_energetics = 'EL_vs_V'"
OUTPUT="filefolder_energetics = 'F_vs_V_${i}'"

INPUT_2="filefolder_energetics, 'EL_vs_V.dat'"
OUTPUT_2="filefolder_energetics, 'F_vs_V_${i}.dat'"

echo $INPUT
echo $OUTPUT

echo $INPUT_2
echo $OUTPUT_2

sed -i "s/${INPUT}/${OUTPUT}/g" InTerSect_G_PT_level.py 
sed -i "s/${INPUT_2}/${OUTPUT_2}/g" InTerSect_G_PT_level.py 
python InTerSect_G_PT_level.py 

cp P_H_analytic_intersection* $ScriptDir/All_P_H_real_intersection_at_Ts

cd $ScriptDir

done

########### Analytical cat-ing:
cd $ScriptDir/All_P_H_real_intersection_at_Ts
cat P_H_analytic_intersection* > All_TEMPERATS_P_H_analytic_intersection.dat

#cat $ScriptDir/EL_level/P_H_analytic_intersection*  >> All_TEMPERATS_P_H_analytic_intersection.dat
cat $ScriptDir/EL_plus_E0_level/P_H_analytic_intersection* >>  All_TEMPERATS_P_H_analytic_intersection.dat
grep -v "Temperature" All_TEMPERATS_P_H_analytic_intersection.dat > templat && mv templat  All_TEMPERATS_P_H_analytic_intersection.dat
 
sort -k1 -n All_TEMPERATS_P_H_analytic_intersection.dat > templat2 && mv templat2 All_TEMPERATS_P_H_analytic_intersection.dat

(echo "# Temperature of Intersection(K)        Pressure of Intersection (GPa)" && cat All_TEMPERATS_P_H_analytic_intersection.dat) > templat3 && mv templat3  All_TEMPERATS_P_H_analytic_intersection.dat

cd $ScriptDir


for i in ${FOLDERS}; do

mkdir -p ET/${i}
cp  InTerSect_ET_level.py ./ET/${i}

cd ./ET/${i}
INPUT="filefolder_energetics = 'EL_vs_V'"
OUTPUT="filefolder_energetics = 'ET_vs_V_${i}'"

INPUT_2="filefolder_energetics, 'EL_vs_V.dat'"
OUTPUT_2="filefolder_energetics, 'ET_vs_V_${i}.dat'"

echo $INPUT
echo $OUTPUT

echo $INPUT_2
echo $OUTPUT_2

sed -i "s/${INPUT}/${OUTPUT}/g" InTerSect_ET_level.py
sed -i "s/${INPUT_2}/${OUTPUT_2}/g" InTerSect_ET_level.py
python InTerSect_ET_level.py

cd $ScriptDir

done

for i in ${FOLDERS}; do

mkdir -p TS/${i}
cp  InTerSect_TS_level.py ./TS/${i}

cd ./TS/${i}
INPUT="filefolder_energetics = 'EL_vs_V'"
OUTPUT="filefolder_energetics = 'TS_vs_V_${i}'"

INPUT_2="filefolder_energetics, 'EL_vs_V.dat'"
OUTPUT_2="filefolder_energetics, 'TS_vs_V_${i}.dat'"

echo $INPUT
echo $OUTPUT

echo $INPUT_2
echo $OUTPUT_2

sed -i "s/${INPUT}/${OUTPUT}/g" InTerSect_TS_level.py
sed -i "s/${INPUT_2}/${OUTPUT_2}/g" InTerSect_TS_level.py
python InTerSect_TS_level.py

cd $ScriptDir

done


for i in ${FOLDERS}; do

mkdir -p H_PT/${i}
cp  InTerSect_H_PT_level.py ./H_PT/${i}

cd ./H_PT/${i}
INPUT="filefolder_energetics = 'EL_vs_V'"
OUTPUT="filefolder_energetics = 'F_vs_V_${i}'"

INPUT_2="filefolder_energetics, 'EL_vs_V.dat'"
OUTPUT_2="filefolder_energetics, 'F_vs_V_${i}.dat'"

INPUT_3="'ET_vs_V', 'ET_vs_V_a', 'ET_vs_V_a.dat'"
OUTPUT_3="'ET_vs_V', 'ET_vs_V_${i}', 'ET_vs_V_${i}.dat'"


echo $INPUT
echo $OUTPUT

echo $INPUT_2
echo $OUTPUT_2

sed -i "s/${INPUT}/${OUTPUT}/g" InTerSect_H_PT_level.py
sed -i "s/${INPUT_2}/${OUTPUT_2}/g" InTerSect_H_PT_level.py
sed -i "s/${INPUT_3}/${OUTPUT_3}/g" InTerSect_H_PT_level.py
python InTerSect_H_PT_level.py

cd $ScriptDir

done

python Phase_Boundary.py

#rm -Rf PDF_summary
#mkdir PDF_summary
#mv make_pdflatex.sh main.tex ./PDF_summary
#cd ./PDF_summary
#bash make_pdflatex.sh


