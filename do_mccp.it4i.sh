#!/bin/bash
# sh do_mccp.sh workdir train.smi test.smi
##args
# set working directory wd
# set training file tr, e.g. Trainset_ratio10.txt
# set testing file te, e.g. test_1.txt
# set output directory od, e.g. testset_1
source /etc/bashrc
#source /home/excape/.bashrc
export PATH=~/libsvm/sklearn-openmp-master/scikit-learn/build/lib.linux-x86_64-3.5/sklearn:~/anaconda3/bin:~/bin:$PATH
export PYTHONPATH=~/anaconda3/lib/python3.5/site-packages:PYTHON$PATH

module load unsupported
module load java

if [ "$#" -lt 3 ]; then
   echo -e "At least 3 arguments, e.g.\nsh do_mccp.sh workingdirectory train.smi test.smi\n"
   exit
fi
wd=$1
#trf=$2
#tef=$3
if [ "$#" -eq 3 ]; then
  filename=$(basename "$3")
  odir="${filename%.*}"
  #extension=extension=$([[ "$filename" = *.* ]] && echo ".${filename##*.}" || echo '')
fi

if [ "$#" -ge 4 ]; then
odir=$4
fi

cd $wd
mkdir -p $odir

# convert to tab spererate file
awk 'BEGIN{FS="[ ,\t]";OFS="\t"} {$1=$1;print $0;}' $3 > $odir/$3.tsv
awk 'BEGIN{FS="[ ,\t]";OFS="\t"} {$1=$1;print $0;}' $2 > $odir/$2.tsv

#columns of training and test files
col1=$(head -1 $odir/$2.tsv | awk '{print NF}')
col2=$(head -1 $odir/$3.tsv | awk '{print NF}') 

echo -e "convert smiles to sdf and compute ECFP fingerprints in 2048bit...\n"
#training file should have a head, e.g. Smiles,id and Label (A or 1: active, N or 0: inactive)
tail -n +2 $odir/$2.tsv > $odir/$2.smi
if [ "$col1" -gt 2 ]
then
  sdfilter -i t=smiles -o t=sdf -d table=$odir/$2.tsv $odir/$2.smi $odir/$2.sdf
  mv sdfilter.trc $odir/$2.sdfilter.trc
  trl=$(head -1 $odir/$2.tsv | cut -f 2)
  java -Xmx6g -jar ~/jCMapperCLI.jar -hs 2048 -d 6 -a DAYLIGHT_INVARIANT_RING -c ECFP -l $trl -k true -ff LIBSVM_SPARSE -f $odir/$2.sdf -o $odir/$2.ECFP.2048bit 
  
  l1=$(wc -l < $odir/$2.ECFP.2048bit)
  l2=$(wc -l < $odir/$2.smi)
  if [ "$l1" -eq "$l2" ]; then
     cut -f 3 $odir/$2.smi > $odir/t1
     paste -d" " $odir/t1 $odir/$2.ECFP.2048bit | cut -d" " -f 1,3- > $odir/$2.ECFP.2048bit.libsvm.binary
     rm $odir/t1
  else
     awk 'BEGIN{FS="[ \t]";OFS=" "} NR==FNR{a[$2]=$3;next} {if($1 in a) print a[$1],$0;}' $odir/$2.tsv $odir/$2.ECFP.2048bit | cut -d" " -f 1,3- > $odir/$2.ECFP.2048bit.libsvm.binary
  fi
else
  echo -e "Training file:\nAt least 3 columns should be provided, first column is compound smiles, second is compound id like AZ123456,third is label of active (A or 1) or inacitve (N or 0)\n"
  exit
fi
rm $odir/$2.smi $odir/$2.tsv $odir/$2.sdf

#test file should have a head, e.g. Smiles, id, label can be empty
tail -n +2 $odir/$3.tsv > $odir/$3.smi
if [ "$col2" -gt 1 ]
then
  sdfilter -i t=smiles -o t=sdf -d table=$odir/$3.tsv $odir/$3.smi $odir/$3.sdf
  mv sdfilter.trc $odir/$3.sdfilter.trc
  tel=$(head -1 $odir/$3.tsv | cut -f 2)
  java -Xmx6g -jar ~/jCMapperCLI.jar -hs 2048 -d 6 -a DAYLIGHT_INVARIANT_RING -c ECFP -l $tel -k true -ff LIBSVM_SPARSE -f $odir/$3.sdf -o $odir/$3.ECFP.2048bit
  if [ "$col2" -gt 2 ] 
  then
     l1=$(wc -l < $odir/$3.ECFP.2048bit)
     l2=$(wc -l < $odir/$3.smi)
     if [ "$l1" -eq "$l2" ]; then
        cut -f 3 $odir/$3.smi > $odir/t2
        paste -d" " $odir/t2 $odir/$3.ECFP.2048bit | cut -d" " -f 1,3- > $odir/$3.ECFP.2048bit.libsvm.binary
        rm $odir/t2
     else
        awk 'BEGIN{FS="[ \t]";OFS=" "} NR==FNR{a[$2]=$3;next} {if($1 in a) print a[$1],$0;}' $odir/$3.tsv $odir/$3.ECFP.2048bit | cut -d" " -f 1,3- > $odir/$3.ECFP.2048bit.libsvm.binary
     fi
  else
    cut -d" " -f 2- $odir/$3.ECFP.2048bit | awk 'BEGIN{FS=" ";OFS=" "} {print 0,$0;}'  > $odir/$3.ECFP.2048bit.libsvm.binary
  fi
else
  echo -e "Test file:\nAt least 2 columns should be provided, first column is compound smiles, second is compound id like AZ123456!\n"
  exit
fi
rm $odir/$3.smi $odir/$3.tsv $odir/$3.sdf

# convert label A as 1 and N as 0
sed -i 's|A|1|' $odir/$2.ECFP.2048bit.libsvm.binary
sed -i 's|N|0|' $odir/$2.ECFP.2048bit.libsvm.binary

sed -i 's|A|1|' $odir/$3.ECFP.2048bit.libsvm.binary
sed -i 's|N|0|' $odir/$3.ECFP.2048bit.libsvm.binary


## or submit to slurm job
echo -e "\nRun MCCP...\n"
time OMP_NUM_THREADS=12 python ~/mccp_svm_openmp_sklearn.it4i.py -i $odir/$2.ECFP.2048bit.libsvm.binary -t $odir/$3.ECFP.2048bit.libsvm.binary -o $odir > $odir/$3.test.log

#post processing
echo -e "\nPost processing...\n"
if [ -s "$odir/mccp_svm_pred_rst.txt" ]
then
   echo -e "ID" > $odir/test_label
   cut -d" " -f 1 $odir/$3.ECFP.2048bit >> $odir/test_label
   paste $odir/test_label $odir/mccp_svm_pred_rst.txt > $odir/mccp_svm_pred_rst.ID.txt
   mv $odir/mccp_svm_pred_rst.ID.txt  $odir/mccp_svm_pred_rst.txt
   rm $odir/test_label $odir/$3.ECFP.2048bit.libsvm.binary $odir/$2.ECFP.2048bit.libsvm.binary $odir/$2.ECFP.2048bit
fi

echo -e "Done!\n"
