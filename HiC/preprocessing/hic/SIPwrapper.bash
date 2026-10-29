#!/bin/bash
version=1.0

hicpath=0
myres=0
help=0
norm=`echo "KR"`
#mychr=0
#juicepath=0
outdir=`echo "SIPfiles"`
strawpath=`echo "straw"`
sizefile=0
matsize=2000
cpu=1
sippath=0
d=6
g=1.5
max=2
min=2
sat=0.01
t=2800
nbzero=6
del=`echo "true"`
fdr=0.05
isDroso=`echo "false"`

menu=`printf "\nThis tool is a wrapper to use straw to dump and prepare files for SIP.\noptions: \n-hic\t\tinput .hic file\n-res\t\tresolution desired. Default is 5000\n-sizes\t\tchromosome size file with two columns corresponding to chromosome and size respectively\n-outdir\t\toutput directory names. Default is SIPfiles.\n-cpu\t\tnumber of threads to use; default is 1\n\n-straw\t\tSpecific path to straw. Default is to use straw in path.\n\n-norm\t\tnormalization scheme to dump\n-matsize\t\tSize of matrix to dump as in SIP manual. Default is 2000 bins.\n-sippath\t\tPath to SIP. Not necessary if you want to run separately.\n\nOptional Parameters only Passed to SIP. Defaults are same as SIP:\n-d\n-g\n-max\n-min\n-sat\n-t\n-nbzero\n-del\n-fdr\n-isDroso"`
while test $# -gt 0
do
        case "$1" in
                -d)
                d=$2
                ;;
                -g)
                g=$2
                ;;
                -max)
                max=$2
                ;;
                -min)
                min=$2
                ;;
                -sat)
                sat=$2
                ;;
                -t)
                t=$2
                ;;
                -nbzero)
                nbzero=$2
                ;;
                -del)
                del=$2
                ;;
                -fdr)
                fdr=$2
                ;;
                -isDroso)
                isDroso=$2
                ;;
                -res)
                myres=$2
                ;;
                -hic)
                hicpath=$2
                ;;
                -straw)
                strawpath=$2
                ;;
                -sizes)
                sizefile=$2
		;;
		-matsize)
		matsize=$2
                ;;
                -outdir)
                outdir=$2
                ;;
		-norm)
                norm=$2
                ;;
                -cpu)
		cpu=$2
		;;
                -sippath)
                sippath=$2
                ;;
        esac
        shift
done
if [ $hicpath == 0 ]
then
echo "You must specify a .hic file!......""$menu"
exit 0
fi
if [ $myres == 0 ]
then
echo "You must specify a resolution!......""$menu"
exit 0
fi
if [ $sizefile == 0 ]
then
echo "You must specify a path to a chromosome size file!......""$menu"
exit 0
fi
echo "SIPprepr_v""$version"" -straw ""$strawpath"" -hic ""$hicpath"" -res ""$res"" -sizes ""$sizefile"" -matsize ""$matsize"" -cpu ""$cpu"" -outdir ""$outdir"  > SIPprep_paramters.txt
mkdir $outdir


totchroms=`cat $sizefile | wc -l`
coordchunk=`expr $matsize \* $myres`
#coordchunk=`echo $matsize | awk -v var=$myres '{print $1*$myres}'`

echo "coordchunk is ""$coordchunk"
totbins=`cat $sizefile | awk -v var=$myres '{print $1"\t"($2/var)**2}' | awk '{sum += $2} END {print sum}'`
stringres=`echo "$myres" | awk '{print ($1/1000)"kb"}'`
echo "Now dumping reads and calculating distance normalized scores for all chromosomes at ""$myres"
mkdir $outdir/$stringres
normout=`echo "$outdir""/""$stringres""/normVector"`
mkdir $normout
task(){
#echo $myiter
mychr=`cat $sizefile | head -n $myiter | tail -1 | cut -f 1`
#numbins=`cat $sizefile | head -n $myiter | tail -1 | cut -f 2 | awk -v var=$myres '{print int($2/var)+1}'`
mysize=`cat $sizefile | head -n $myiter | tail -1 | cut -f 2`
tmpnumbins=`expr $mysize / $myres`
numbins=`expr $tmpnumbins + 1`
newoutdir=`echo "$outdir""/""$stringres""/""$mychr"`
mkdir $newoutdir
echo "getting distance normalized reads for ""$mychr"
efile=`echo "expecteds"`
ofile=`echo "observeds"`


$strawpath observed $norm $hicpath $mychr $mychr BP $myres | grep -v WARN | grep -i -v Na | awk -v var=$myres -v ms=$coordchunk '{if ($2-$1 <= ms) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' | sort -k 1b,1b --stable > $newoutdir/$ofile &
wait

errorfile=`echo "errorlog"`
cat $newoutdir/$ofile | groupBy -i stdin -g 1 -c 4 -o sum | awk '{print $1"\t"$2"\t"NR-1}' | awk -v var=$numbins '{if ((var-$3) > 0) print $1"\t"$2/(var-$3)}' > $newoutdir/$efile
nfile=`echo "distnorms"`
wait
join -j 1 $newoutdir/$ofile $newoutdir/$efile | sed 's/ /\t/g' | awk '{print $2"\t"$3"\t"$4-$5"\t"($4+1)/($5+1)}' > $newoutdir/$nfile 
wait

#rm $newoutdir/$ofile
#rm $newoutdir/$efile

echo "chunking ""$mychr"
#curchunk=`echo "$coordchunk"`
#numints=`cat $sizefile | head -n $myiter | tail -1 | awk -v var=$matsize -v mres=$myres '{print int($2/(var*mres))+1}'`
tmpnumints=`expr $mysize / $coordchunk`
numints=`expr $tmpnumints + 1`
echo "number of chunks for ""$mychr"" is ""$numints"

for (( start = 0; start <= $numints; start++ ))
do
curstart=`expr $start \* $coordchunk`
tmpstart=`expr $start + 1`
tmpend=`expr $tmpstart \* $coordchunk`
curend=`expr $tmpend - 1`
outchunk=`echo "$newoutdir""/""$mychr""_""$curstart""_""$curend"".txt"`
cat $newoutdir/$nfile | awk -v st=$curstart -v en=$curend '{if (($1 >= st) && ($1 <= en) && ($2 >= st) && ($2 <= en)) printf "%s\t%s\t%.10f\t%.10f\n", $1,$2,$3,$4}' > $outchunk
done
wait
rm $newoutdir/$nfile
rm $newoutdir/$efile
echo "calculating norm vector for ""$mychr"
myVCouttmp=`echo "$normout""/""$mychr"".normtmp"`
lefttmp=`echo "$normout""/""$mychr"".left"`
righttmp=`echo "$normout""/""$mychr"".right"`
zerotmp=`echo "$normout""/""$mychr"".zeros"`
myVCout=`echo "$normout""/""$mychr"".norm"`
cat $newoutdir/$ofile | awk '{if ($3 > $2) print $3"\t"$4}' | groupBy -i stdin -g 1 -c 2 -o sum | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o sum > $righttmp
cat $newoutdir/$ofile | awk '{print $2"\t"$4}' | groupBy -i stdin -g 1 -c 2 -o sum | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o sum > $lefttmp
cat $sizefile | awk -v var=$myres '{for (i=0;i<=$2;i+=var) print i"\t"0}' > $zerotmp
wait
cat $lefttmp $rightmp | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o sum > $myVCouttmp
wait
mymed=`cat $myVCouttmp | awk '{print 1"\t"$2}' | groupBy -i stdin -g 1 -c 2 -o median | cut -f 2`
cat $zerotmp $myVCouttmp | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o sum | awk -v var=$mymed '{print $2/var}' > $myVCout
wait

rm $zerotmp
rm $newoutdir/$ofile
rm $lefttmp
rm $righttmp
rm $myVCouttmp
}

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    
    printf '%.3d' $? >&3
    )&
}


open_sem $cpu

for (( myiter=1; myiter<=$totchroms; myiter++ ))
do
    run_with_lock task $myiter
done

wait
echo "Done Dumping!"

if [ $sippath == 0 ]
then
exit 0
fi

java -jar $sippath processed $outdir $sizefile $outdir -res $myres -mat $matsize -cpu $cpu -d $d -t $t -min $min -max $max -g $g -sat $sat -nbZero $nbzero -del $del -isDroso $isDroso -fdr $fdr

