for i in data/hic/INTS11auxSIP/ data/hic/INTS1noauxSIP/ data/hic/INTS11noauxSIP/
do
for j in ctcf_loops_stable_raw2fold.bedpe ctcf_loops_down_raw2fold.bedpe ctcf_loops_up_raw2fold.bedpe
do
degron=`echo "$i"  | sed 's/\//\t/g' | awk '{print "filt_"$3}'`
java -jar tools/SIPMeta_1.3.jar simple $j $i tools/bullseye.py 21 2000 -min 0 -max 5 -res 25000 -prefix $degron -c Red
#java -jar tools/SIPMeta_1.3.jar simple ctcf_loops_stable_raw2fold.bedpe data/hic/INTS11auxSIP/ tools/bullseye.py 21 2000 -min 0 -max 5 -res 25000 -prefix filt_tot_ints11 -c Red
#java -jar tools/SIPMeta_1.3.jar simple ctcf_loops_stable_raw2fold.bedpe data/hic/INTS11auxSIP/ tools/bullseye.py 21 2000 -min 0 -max 5 -res 25000 -prefix filt_tot_ints11 -c Red
done
done

