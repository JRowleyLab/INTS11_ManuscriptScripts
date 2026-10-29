python tools/Scaled2DMetaplotter2_straw.py -i data/hic/INT1noAux.hic -a over10kb_6anchors.txt -b 20 -o over10kb_ctl_1kb -r 1000 -j straw -f &
python tools/Scaled2DMetaplotter2_straw.py -i data/hic/INTS11aux.hic -a over10kb_6anchors.txt -b 20 -o over10kb_aux_1kb -r 1000 -j straw -f &
python tools/Scaled2DMetaplotter2_straw.py -i data/hic/INTS11noAux.hic -a over10kb_6anchors.txt -b 20 -o over10kb_noaux_1kb -r 1000 -j straw -f &
wait
