import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys, getopt
import time
import argparse
import seaborn as sns

parser = argparse.ArgumentParser(description='Tool used to create scaled 1D metaplots (similar to Scaled2DMetaplotter)\n\n', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i", "--bwfile", dest="bw_file", help="Input .bw file\n\n", metavar="bw_file", required=True)
parser.add_argument("-a", "--anchorfile", dest="anchor_file", help="Input anchor file. Format: Chromosome Anchor(s) Strand\n\n", metavar="anchor_file", required=True)
parser.add_argument("-b", "--bins", dest="numbins", default=10, 
                    help="Number of pixels to plot between each anchor, default is 10 pixels. (Can be a list)\n\n", metavar="bins")
parser.add_argument("-o", "--output", dest="output", default="meta_out", help="Prefix for output scaled matrices, default is meta_out\n\n", metavar="output_prefix")
parser.add_argument('-f', "--fig", dest='createfig', action='store_true', help="Set this to create figures, otherwise will save scaled average and median matrices for later figure creation and customization in python using numpy.loadtxt \n\n")
parser.add_argument('-y', "--ylims", dest='ylims', default='0,0', help="Two comma separated values for the desired min and max on the y-axis when creating figures with the -f option. \n\n", metavar="ylims")
parser.add_argument('-l', "--linewidth", dest='lw', default="3", help="Linewidth for the plot. Default is 3. \n\n", metavar="lw")
parser.add_argument('-c', "--color", dest='color', default="black", help="Matplotlib color to use with the -f option. Default is black\n\n",metavar="color")
args = parser.parse_args()

#passing arguments
bwfile = str(args.bw_file)
anchorfile = str(args.anchor_file)

binstr = args.numbins
if ',' in binstr:
    binlist = [int(x) for x in binstr.split(',')]
else:
    binlist = None
    binnum = int(binstr)

outputname = str(args.output)
inputmm = str(args.ylims)
myminmax = inputmm.split(",")
mymin = float(myminmax[0])
mymax = float(myminmax[1])

#initialized variables
bw = pyBigWig.open(bwfile)
anchorlist=[]
linenum=0


with open(anchorfile, 'r') as f:
    for i, line in enumerate(tqdm(f)):
        li=line.split()
        #listpos=0
        strand = str(li[-1])
        if linenum == 0:
            numintervals = len(li) - 3
            if binlist:
                if len(binlist) != numintervals:
                    print(f"Error: You specified {len(binlist)} bins, but file implies {numintervals} regions.")
                    sys.exit(1)
                totalbins = sum(binlist)
            else:
                totalbins = numintervals * binnum

            mymat = [[] for b in range(totalbins)]
            linenum+=1

        else:
            newnumintervals = len(li) - 3
            if newnumintervals != numintervals:
                print("inconsistent number of columns on! exiting at line: ")
                print(linenum)
                break
            linenum+=1
        if strand == "-":
            listpos=totalbins-1
        else:
            listpos=0

        #ensure distances are greater than 1
        distancechecker=[]
        for i in range(0,numintervals):
            distancechecker.append(int(li[i+2])-int(li[i+1]))
        mindist = min(distancechecker)
        if binlist:
            minbins = min(binlist)
        else:
            minbins = binnum

        if mindist <= minbins:
            continue

        for i in range(numintervals):
            intervaldistance = int(li[i+2]) - int(li[i+1])
            nbins = binlist[i] if binlist else binnum
            truebinsize = float(intervaldistance) / float(nbins)

            for j in range(nbins):
                start = int(li[i+1])+ round(float(truebinsize)*j)
                end = int(li[i+1]) + round(float(truebinsize)*(j+1))
                try:
                    mysignal = float(bw.stats(li[0], int(start), int(end))[0])
                except:
                    mysignal = 0
                mymat[listpos].append(mysignal)
                #print(listpos)
                if strand == "-":
                    listpos-=1
                else:
                    listpos+=1
avelist=[]
for x in mymat:
    avelist.append(np.nanmean(x))
medlist=[]
for x in mymat:
    medlist.append(np.nanmedian(x))

with open(str(outputname)+"ave_med.txt", 'w') as outie:
    outie.write("averages\t")
    for item in avelist:
        outie.write(str(item) + "\t")
    outie.write('\n')
    outie.write("medians\t")
    for item in medlist:
        outie.write(str(item) + "\t")
    outie.write('\n')
myxmat = np.transpose(np.array(mymat))
np.savetxt(str(outputname)+"mat.txt", myxmat)


if args.createfig:
    lw = float(args.lw)
    mycolor = str(args.color)
    avefig=str(outputname)+"ave.png"
    medfig=str(outputname)+"med.png"
    plt.plot(avelist, linewidth=lw, color=mycolor)
    if mymax > mymin:
        plt.ylim(mymin,mymax)
    plt.savefig(avefig)
    plt.close()
    plt.plot(medlist, linewidth=lw, color=mycolor)
    if mymax > mymin:
        plt.ylim(mymin,mymax)
    plt.savefig(medfig)
    plt.close()

    myheatmatrix = sns.heatmap(myxmat, cmap="Reds")
    figure=myheatmatrix.get_figure()
    figure.savefig(str(outputname)+"heatmap.png")

