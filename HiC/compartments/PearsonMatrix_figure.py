#!/usr/bin/env python
# coding: utf-8

import sys, getopt
import math
import matplotlib.pyplot as plt
import time
#import strawC
import hicstraw
import argparse
import numpy as np
import seaborn as sns
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Tool used to create images of the pearson correlation matrix \n\n', formatter_class=argparse.RawTextHelpFormatter)

#parser.add_argument("-p", "--periphery_hicfile", dest="pfile", help="Input periphery .hic file\n\n", metavar="pfile", required=True)
parser.add_argument("-i", "--input_hicfile", dest="ifile", help="Input .hic file\n\n", metavar="ifile", required=True)
parser.add_argument("-p", "--position", dest="chrposition", help="Region to plot. Format is chr:start:end, e.g. 2L:10000000:15000000\n\n", metavar="chrposition", required=True)
#parser.add_argument("-n", "--numcomps", dest="numcomps", default=10, help="Number of distances to plot, default is 10\n\n", metavar="numcomps", required=False)
parser.add_argument("-o", "--output", dest="output", default="PearsonPlot", help="Prefix for output figure, default is PearsonPlot\n\n", metavar="output_prefix", required=False)
parser.add_argument('-r', "--res", dest='res', help="Resolution to process the Hi-C file\n\n", metavar="res", required=True)
#parser.add_argument('-j', "--juicer", dest='juice', help="path to juicer_tools\n\n", metavar="juice", required=True)
#parser.add_argument('-c', "--chrom", dest='chrom', help="Chromosome to process\n\n", metavar="chrom", required=True)
#parser.add_argument('-t', "--thresh", dest='threshold', default=10, help="Filter all distance normalized values above this number. i.e. filter outliers, default is 10. Disable filtering by setting this option to 0.\n\n", metavar="thresh", required=False)
#parser.add_argument('-f', "--fig", dest='createfig', action='store_true', help="Set this to create figures, otherwise will save scaled average and median matrices for later figure creation and customization in python using numpy.loadtxt \n\n")
#parser.add_argument('-d', "--distance", dest='distthresh', default=0, help="Optional distance threshold. Set this to only examine interactions within a certain distance of the diagonal. Default is 0, which does not filter for distance\n\n", required=False, metavar="distthresh")
#parser.add_argument('-f', "--filtdistance", dest='distfilt', default=0, help="Optional upper distance threshold. Set this to only examine interactions further than a certain distance from the diagonal. Default is 0, which doesn't filter by distance.\n\n", required=False, metavar="upperdist") 
#parser.add_argument('-c', "--color", dest='colorscale', default="Reds", help="Matplotlib color gradient to use with the -f option. Default is RdBu_r\n\n",metavar="colorscale")
parser.add_argument('-l', "--lower", dest='lower', default=-1, help="Lower bounds in colorscale, default is -1", required=False, metavar="lower")
parser.add_argument('-u', "--upper", dest='upper', default=1, help="Upper bounds in colorscale, default is 1", required=False, metavar="upper") 
#parser.add_argument('-d', "--decay", dest='decay', action='store_true', help="Set this to include diagonal decay, i.e. do not distance normalize.\n\n")
args = parser.parse_args()

#pfile = str(args.pfile)
#ifile = str(args.ifile)

#eigenfile = str(args.eigen_file)

#numcomps = int(args.numcomps)
myposition = str(args.chrposition)
myhic = str(args.ifile)
myres = int(args.res)
#myjuicer = str(args.juice)
#color = str(args.colorscale)
mymin = float(args.lower)
mymax = float(args.upper)
#outfile = str(args.output) + "_" + str(myposition) + ".png"

mycoords=myposition.split(":")
outfile = str(args.output) + "_" + str(mycoords[0]) + "_" + str(mycoords[1]) + "_" + str(mycoords[2]) + ".png"
mystart=int(mycoords[1])
offset=int(float(mystart)/float(myres))
myend=int(mycoords[2])
mylength=myend-mystart
numbins = int(float(mylength)/float(myres)+1)
myexp = np.zeros(numbins)
newexp = []
#matsize = int((in)/(myres+1)
dense = np.zeros((numbins,numbins))
resultstmp = hicstraw.straw('observed','VC_SQRT', myhic, mycoords[0], mycoords[0], 'BP', myres)
for i in range(len(resultstmp)):
    dist = int(float(resultstmp[i].binY - resultstmp[i].binX)/float(myres))
    if dist < numbins:
        score = float(resultstmp[i].counts)
        myexp[dist] = myexp[dist] + score
    else:
        continue        
    #myexp[dist] = myexp[dist] + score
for j in range(len(myexp)):
    newexp.append(myexp[j]/((numbins)-j))
result = hicstraw.straw('observed','VC_SQRT', myhic, myposition, myposition, 'BP', myres)
mysum=0
for i in range(len(result)):
    score = result[i].counts
    dist = int(float(int(result[i].binY) - int(result[i].binX))/float(myres))
    newscore= (float(score)+1)/(float(myexp[dist])+1)
    left = int(float(int(result[i].binX))/float(myres))-offset
    right = int(float(int(result[i].binY))/float(myres))-offset
    dense[left,right] = newscore
    dense[right,left] = newscore
np.seterr(invalid='ignore')
mycor = np.corrcoef(dense)
f, ax = plt.subplots(figsize=(18, 12))
am = sns.heatmap(mycor,square=True,cmap="bwr",vmin=mymin,vmax=mymax)
f.tight_layout()
avefig = am.get_figure()
avefig.savefig(outfile)


#expecteds = ! java -jar  dump expected VC_SQRT $interactionfile $mychr BP $intres | grep -v HiC | grep -v WARN



#mychrom = str(args.chrom)

#myfilt = float(args.threshold)
#mydist = int(args.distthresh)
#myupfilt = int(args.distfilt)
#mybindist = int(float(mydist)/float(resolution))
#myubindist = int(float(myupfilt)/float(resolution))
#mylow = float(args.lower)
#myup = float(args.upper)
#if mylow == "None":
#    autocheck = 1
#else:
#    mylow = float(mylow)
#if myup == "None":
#    autocheck2 = 1
#else:
#    myup = float(myup)

#First we get the quantiles
"""distnorm="no"
chrlist = [0]
complist = [0]
compchrlist = [0]
compstartlist =[0]
i=0
myslop=0
comps = open(eigenfile,'r')
for line in comps:
    li = line.split()
    if str(mychrom) != str(li[0]):
        continue
    if float(li[3]) > 0:
        compstatus = "A"
    elif float(li[3]) < 0:
        compstatus = "B"
    else:
        continue
    #print(i)
    #print(complist[i-1])
    #print(compstatus)
    #initialize 
    if i == 0:
        complist[i] =(compstatus)
        compchrlist[i]=str(li[0])
        compstartlist[i]=str(li[1])
        i+=1
        myslop=1
        complist.append(0)
        compchrlist.append(0)
        compstartlist.append(0)
        continue
    elif compstatus == complist[i-1]:
        myslop=0
        continue
    elif myslop == 0:
        complist[i] =(compstatus)
        compchrlist[i]=str(li[0])
        compstartlist[i]=str(li[1])
        complist.append(0)
        compchrlist.append(0)
        compstartlist.append(0)
        i+=1
        myslop=1
    else:
        i-=1
        myslop=0
for i in range(0, len(complist)):
    if complist[i] == 0:
        break
    else:
        myposition = i
complist = complist[0:myposition-1]
compchrlist = compchrlist[0:myposition-1]
compstartlist = compstartlist[0:myposition-1]

if distnorm=="yes":
    largestbin=0

    result = strawC.straw('NONE', ifile, str(mychrom), str(mychrom), 'BP',myres)
    for i in range(len(result)):
        #print(result[i].counts)
        if int(result[i].binY) > largestbin:
            largestbin = int(result[i].binY)

    mynumbins = int(float(largestbin)/float(myres))
    totpossible=[]
    for eachitem in range(mynumbins+1,1,-1):
        totpossible.append(eachitem)
    expecteds=[0.0001]*mynumbins
    result = strawC.straw('NONE', ifile, str(mychrom), str(mychrom), 'BP',myres)
    for i in range(len(result)):
        distance = int(float(int(result[i].binY)-int(result[i].binX))/float(myres))
        try:
            myscore = float(result[i].counts) / float(totpossible[distance])
        except IndexError as error:
            #print(distance)
            #print(len(totpossible))
            myscore = float(result[i].counts)
        expecteds[distance] += float(myscore)
else:
    print("...")


myboxes=[]
myboxes2=[]
myAAboxes=[]
myABboxes=[]
myBBboxes=[]
myrelAAboxes=[]
myrelBBboxes=[]
myboxcompA=[]
myboxcompB=[]
myboxcompAB=[]
myrelboxcomp=[]
ncomps=numcomps
for d in range(1,ncomps+1):
    print("Processing Comps with Separation of " + str(d))
    AA=[]
    AB=[]
    BB=[]
    relAA = []
    relBB = []
    AA2=[]
    AB2=[]
    BB2=[]
    relAA2 = []
    relBB2 = []
    AA3=[]
    AB3=[]
    BB3=[]
    relAA3=[]
    relBB3=[]
    for h in tqdm(range(0,(len(complist)-((d*2)+2)))):
        if compchrlist[h] == compchrlist[h+(d*2)+2]:
            if float(int(compstartlist[h+1]) - int(compstartlist[h]))/float(myres) >= 1000:
                continue
            loc1 = str(mychrom) + ":" + str(compstartlist[h]) + ":" + str(compstartlist[h+1])
            samecomp = str(mychrom) + ":" + str(compstartlist[h+(d*2)]) + ":" + str(compstartlist[h+(d*2)+1])
            diffcomp = str(mychrom) + ":" + str(compstartlist[h+(d*2)+1]) + ":" + str(compstartlist[h+(d*2)+2])
            numbinsi = float(int(compstartlist[h+1]) - int(compstartlist[h]))/float(myres)
            numbinsj = float(int(compstartlist[h+2+d]) - int(compstartlist[h+1+d]))/float(myres)
            numbinsk = float(int(compstartlist[h+3+d]) - int(compstartlist[h+2+d]))/float(myres)
            posssame = numbinsj*numbinsi
            possdiff = numbinsk*numbinsi
            mysum=0
            result = strawC.straw('NONE', ifile, loc1, samecomp, 'BP', myres)
            mysum=0
            if len(result) == 0:
                continue
            for i in range(len(result)):
                score = result[i].counts
                
                if distnorm=="yes":
                    dist = int(float(int(result[i].binY) - int(result[i].binX))/float(myres))
                    try:
                        newscore= (float(score)+1)/(float(expecteds[dist])+1)
                    except IndexError:
                        print(dist)
                else:
                    newscore=score
                mysum+=newscore
            try:
                mymean = mysum/posssame
            except ZeroDivisionError:
                #print(len(result))        
                continue
            #for i in range(len(result)):
            #    mysum += result[i].countsprint("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
            if complist[h] == "A":
                AA.append(math.log(mymean))
            else:
                BB.append(math.log(mymean))
                
            result = strawC.straw('NONE', pfile, loc1, samecomp, 'BP', myres)
            mysum=0
            if len(result) == 0:
                continue
            for i in range(len(result)):
                score = result[i].counts
                
                if distnorm=="yes":
                    dist = int(float(int(result[i].binY) - int(result[i].binX))/float(myres))
                    try:
                        newscore= (float(score)+1)/(float(expecteds[dist])+1)
                    except IndexError:
                        print(dist)
                else:
                    newscore=score
                mysum+=newscore
            try:
                mymean2 = mysum/posssame
            except ZeroDivisionError:
                #print(len(result))        
                continue
            #for i in range(len(result)):
            #    mysum += result[i].countsprint("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
            if complist[h] == "A":
                AA2.append(math.log(mymean2))
                AA3.append(math.log(mymean/mymean2)*-1)
            else:
                BB2.append(math.log(mymean2))
                BB3.append(math.log(mymean/mymean2)*-1)
            
                
            result = strawC.straw('NONE', ifile, loc1, diffcomp, 'BP', myres)

            mysum=0
            if len(result) == 0:
                continue
            for i in range(len(result)):
                score = result[i].counts
                if distnorm=="yes":
                    dist = int(float(int(result[i].binY) - int(result[i].binX))/float(myres))
                    newscore= (float(score)+1)/(float(expecteds[dist])+1)
                else:
                    newscore=score
                mysum+=newscore
            try:
                mymeanAB = mysum/possdiff
            except ZeroDivisionError:
                #print(len(result))
                continue
            AB.append(math.log(mymeanAB))
            AB3.append(math.log(mymeanAB)*-1)
            if complist[h] == "A":
                relAA.append(math.log(mymean/mymeanAB))
            else:
                relBB.append(math.log(mymean/mymeanAB))
            
            result = strawC.straw('NONE', pfile, loc1, diffcomp, 'BP', myres)

            mysum=0
            if len(result) == 0:
                continue
            for i in range(len(result)):
                score = result[i].counts
                if distnorm=="yes":
                    dist = int(float(int(result[i].binY) - int(result[i].binX))/float(myres))
                    newscore= (float(score)+1)/(float(expecteds[dist])+1)
                else:
                    newscore=score
                mysum+=newscore
            try:
                mymeanAB2 = mysum/possdiff
            except ZeroDivisionError:
                #print(len(result))
                continue
            AB2.append(math.log(mymeanAB2))
            AB3.append(math.log(mymeanAB2))
            if complist[h] == "A":
                relAA2.append(math.log(mymean2/mymeanAB2))
                relAA3.append(math.log((mymean/mymeanAB)/(mymean2/mymeanAB2))*-1)
            else:
                relBB2.append(math.log(mymean/mymeanAB))
                relBB3.append(math.log((mymean/mymeanAB)/(mymean2/mymeanAB2))*-1)
    myAAboxes.append(AA)
    myAAboxes.append(AA2)
    myABboxes.append(AB)
    myABboxes.append(AB2)
    myBBboxes.append(BB)
    myBBboxes.append(BB2)
    myrelAAboxes.append(relAA)
    myrelAAboxes.append(relAA2)
    myrelBBboxes.append(relBB)
    myrelBBboxes.append(relBB2)
    myboxcompA.append(AA3)
    myboxcompAB.append(AB3)
    myboxcompB.append(BB3)
    #myboxcompAB.append(AB3)
    myrelboxcomp.append(relAA3)
    myrelboxcomp.append(relBB3)

fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
 
# Creating axes instance
bp = ax.boxplot(myboxcompA, patch_artist = True, showfliers=False)
colors=[]
for i in range(0,ncomps):
    
    colors.append('#FF6600')
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)
plt.axhline(0, c='blue')
ax.set_ylim(-4,4)
outA=str(outfile) + "_A.png"
fig.savefig(outA)

fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
 
# Creating axes instance
bp = ax.boxplot(myboxcompB, patch_artist = True, showfliers=False)
colors=[]
for i in range(0,ncomps):
    
    colors.append('#0099ff')
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)
plt.axhline(0, c='blue')
ax.set_ylim(-4,4)

outB=str(outfile) + "_B.png"
fig.savefig(outB)


#myeig = pd.read_table(eigenfile,delimiter='\t',names=["chrom","start","end","eig"]).sort_values("eig")
#myeig = myeig.reset_index(drop=True)
#eigsize = (len(myeig))

# estimate expecteds
chromtable = myeig.groupby(['chrom']).max().sort_values("end", ascending=False)
chromlist=list(chromtable.index)
sizelist = chromtable['end'].tolist()
bigsize = sizelist[0]
bigchrom = chromlist[0]
numbins = int(float(bigsize)/float(resolution)+1)
myexp = np.zeros(numbins)
newexp = []
resultstmp = strawC.straw('NONE', interactionfile, bigchrom, bigchrom, 'BP', resolution)
for i in range(len(resultstmp)):
    dist = int(float(resultstmp[i].binY - resultstmp[i].binX)/float(resolution))
    score = float(resultstmp[i].counts)
    myexp[dist] = myexp[dist] + score
for j in range(len(myexp)):
    newexp.append(myexp[j]/((numbins)-j))

#get quantile
qmap = np.zeros((bins,bins))
countmap= np.ones((bins,bins))
print
for c in range(0,len(chromlist)):
       
    curchrom = chromlist[c]
    myeigtmp = myeig[myeig["chrom"] == curchrom]
    myeigtmp.reset_index(inplace=True)
    myeigsize = int(float(bigsize)/float(resolution))+1
    mytmplist = [-1]*myeigsize
    for index, row in myeigtmp.iterrows():
        mypos = row['start']
        myindx = row['index']
        mybin = int(float(mypos)/float(resolution))
        mytmplist[mybin] = int(myindx)
        
    resultstmp = strawC.straw('NONE', interactionfile, curchrom, curchrom, 'BP', resolution)
    for i in tqdm(range(len(resultstmp)),desc="Calculating chromosome " + str(chromlist[c])):
        dist = int(float(resultstmp[i].binY - resultstmp[i].binX)/float(resolution))
        if (mybindist > 0):
            if (dist > mybindist):
                continue
            else:
                score = float(resultstmp[i].counts+1)/float(newexp[dist]+1)
        elif (myubindist > 0):
            if (dist < myubindist):
                continue
            else:
                score = float(resultstmp[i].counts+1)/float(newexp[dist]+1)
        else:
            score = float(resultstmp[i].counts+1)/float(newexp[dist]+1)
        if (score > myfilt) and (myfilt > 0):
            continue
        Xpos = int(float(resultstmp[i].binX)/float(resolution))
        Ypos = int(float(resultstmp[i].binY)/float(resolution))        
        if (mytmplist[Xpos] > 0):
            try:
                percentileX = int(math.floor((mytmplist[Xpos])/float(eigsize)*bins))
            except IndexError:
                continue
        else:
            continue
        if (mytmplist[Ypos] > 0):
            try:
                percentileY = int(math.floor((mytmplist[Ypos])/float(eigsize)*bins))
            except IndexError:
                continue
        else:
            continue
        qmap[percentileX,percentileY] = qmap[percentileX,percentileY] + score
        countmap[percentileX,percentileY] = countmap[percentileX,percentileY] + 1

#save matrix
outmat = str(outfile) + "_mat.txt"
normmap = np.divide(qmap,countmap)
np.savetxt(outmat, normmap)

#save figure
if (mylow == 0) and (myup == 0):
    f, ax = plt.subplots(figsize=(9, 6))
    am = sns.heatmap(np.divide(qmap,countmap),square=True,cmap=color,vmin=None,vmax=None)
    f.tight_layout()
    avefig = am.get_figure()
    avefig.savefig(str(outfile) + "_heat.png")
else:
    f, ax = plt.subplots(figsize=(9, 6))
    am = sns.heatmap(np.divide(qmap,countmap),square=True,cmap=color,vmin=mylow,vmax=myup)
    f.tight_layout()
    avefig = am.get_figure()
    avefig.savefig(str(outfile) + "_heat.png")


"""
