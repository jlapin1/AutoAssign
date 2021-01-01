#! /usr/bin/python3
import struct
import numpy as np
import os
import sys
import math
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator

class params():
    cmax = 0
    cmin = 0
    rmax = 0
    rmin = 0
    csz = 0
    rsz = 0
    caxis = 0
    raxis = 0
    grid = 0
    filename = 'file'
    def set_flce(self, floor, ceiling):
        rows,cols = np.shape(self.grid)
        for i in range(rows):
            for j in range(cols):
                if self.grid[i,j] < floor:
                    self.grid[i,j] = floor
                if self.grid[i,j] > ceiling:
                    self.grid[i,j] = ceiling

class DynamicUpdate():
    
    def __init__(self):
        self.fig = 0
        
    # on_launch(), on_running(), and peakshift() are used together
    def on_launch(self, out, peaks, levels, invy):
        #Set up plot
        self.figure = plt.figure(figsize=(10,10))
        self.ax = plt.subplot(111)
        X,Y = np.meshgrid(out.caxis, out.raxis)
        self.ax.contour(X,Y,out.grid,levels, linewidths=0.5, cmap='Blues')
        self.lines, = self.ax.plot(peaks[:,0],peaks[:,1], 'o', markerfacecolor='none', markeredgecolor='k', markersize=5)
        self.ax.invert_xaxis()
        
        self.ax.xaxis.set_major_locator(MultipleLocator(5))
        self.ax.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax.xaxis.set_tick_params(top=True)
        if max(out.raxis)<20:
            self.ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            self.ax.yaxis.set_major_locator(MultipleLocator(1))
        else:
            self.ax.yaxis.set_major_locator(MultipleLocator(5))
            self.ax.yaxis.set_minor_locator(MultipleLocator(1))
        
        if invy==True:
            self.ax.invert_yaxis()
        txt = []
        for m in range(peaks.shape[0]):
            txt.append(self.ax.text(peaks[m,0],peaks[m,1],'%d'%(m), size=8, color='red'))
        self.ax.set_title('%s'%(out.filename))
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
        plt.show(block=False)
        return txt

    def on_running(self,peaks):
        #Update data (with the new _and_ the old points)
        self.lines.set_xdata(peaks[:,0])
        self.lines.set_ydata(peaks[:,1])
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
    
    def peakshift(self, out, Peaks, levels, invy):
        peaks = copy.copy(Peaks)
        
        tick=0
        while tick==0:
            txt = self.on_launch(out, peaks, levels, invy)
            self.figure.canvas.draw() # To update, this must be executed outside subroutine that does the plotting
            self.figure.canvas.flush_events()
            print("\nSpectrum is normalized such that the most intense point is -1, near 0 is noise")
            I = 'Reset contour levels (y,n)?'
            log.append([I, input("\n"+I+"\n>>> ")])
            ans = log[-1][1]
            if ans.upper()=='Y':
                I = 'Input low,high'
                log.append([I, input('\n'+I+'\n>>> ')])
                [low,high] = [float(m) for m in log[-1][1].split(",")]
                levels = np.linspace(low, high, 50)
                plt.close()
            else:
                tick=1
        tick=0
        while tick==0:
            I = "Do you want to shift any peaks (y,n)?"
            log.append([I,input('\n'+I+'\n>>> ')])
            ans = log[-1][1]
            if ans.upper()=='Y':
                I = "Index (<0 for all)"
                log.append([I,input('\n'+I+'\n>>> ')])
                ans1 = int(log[-1][1])
                I = "Enter x,y"
                log.append([I,input('\n'+I+'\n>>> ')])
                ans2 = np.array([float(m) for m in log[-1][1].split(',')])
                if ans1>=0:
                    peaks[ans1]+=ans2
                    txt[ans1].set_position((txt[ans1].get_position()[0]+ans2[0], txt[ans1].get_position()[1]+ans2[1]))
                else:
                    peaks+=ans2
                    for m in txt:
                        m.set_position((m.get_position()[0]+ans2[0],m.get_position()[1]+ans2[1]))
            else:
                tick=1
            self.on_running(peaks)
            self.figure.savefig("%s.png"%(out.filename), type='png', dpi=200)
        return peaks
    
    def plotoverlay(self, out1, out2, invy):
        self.figure = plt.figure(figsize=(10,10))
        self.ax = plt.subplot(111)
        self.ax.invert_xaxis()
        if invy==True:
            self.ax.invert_yaxis()
        self.plotov(out1, 'Blues')
        self.plotov(out2, 'Reds')
        self.figure.canvas.draw() # To update, this must be executed outside subroutine that does the plotting
        self.figure.canvas.flush_events()
    def plotov(self, out, cmap):
        X,Y = np.meshgrid(out.caxis, out.raxis)
        self.ax.contour(X,Y,out.grid, np.linspace(-1,-0.1,50),linewidths=0.5,cmap=cmap)
        self.figure.canvas.draw() # To update, this must be executed outside subroutine that does the plotting
        self.figure.canvas.flush_events()
        plt.show(block=False)

def readraw(filename, rnyquist, cnyquist, refr, refc, rmin, rmax, cmin, cmax):
    out = params()
    out.filename = filename
    with open(filename,'rb') as f:
        header = struct.unpack('f'*512, f.read(4*512))
        m = int(header[99])
        n = int(header[219])
        spectrum1 = struct.unpack('f'*n*m,f.read(4*n*m))
        S1 = np.reshape(spectrum1, [m,n], order='F')
    CSA10 = refr+np.linspace(rnyquist, -rnyquist, n)#switched n&m
    CSA20 = refc+np.linspace(cnyquist, -cnyquist, m)#switched n&m
    point12 = int(round(1 + ( ((n-1)/2) * (1+((refr-rmin)/rnyquist)) ) )) # magic number
    point11 = int(round(1 + ( ((n-1)/2) * (1+((refr-rmax)/rnyquist)) ) )) # magic number
    point22 = int(round(1 + ( ((m-1)/2) * (1+((refc-cmin)/cnyquist)) ) )) # magic number
    point21 = int(round(1 + ( ((m-1)/2) * (1+((refc-cmax)/cnyquist)) ) )) # magic number
    C1 = S1[point21-1:point22,point11-1:point12]
    C1 = C1.transpose()
    C1 /= np.max(C1)
    C1 *= -1
    
    out.cmax = CSA20[point21-1]
    out.cmin = CSA20[point22-1]
    out.rmax = CSA10[point11-1]
    out.rmin = CSA10[point12-1]
    out.rsz,out.csz = np.shape(C1)
    out.caxis = np.linspace(out.cmax, out.cmin, out.csz)
    out.raxis = np.linspace(out.rmax, out.rmin, out.rsz)
    out.grid = C1
    return out

def genspec(out, peaks, mamp, msigx, msigy):
    spac = 30
    xspac = (out.caxis[0]-out.caxis[-1])/(out.csz-1)
    yspac = (out.raxis[0]-out.raxis[-1])/(out.rsz-1)
    specout = np.zeros((out.rsz, out.csz))
    for k in range(peaks.shape[0]):
        xind = int(round((out.caxis[0]-peaks[k,0])/xspac))
        yind = int(round((out.raxis[0]-peaks[k,1])/yspac))
        for i in range(yind-spac, yind+spac+1):
            for j in range(xind-spac, xind+spac+1):
                if (i>=0) & (i<out.rsz) & (j>=0) & (j<out.csz):
                    specout[i,j] += -mamp*np.exp( -1*((peaks[k,1]-out.raxis[i])/(np.sqrt(2)*msigy))**2 - ((peaks[k,0]-out.caxis[j])/(np.sqrt(2)*msigx))**2)
    return specout

def plot(out, Peaks, invy):
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(111)
    X,Y = np.meshgrid(out.caxis, out.raxis)
    ax.invert_xaxis()
    if invy==True:
        ax.invert_yaxis()
    ax.contour(X,Y,out.grid,np.linspace(np.min(out.grid),-0.1,50),linewidths=0.5,cmap='Blues')
    ax.plot(Peaks[:,0],Peaks[:,1],'ko', markerfacecolor='none')
    fig.canvas.draw() # To update, this must be executed outside subroutine that does the plotting
    fig.canvas.flush_events()
    plt.show(block=False)
    
def xpeaks(pk):
    length = len(pk[0,:])
    xpk = np.zeros((2, 2*(length-1)))
    for i in range(length-1):
        xpk[0][i] = pk[0][i]
        xpk[1][i] = pk[1][i+1]
        xpk[0][i+length-1] = pk[0][i+1]
        xpk[1][i+length-1] = pk[1][i]
    return xpk

def Expk(out, peak1, peak2):
    cspac = out.caxis[0]-out.caxis[1]
    rspac = out.raxis[0]-out.raxis[1]
    cind = round((out.caxis[0]-peak1[0])/cspac).astype('int')
    rind = round((out.raxis[0]-peak2[1])/rspac).astype('int')
    e = out.grid[rind, cind]
    return e

def eng3(out, length, xpks):
    xlength = 2*(length-1)
    rE = np.zeros((xlength))
    rspac = (out.rmax-out.rmin)/float(out.rsz-1)
    cspac = (out.cmax-out.cmin)/float(out.csz-1)
    E = 0
    for i in range(xlength):
        rind = round((out.rmax-xpks[0,i])/rspac)
        cind = round((out.cmax-xpks[1,i])/cspac)
        rE[i] = out.grid[int(rind), int(cind)]
        E += rE[i]
    return E, rE

def writeinput(filename, aE, pkis, res=0):

    if res==0:
        res=len(pks[0,:])*'X'
    with open(filename, 'w') as f:
        # amt
        f.write('%d\n'%(len(aE)))
        
        # Length
        f.write('%d\n'%(len(pkis)))
        
        # Peaks
        hold = [str(i) for i in pkis]
        f.write(",".join(hold)+"\n")
        
        # Residues
        f.write(",".join(res)+'\n')
        
        # Grid
        for i in range(len(aE)):
            for j in range(len(pkis)):
                hold = [str(k) for k in aE[i][j,:]]
                f.write(",".join(hold)+"\n")

def writespec(filename, outob):
    
    with open("Potgrids/"+filename, "w") as f:
        
        # Grid size
        f.write("%d,%d\n"%(outob.rsz,outob.csz))
        
        # Row max/min, col max/min
        f.write("%f,%f,%f,%f\n"%(outob.rmax,outob.rmin,outob.cmax,outob.cmin))
        
        for i in range(outob.rsz):
            f.write(",".join([str(j) for j in outob.grid[i]])+"\n")

def neighbors(ind1, ind2, scale, order1, order2, pks1, pks2, out1, out2):
    xpk11 = np.array([pks1[0,ind1],pks1[1,ind2]])
    xpk12 = np.array([pks1[0,ind2],pks1[1,ind1]])
    xpk21 = np.array([pks2[0,ind1],pks2[1,ind2]])
    xpk22 = np.array([pks2[0,ind2],pks2[1,ind1]])
    
    rspac = (out1.rmax-out1.rmin)/(out1.rsz-1)
    cspac = (out1.cmax-out1.cmin)/(out1.csz-1)
    rind = round((out1.raxis[0]-xpk11[0])/rspac)
    cind = round((out1.caxis[0]-xpk11[1])/cspac)
    out1.grid[int(rind),int(cind)]*=scale
    rind = round((out1.raxis[0]-xpk12[0])/rspac)
    cind = round((out1.caxis[0]-xpk12[1])/cspac)
    out1.grid[int(rind),int(cind)]*=scale
    
    rspac = (out2.rmax-out2.rmin)/(out2.rsz-1)
    cspac = (out2.cmax-out2.cmin)/(out2.csz-1)    
    rind = round((out2.raxis[0]-xpk21[0])/rspac)
    cind = round((out2.caxis[0]-xpk21[1])/cspac)
    out2.grid[int(rind),int(cind)]*=scale
    rind = round((out2.raxis[0]-xpk22[0])/rspac)
    cind = round((out2.caxis[0]-xpk22[1])/cspac)
    out2.grid[int(rind),int(cind)]*=scale

def readtab(fn):
    f = open(fn, 'r')
    header = []
    for m in range(6):
        header.append(f.readline())
    out = [];outall = []
    for line in f:
        buffer = line.split()
        outall.append(buffer)
        out.append([float(buffer[5]),abs(float(buffer[6]))])
    f.close()
    return np.array(out)
        
def readft(filename):
    f = open(filename, "r")
    cnyq=0;rnyq=0;cnyqd=1;rnyqd=1;refc=0;refr=0
    for line in f:
        ln = line.split()
        if len(ln)>0:
            if ln[0]=='-xSW':
                cnyq = float(ln[1])
                rnyq = float(ln[3])
            if ln[0]=='-xOBS':
                cnyqd = float(ln[1])
                rnyqd = float(ln[3])
                if rnyqd==0:
                    rnyqd = 1000
                    print("yOBS is zero -> defaulting to 1000")
            if ln[0]=='-xCAR':
                refc = float(ln[1])
                refr = float(ln[3])
    return cnyq/cnyqd/2, rnyq/rnyqd/2, refc, refr
#################### end functions/subroutines ################################
plt.close('all')

log = []

I = "Type path to working directory"
log.append([I,input(I+'\n>>> ')])
os.chdir(log[-1][1])

# Read in peaks
files = os.listdir("Order/")
for m in files:
    if m[-3:]=='tab':
        nh = readtab("Order/"+m)
        break
sort = np.argsort(nh[:,0])
nh = nh[sort]
nn = np.column_stack((nh[:,0],nh[:,0]))
csmin = min(nh[:,0])-5
csmax = max(nh[:,0])+5
nhmin = min(nh[:,1])-0.5
nhmax = max(nh[:,1])+0.5

# Read in spectra
files = os.listdir('Spectra/')
fid=[]
fnd = {
       'NNY':False,
       'NNX':False,
       'NHY':False,
       'NHX':False
       }
for m in files:
    if m[-3:]=='.ft':
        if m[-6:-3]=='NNY':
            fidNNY = m
            fnd['NNY'] = True
            print('Found NN noXC')
            cnyq,rnyq,refc,refr = readft('Spectra/NNY.com')
            NNYo = readraw('Spectra/'+m, rnyq, cnyq, refr, refc, csmin, csmax, csmin, csmax)
        elif m[-6:-3]=='NNX':
            fidNNX = m
            fnd['NNX'] = True
            print('Found NN XC')
            cnyq,rnyq,refc,refr = readft('Spectra/NNX.com')
            NNXo = readraw('Spectra/'+m, rnyq, cnyq, refr, refc, csmin, csmax, csmin, csmax)
        elif m[-6:-3]=='NHY':
            fidNHY = m
            fnd['NHY'] = True
            print('Found NH noXC')
            cnyq,rnyq,refc,refr = readft('Spectra/NHY.com')
            NHYo = readraw('Spectra/'+m, rnyq, cnyq, refr, refc, nhmin, nhmax, csmin, csmax)
        elif m[-6:-3]=='NHX':
            fidNHX = m
            fnd['NHX'] = True
            print('Found NH XC')
            cnyq,rnyq,refc,refr = readft('Spectra/NHX.com')
            NHXo = readraw('Spectra/'+m, rnyq, cnyq, refr, refc, nhmin, nhmax, csmin, csmax)
        fid.append(m)

# Read in sequence
for m in os.listdir("Order"):
    if m[-3:]=='seq':
        with open("Order/"+m,'r') as f:
            f.readline()
            resraw = f.readline().strip()
        break
resarr = np.array(" ".join(resraw).split())

# Read in selective labeling file
isc = [];iscz=[];sl=[]
with open("Order/labs.sl","r") as f:
    for line in f:
        line=line.strip()
        if len(line)==1: # Line is a header for the amino acid label
            if line!='Z':
                sl.append(line) # add amino acid to variable
                cnt = 0
                for m in range(sum(resarr==line)): # loop through number of "line" residues in sequence
                    pk = [float(n) for n in f.readline().split()]
                    while cnt < len(resarr): # Find indices of all <line> residues in sequence
                        if resarr[cnt]==line: # Found "line" residue at index "cnt"
                            isc.append([cnt,line,pk[0],pk[1]]) # Save index, amino acid type, and shift/coup
                            cnt+=1
                            break
                        cnt+=1
            elif line=='Z':
                line2 = f.readline().split()
                while line2!=['END']:
                    pk = [float(n) for n in line2[1:]]
                    iscz.append([ int(line2[0]), pk[0], pk[1]])
                    line2 = f.readline().split()

# Create a starting permutation
restyp = list(len(resraw)*'X') # create amino acid list of all Xs
NNH = -1*np.ones((nn.shape[0]),'int')
for m in range(len(isc)):
    # Find current index of desired shift/coup
    low=[0,1E10]
    for n in range(nn.shape[0]):
        dif = np.sqrt((isc[m][2]-nh[n,0])**2 + (isc[m][3]-nh[n,1])**2)
        if dif<low[1]:
            low = [n,dif]
    # Make switch
    NNH[isc[m][0]] = low[0]
    restyp[isc[m][0]] = isc[m][1] # Set the index to the amino acid type
for m in iscz:
    # Find current index of desired shift/coup
    low=[0,1E10]
    for n in range(nn.shape[0]):
        dif = np.sqrt((m[1]-nh[n,0])**2 + (m[2]-nh[n,1])**2)
        if dif<low[1]:
            low = [n,dif]
    # Make switch
    NNH[m[0]] = low[0]
    restyp[m[0]] = 'Z'
# Put the remaining unused indices into the "-1" designation 
for m in np.arange(nn.shape[0]):
    if m not in NNH:
        for n,o in enumerate(NNH):
            if o==-1:
                NNH[n] = m
                break
restyp = "".join(restyp)

permtot=1
hold,counts = np.unique(list(restyp),return_counts=True)
for i in counts:
    permtot*=math.factorial(i)
print("\n%s\n%d total order permutations"%(restyp,permtot))

# Editing peak positions
print("\nEdit peak positions")
plt.close('all')
levels1 = np.linspace(-1, -0.05, 50)
d = DynamicUpdate() # Functions contained in a class object just because
nn2 = d.peakshift(NNXo, nn, levels1, invy=True) # NN plot
nh2 = d.peakshift(NHYo, nh, levels1, invy=False) # NH plot
#nh3 = d.peakshift(NHXo, nh2, levels2, invy=False)
plt.close('all')


print("\nEdit NN spectra")
syn1 = copy.deepcopy(NNXo)
pki = np.zeros((nn.shape[0]))
for m in range(nn.shape[0]):
    cind = int(round((syn1.caxis[0]-nn[m][0])/(syn1.caxis[0]-syn1.caxis[1])))
    rind = int(round((syn1.raxis[0]-nn[m][1])/(syn1.raxis[0]-syn1.raxis[1])))
    pki[m] = syn1.grid[rind,cind]
tick=0
while tick==0:
    print("\nIntensity |avg| std: %.3f %.3f\nsig: ~1 for ppm, ~0.15 for kHz"%(abs(np.mean(pki)), np.std(pki)))
    I = "Input gaussian amp,sigx,sigy"
    log.append([I,input('\n'+I+'\n>>> ')])
    [a,sx,sy] = [float(m) for m in log[-1][1].split(",")]
    syn1.grid = genspec(syn1, nn2, a, sx, sy)
    d.plotoverlay(NNXo, syn1, invy=True)
    input("\nPress enter to continue to edited plot")
    
    outnn = copy.deepcopy(NNXo)
    outnn.grid -= syn1.grid
    outnn.set_flce(np.min(outnn.grid), 0)
    outnn.grid /= abs(np.min(outnn.grid))
    
    levels1 = np.linspace(np.min(outnn.grid), -0.1, 50)
    plot(outnn, nn2, invy=True)
    plt.gcf().canvas.draw() # To update, this must be executed outside subroutine that does the plotting
    plt.gcf().canvas.flush_events()
    
    I = "Happy (y/n)?"
    log.append([I,input('\n'+I+'\n>>> ')])
    ans = log[-1][1]
    if ans.upper()=='Y':
        tick=1
    else:
        plt.close()
plt.close('all')

print("\nEdit NH spectra")
syn2 = copy.deepcopy(NHXo)
pki = np.zeros((nn.shape[0]))
for m in range(nn.shape[0]):
    cind = int(round((syn2.caxis[0]-nh[m][0])/(syn2.caxis[0]-syn2.caxis[1])))
    rind = int(round((syn2.raxis[0]-nh[m][1])/(syn2.raxis[0]-syn2.raxis[1])))
    pki[m] = syn2.grid[rind,cind]
tick=0
while tick==0:
    if fnd['NHY']==True:
        I = "Subtract real spectrum (r) or synthetic spectrum (s)"
        log.append([I,input('\n'+I+'\n>>> ')])
        ans = log[-1][1]
    else:
        ans = 'S'
    
    if ans.upper()=='R':
        d.plotoverlay(NHXo, NHYo, invy=False)
        I = "\nPress enter to continue to edited plot"
        log.append([I,input(I)])
        outnh = copy.deepcopy(NHXo)
        if sum((np.array(outnh.grid.shape)-np.array(NHYo.grid.shape))==[0,0])!=0:
            print("\nSpectra do not have the same dimensions and cannot be subtracted")
            continue
        outnh.grid -= NHYo.grid
    else:
        print("\nIntensity |avg| std: %.3f %.3f\nsig: ~1 for ppm, ~0.15 for kHz"%(abs(np.mean(pki)), np.std(pki)))
        I = "Input gaussian amp,sigx,sigy"
        log.append([I,input('\n'+I+'\n>>> ')])
        [a,sx,sy] = [float(m) for m in log[-1][1].split(",")]
        syn2.grid = genspec(syn2, nh2, a, sx, sy)
        d.plotoverlay(NHXo, syn2, invy=False)
        input("\nPress enter to continue to edited plot")
        outnh = copy.deepcopy(NHXo)
        outnh.grid -= syn2.grid
        outnh.set_flce(np.min(outnh.grid), 0)
        outnh.grid /= abs(np.min(outnh.grid))
    
    levels1 = np.linspace(np.min(outnh.grid), -0.1, 50)
    plot(outnh, nh2, invy=False)
    plt.gcf().canvas.draw() # To update, this must be executed outside subroutine that does the plotting
    plt.gcf().canvas.flush_events()
    
    I = "Happy (y/n)?"
    log.append([I, input('\n'+I+'\n>>> ')])
    ans = log[-1][1]
    if ans.upper()=='Y':
        tick=1
    else:
        plt.close()
plt.close('all')

# Set crosspeak energies
aEnn = np.zeros((nn2.shape[0], nn2.shape[0]))
aEnh = np.zeros((nh2.shape[0], nh2.shape[0]))
for i in range(nn.shape[0]):
    for j in range(nn.shape[0]):
        aEnn[i,j] = Expk(outnn, nn2[j], nn2[i])
        aEnh[i,j] += Expk(outnh, nh2[j], nh2[i])
aE = [aEnn, aEnh]

writeinput('Potgrids/input.csv', aE, NNH, res=restyp)
print("\nWrote input.csv to directory Potgrids/")

writespec("specnn.csv", outnn)
print("\nWrote specnn.csv to directory Potgrids/")
writespec("specnh.csv", outnh)
print("\nWrote specnh.csv to directory Potgrids/")

f = open("Order/pkey.csv","w")
f.write(",".join([str(i) for i in nn2[:,0]])+"\n")
f.write(",".join([str(i) for i in nn2[:,1]])+"\n")
f.write(",".join([str(i) for i in nh2[:,0]])+"\n")
f.write(",".join([str(i) for i in nh2[:,1]])+"\n")
f.close()
print("\nWrote pkey.csv to directory Order/")

from datetime import datetime
with open("Output/speced.log","w") as f:
    t = datetime.now()
    f.write("Timestamp: %4d/%02d/%02d | %02d:%02d:%02d\n"%(t.year, t.month, t.day, t.hour, t.minute, t.second))
    for m in log:
        f.write('P: %s\nA: %s\n'%(m[0],m[1]))
print("\nWrote speced.log to directory Output/\n")
