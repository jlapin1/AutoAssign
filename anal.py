#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MultipleLocator

class inputobject:
    
    def __init__(self, fn):
        with open("Potgrids/"+fn, "r") as f:
            
            # amt
            self.amt = int(f.readline())
            
            # Length
            self.len = int(f.readline())
            
            # Starting permutation
            f.readline()
            
            # Residues
            self.res = f.readline().split(",")
            self.res[-1] = self.res[-1].strip()
            
            # Grid
            self.grid = np.zeros((self.amt, self.len, self.len))
            for i in range(self.amt):
                for j in range(self.len):
                    self.grid[i][j] = [float(k) for k in f.readline().split(",")]
        
        self.readkey()
        self.readspec(["specnn.csv","specnh.csv"])
        
    def readkey(self):
        with open("Order/pkey.csv","r") as f:
            self.key = np.zeros((self.amt, 2, self.len))
            for i in range(self.amt):
                for j in range(2):
                    self.key[i][j] = [float(k) for k in f.readline().split(",")]
    
    def readspec(self, fnames):
        self.spec = []
        self.raxis = []
        self.caxis = []
        self.rsz = np.zeros((self.amt), dtype='int')
        self.csz = np.zeros((self.amt), dtype='int')
        self.rmax = np.zeros((self.amt))
        self.rmin = np.zeros((self.amt))
        self.cmax = np.zeros((self.amt))
        self.cmin = np.zeros((self.amt))
        for i in range(self.amt):
            with open("Potgrids/"+fnames[i], "r") as f:
                [self.rsz[i], self.csz[i]] = [int(m) for m in f.readline().split(",")]
                [self.rmax[i], self.rmin[i], self.cmax[i], self.cmin[i]] = [float(m) for m in f.readline().split(",")]
                self.spec.append( np.zeros((self.rsz[i], self.csz[i])) )
                self.raxis.append(np.linspace(self.rmax[i], self.rmin[i], self.rsz[i]))
                self.caxis.append(np.linspace(self.cmax[i], self.cmin[i], self.csz[i]))
                for m in range(self.rsz[i]):
                    self.spec[i][m] = [float(n) for n in f.readline().split(",")] 

class outputobject:
    
    def __init__(self,fn,io):
        self.io = io
        with open("Output/"+fn, "r") as f:
            self.top = int(f.readline())
            self.perms = []
            for line in f:
                self.perms.append([int(m) for m in line.split(",")])
        self.perms = np.array(self.perms)
        
        self.pks = np.zeros((self.top, io.amt, 2, io.len))
        for i in range(self.top):
            for j in range(io.amt):
                for k in range(2):
                    self.pks[i][j][k] = io.key[j][k][self.perms[i]]
        
        self.xpeaks()
        self.eng3()
    
    def xpeaks(self):
        length = self.pks.shape[-1]
        self.xpks = np.zeros(( self.top, self.io.amt, 2, 2*(length-1) ))
        for i in range(self.top):
            for j in range(self.io.amt):
                 for k in range(length-1):
                    self.xpks[i][j][0][k] = self.pks[i][j][0][k]
                    self.xpks[i][j][1][k] = self.pks[i][j][1][k+1]
                    self.xpks[i][j][0][k+length-1] = self.pks[i][j][0][k+1]
                    self.xpks[i][j][1][k+length-1] = self.pks[i][j][1][k]
    
    def eng3(self):
        xlength = self.xpks[0][0].shape[1]
        self.rE = np.zeros(( self.top, self.io.amt, xlength ))
        self.E = np.zeros(( self.top, self.io.amt+1 ))
        for i in range(self.top):
            for j in range(self.io.amt):
                rspac = self.io.raxis[j][0]-self.io.raxis[j][1]
                cspac = self.io.caxis[j][0]-self.io.caxis[j][1]
                for k in range(xlength):
                    rind = round((self.io.rmax[j]-self.xpks[i][j][1,k])/rspac)
                    cind = round((self.io.cmax[j]-self.xpks[i][j][0,k])/cspac)
                    self.rE[i][j][k] = self.io.spec[j][int(rind), int(cind)]
                    self.E[i][j] += self.rE[i][j][k]
                self.E[i][-1] += self.E[i][j]
    
    def plot2(self, ind=-1, start=0, block=True):
        
        # Plotting specs
        fig = plt.figure(figsize=(8, 12))
        # Create a grid of subplots (2 rows, 1 column)
        wid = gs.GridSpec(2,1,width_ratios=[1])
        
        ax1 = plt.subplot(wid[0,:]) 
        ax2 = plt.subplot(wid[1,:])
    
        levels1 = np.linspace(self.io.spec[0].min(), -0.1, 50)
        levels2 = np.linspace(self.io.spec[1].min(), -0.1, 50)
        
        raxis1 = self.io.raxis[0]
        caxis1 = self.io.caxis[0]
        C1, R1 = np.meshgrid(caxis1, raxis1)  # X/Y have the x/y-coords for every pt, separately
        raxis2 = self.io.raxis[1]
        caxis2 = self.io.caxis[1]
        C2, R2 = np.meshgrid(caxis2, raxis2)
        
        cnt1 = ax1.contour(C1, R1, self.io.spec[0], levels1, linewidths=0.5, cmap='Blues')
        cnt2 = ax2.contour(C2, R2, self.io.spec[1], levels2, linewidths=0.5, cmap='Blues')
    
        ax1.plot(self.pks[0][0][0], self.pks[0][0][1], 'o', markersize=5, markerfacecolor='none', markeredgecolor='k')
        ax2.plot(self.pks[0][1][0], self.pks[0][1][1], 'o', markersize=5, markerfacecolor='none', markeredgecolor='k')
        
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax2.invert_xaxis()
    
        wid.update(hspace=0) # upper plot on top of lower plot
        
        ax1.xaxis.set_major_locator(MultipleLocator(5));ax1.xaxis.set_minor_locator(MultipleLocator(1))
        ax1.yaxis.set_major_locator(MultipleLocator(5));ax1.yaxis.set_minor_locator(MultipleLocator(1))
        ax1.tick_params(axis='x', which='both', direction='in')
        
        ax2.xaxis.set_major_locator(MultipleLocator(5));ax2.xaxis.set_minor_locator(MultipleLocator(1))
        ax2.yaxis.set_major_locator(MultipleLocator(1));ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        ax1.set_ylabel('$^{15}N$ Chemical Shift (ppm)')
        ax2.set_xlabel('$^{15}N$ Chemical Shift (ppm)')
        ax2.set_ylabel('$^{1}H-^{15}N$ Dipolar Coupling (kHz)')
        
        if ind>=0:
            ax1.plot(self.xpks[ind][0][0], self.xpks[ind][0][1], 'rx')
            ax2.plot(self.xpks[ind][1][0], self.xpks[ind][1][1], 'rx')
            for i in range(self.io.len):
                ax1.text(self.pks[ind][0][0][i], self.pks[ind][0][1][i], '%c%d'%(self.io.res[i],start+i))
                ax2.text(self.pks[ind][1][0][i], self.pks[ind][1][1][i], '%c%d'%(self.io.res[i],start+i))
            
            ax1.set_title('E = %f (%f, %f)'%(self.E[ind][-1],self.E[ind][0],self.E[ind][1]))
        
        if block==True:
            print("\nClose plot to continue")
        plt.show(block=block)

class filterer():
    def __init__(self, out):
        self.out = out
        self.inds = np.arange(self.out.top)
    
    def selectxpk(self, spec, bounds, typ, num):
        xlength = 2*(self.out.pks.shape[-1]-1)
        inds_out = []
        for i in range(len(self.inds)):
            tick=0
            for j in range(xlength):
                if (self.out.xpks[self.inds[i]][spec][0][j]>bounds[0])&\
                   (self.out.xpks[self.inds[i]][spec][0][j]<bounds[1])&\
                   (self.out.xpks[self.inds[i]][spec][1][j]>bounds[2])&\
                   (self.out.xpks[self.inds[i]][spec][1][j]<bounds[3]):
                        tick+=1
            if typ==0: # at least
                if tick>=num:
                    inds_out.append(self.inds[i])
            elif typ==1: # at most
                if tick<=num:
                    inds_out.append(self.inds[i])
            elif typ==2:
                if tick==num: # exactly
                    inds_out.append(self.inds[i])
            else: # range
                if (len(num)!=2):
                    print('"num" must be a list of 2 numbers')
                if (tick>=num[0]) & (tick<=num[1]):
                    inds_out.append(self.inds[i])
        return inds_out
    
    def orfilt(self, spec, cut, tol):
        return [self.inds[i] for i,x in enumerate(np.sum(self.out.rE[self.inds,spec]>cut, axis=1)<=tol) if x]
    
    def process(self):
        svspec = []
        svbound = []
        svtyp = []
        svnum = []
        svlen = []
        tick = 0
        while tick==0:
            
            ans = int(input("\nPlot # (0-%d; <0 for none)\n>>> "%(len(self.inds)-1)))
            if ans>=0:
                plt.close('all')
                self.out.plot2(ind=self.inds[ans], block=False)
            
            ans = int(input("\nFilter type\n0: Orphan xpeaks\n1: Select xpeaks\n>>> "))
            
            spec = int(input("\nWhich spectrum (#: 0 for top, 1 for bottom)?\n>>> "))
            if ans==0:
                cut = float(input("\nCutoff, above which a violation is tallied (-1-0)\n>>> "))
                
                v = np.sum(self.out.rE[self.inds, spec]>cut, axis=1)
                print("\nViolations Min/Median/Max: %d/%d/%d"%(np.min(v),np.median(v),np.max(v)))
                tol = int(input("\nMaximum # of violations (>=0)\n>>> "))
                ilist = self.orfilt(spec, cut, tol)
            elif ans==1:
                I = "Bounds: xmin,xmax,ymin,ymax or cw (current window)"
                log = input("\n"+I+"\n>>> ")
                ax = plt.gcf().get_axes()
                log = (list(list(np.sort(ax[spec].get_xlim()))+list(np.sort(ax[spec].get_ylim()))) if log=='cw' else log.split(","))
                bounds = [float(m) for m in log]
                typ = int(input("\n# of xpeaks\n0: at least\n1: at most\n2: exactly\n<0: range\n>>> "))
                if typ<0:
                    num = [int(m) for m in input("\n# of peaks: min,max\n>>> ").split(",")]
                else:
                    num = int(input("\n# of xpeaks (#)\n>>> "))
                
                ilist = self.selectxpk(spec, bounds, typ, num)
                svspec.append(spec);svbound.append(bounds);svtyp.append(typ);svnum.append(num);svlen.append(len(ilist))
            
            if len(ilist)==0:
                print("\nNo solutions. Reverting back to previous indices.")
                if ans==1:
                    del(svspec[-1],svbound[-1],svtyp[-1],svnum[-1],svlen[-1])
            else:
                self.inds = ilist
            print("\n%d solutions"%(len(self.inds)))
            ans = input("\nKeep filtering (y) or quit process (n)?\n>>> ")
            if ans.upper()=='Y':
                continue
            else:    
                print("\nTop solutions remaining:")
                print(self.inds)
                tick=1
            
            
def plotprof():
    with open("Output/E_profile.csv","r") as f:
        no_incr = int(f.readline())
        Esaveall = np.zeros((no_incr+1, 4))
        for i in range(no_incr+1):
            Esaveall[i] = [float(j) for j in f.readline().split(',')]
    incr = np.arange(no_incr+1)
    fig = plt.figure()
    ax = plt.subplot(111)
    ax2 = ax.twinx()
    ax2.plot(incr,Esaveall[:,0],'k',linewidth=0.5,label='W (~ $T^{-1}$)')
    ax.plot(incr,Esaveall[:,1],'r',linewidth=0.4,label='$E_{tot}$')
    ax.plot(incr,Esaveall[:,2],'b',linewidth=0.4,label='$E_{hom}$')
    ax.plot(incr,Esaveall[:,3],'g',linewidth=0.4,label='$E_{het}$')
    ax.set_xlabel('increment')
    ax2.set_ylabel('W (~ $T^{-1}$)')
    ax.set_ylabel('E (arb)')
    ax2.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=1,
               ncol=3, borderaxespad=0.)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, borderaxespad=0.)
    #fig.savefig('Pictures/E_profile.pdf', format='pdf', dpi=100)
    plt.show()

io = inputobject("input.csv")
out = outputobject("output.csv", io)

f = filterer(out)

tick=0
while tick==0:
    ans = int(input("\nChoose an option\n0: Plot solution\n1: Filter solutions\n2: Exit\n>>> "))
    if ans==0:
        ans2 = int(input("\nSolution index (0-%d)\n>>> "%(out.top-1)))
        out.plot2(ind=ans2, block=True)
    elif ans==1:
        f = filterer(out)
        f.process()
    else:
        tick=1
        print()
