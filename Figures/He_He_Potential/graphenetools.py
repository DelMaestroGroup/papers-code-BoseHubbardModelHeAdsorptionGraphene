#!/usr/bin/env python3
# Filename: graphenetools.py
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from scipy.special import kn
class GrapheneTools:
    """Contains useful functions for working with graphene"""

    class const:
        """Hold constants"""
        def __init__(self):
            pass

    class lattice:
        """Hold lattice vectors"""
        def __init__(self):
            pass

        class basis:
            """Hold basis vectors"""
            def __init__(self):
                pass

    def __init__(   self, 
                    strain = None,
                    a0 = None,
                    poisson = None,
                    datadir = './data/',
                    plotdir = './plots/',
                    notebook = False,
                    verbose = False):
        """Useful tools"""

        self.datadir = datadir
        self.plotdir = plotdir
        self.notebook = notebook
        self.verbose = verbose

        if not os.path.exists(datadir):
            os.makedirs(datadir)

        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

        self.const.strain = 0.00 if strain is None else strain
        self.const.a0 = 1.42 if a0 is None else a0
        self.const.poisson = .165 if poisson is None else poisson

        self.setStrain()


    def setStrain(self,strain = None, a0 = None, poisson = None):
        """Change strain and lattice vectors"""
        self.const.strain = self.const.strain if strain is None else strain
        self.const.a0 = self.const.a0 if a0 is None else a0
        self.const.poisson = self.const.poisson if poisson is None else poisson

        #Lattice vectors
        self.lattice.a1x = (np.sqrt(3.)*self.const.a0/8.)*(4.+self.const.strain-(3.*self.const.strain*self.const.poisson))
        self.lattice.a1y = (3.*self.const.a0/8.)*(4.+(3.*self.const.strain)-(self.const.strain*self.const.poisson))
        self.lattice.a2x = -self.lattice.a1x
        self.lattice.a2y = self.lattice.a1y

        #basis vectors
        self.lattice.basis.b1x = self.lattice.a2x
        self.lattice.basis.b1y = -0.5 * (1. + self.const.strain) * self.const.a0
        self.lattice.basis.b2x = 0.
        self.lattice.basis.b2y = -.5*(1. + self.const.strain) * self.const.a0-(self.lattice.a1x/np.sqrt(3.))

        # reciprocal lattice vectors
        self.lattice.g1x = (8.*np.pi*np.sqrt(3.)*(4. + (3.*self.const.strain) - (self.const.strain*self.const.poisson))
							/(3.*self.const.a0
							*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))
							*(4.+(3*self.const.strain)-(self.const.strain*self.const.poisson))));
        self.lattice.g1y = (8.*np.pi*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))
							/(3.*self.const.a0
							*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))
							*(4.+(3*self.const.strain)-(self.const.strain*self.const.poisson))));
        self.lattice.g2x = -self.lattice.g1x;
        self.lattice.g2y = self.lattice.g1y;
        
        # area of unit cell
        self.lattice.A = np.abs((self.lattice.a1x*self.lattice.a2y) - (self.lattice.a1y*self.lattice.a2x));        

    def hexagonal_lattice(  self,
                            numx = 256,
                            numy = 256,
                            strain = None,
                            a0 = None,
                            poisson = None):
        """Returns x and y coordinates for a hexagonal lattice."""
        strain = self.const.strain if strain is None else strain
        a0 = self.const.a0 if a0 is None else a0
        poisson = self.const.poisson if poisson is None else poisson
        ax = a0 * (1.00 + (strain * (1.0 - (3.0 * poisson)) / 4.0))
        ay = a0 * (1.00 + strain)
        d0 = np.sqrt(3.0) * ax
        superxpoints=np.arange(0,numx,1)
        superypoints=np.arange(0,numy,1)
        transx = (superxpoints-numx/2)*d0
        transy = (ax + ay +ay)*(superypoints-numy/2)
        supertestx=np.hstack((transx[:,np.newaxis],transx[:,np.newaxis],d0/2+transx[:,np.newaxis],d0/2+transx[:,np.newaxis])*numy)
        supertesty=np.vstack((np.hstack(((((ay+ax)/2) + transy)[:,np.newaxis],((-(ay+ax)/2) + transy)[:,np.newaxis],((ay/2) + transy)[:,np.newaxis],((-ay/2) + transy)[:,np.newaxis])).tolist())*numx)
        return np.array(supertestx).T.ravel('F'), np.array(supertesty).T.ravel('F')

    def pimcCommand(    self,
                        x = 4,
                        y = 2,
                        T = 0.5,
                        t = 0.003,
                        E = 100000,
                        S = 1000000,
                        u = -135.,
                        z = 40.,
                        W = 29,
                        N = None,
                        M = 8,
                        epsilon = 16.247,
                        sigma = 2.739,
                        primitive = False,
                        internal = 'szalewicz',
                        external = 'graphenelut'):
        """ Happacher: x = {40,42,40}, y = {20,22,24}"""
        N = x * y / 6 if N is None else N
        Lx = x*self.lattice.a1x
        Ly = y*self.lattice.a1y
        strain = self.const.strain
        if primitive:
                    return 'pimc.e -T %f -N %d -W %d -t %f -M %d -C 1.0 -I %s -X %s --lj_sigma %f --lj_epsilon %f -E %d -S %d -l 7 -u %f --relax --action=primitive --Lx %f --Ly %f --Lz %f --strain %f' % (T,N,W,t,M,internal,external,sigma,epsilon,E,S,u,Lx,Ly,z, strain)

        return 'pimc.e -T %f -N %d -W %d -t %f -M %d -C 1.0 -I %s -X %s --lj_sigma %f --lj_epsilon %f -E %d -S %d -l 7 -u %f --relax --Lx %f --Ly %f --Lz %f --strain %f' % (T,N,W,t,M,internal,external,sigma,epsilon,E,S,u,Lx,Ly,z, strain)


    def pimcCommandBatch(   self,
                            umin=-135,
                            umax=-25,
                            numu=20, **kwargs):
        """Make a batch file for pimc"""
        uu = np.linspace(umin,umax,numu)

        fout = open(self.datadir + 'graphenerun-%d.sh' % int(time.time()*1000),'w')
        #fout.write('#!/bin/bash\n')
        for i, u in enumerate(uu):
            fout.write(self.pimcCommand(u=u,**kwargs) + '\n')
        fout.close()

    def pimcCommandBatch2(   self,
                            tmin=0.0025,
                            tmax=0.0200,
                            num=250,
                            dups=1,
                            timeStart = None,
                            **kwargs):
        """Make a batch file for pimc"""#FIXME need to make just one batch command
        timeStart = int(time.time()*1000) if timeStart is None else timeStart
        tt = np.linspace(tmin,tmax,num)
        fout = open(self.datadir + 'graphenerun-%d.sh' % int(time.time()*1000),'w')
        #fout.write('#!/bin/bash\n')
        for i, t in enumerate(tt):
            for j in range(dups):
                dTime = int(time.time()*1000) - timeStart 
                fout.write(self.pimcCommand(t=t,**kwargs) + ' -p {}\n'.format(10000*dTime))
                time.sleep(1e-3)
        fout.close()
        
    def pimcCommandBatch2a(   self,
                            tmin=0.0025,
                            tmax=0.0200,
                            num=250,
                            dups=1,
                            timeStart = None,
                            **kwargs):
        """Make a batch file for pimc"""#FIXME need to make just one batch command
        timeStart = int(time.time()*1000) if timeStart is None else timeStart
        tt = np.linspace(tmin,tmax,num)
        fout = open(self.datadir + 'graphenerun-%d.sh' % int(time.time()*1000),'w')
        #fout.write('#!/bin/bash\n')
        for i, t in enumerate(tt):
            for j in range(dups):
                dTime = int(time.time()*1000) - timeStart 
                fout.write(self.pimcCommand(t=t,**kwargs) + ' --canonical --gaussian_window_width 0.5 --window 1 -p {}\n'.format(10000*dTime))
                time.sleep(1e-3)
        fout.close()
        
    def pimcCommandBatch2b(   self,
                            tmin=0.0025,
                            tmax=0.0200,
                            num=250,
                            dups=1,
                            timeStart = None,
                            **kwargs):
        """Make a batch file for pimc"""#FIXME need to make just one batch command
        timeStart = int(time.time()*1000) if timeStart is None else timeStart
        tt = np.linspace(tmin,tmax,num)
        fout = open(self.datadir + 'graphenerun-%d.sh' % int(time.time()*1000),'w')
        #fout.write('#!/bin/bash\n')
        for i, t in enumerate(tt[0:3]):
            for j in range(dups):
                dTime = int(time.time()*1000) - timeStart 
                fout.write(self.pimcCommand(t=t,**kwargs) + ' --canonical --gaussian_window_width 0.5 --window 1 -p {}\n'.format(10000*dTime))
                time.sleep(1e-3)
        fout.close()

    def pimcCommandBatch3(  self,
                            dups=125,
                            timeStart = None,
                            **kwargs):
        """Make a batch file for pimc"""#FIXME need to make just one batch command
        strAdded = ' -e "superfluid fraction" -e "energy" -e "number particles" -e "diagonal fraction"'
        timeStart = int(time.time()*1000) if timeStart is None else timeStart
        fout = open(self.datadir + 'graphenerun-%d.sh' % int(time.time()*1000),'w')
        #fout.write('#!/bin/bash\n')
        for j in range(dups):
            dTime = int(time.time()*1000) - timeStart 
            fout.write(self.pimcCommand(**kwargs) + strAdded + ' -p {}\n'.format(10000*dTime))
            time.sleep(1e-3)
        fout.close()
    
    def plotLattice(self, x = 4, y = 2,dim=None,):
        """Check the lattice"""
        Lx = x*self.lattice.a1x
        Ly = y*self.lattice.a1y

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xpoints, ypoints = self.hexagonal_lattice(x,y)
        ax.plot(xpoints, ypoints,'ro')

        ax.set_xlim((-Lx/2.,Lx/2.))
        ax.set_ylim((-Ly/2.,Ly/2.))
        if dim:
            ax.set_xlim((dim[0],dim[1]))
            ax.set_ylim((dim[2],dim[3]))
            ax.hlines(-Ly/2.,-Lx/2.,Lx/2.)
            ax.hlines(Ly/2.,-Lx/2.,Lx/2.)
            ax.vlines(-Lx/2.,-Ly/2.,Ly/2.)
            ax.vlines(Lx/2.,-Ly/2.,Ly/2.)
        ax.set_xlabel(r'$x~[\si{\angstrom}]$')
        ax.set_ylabel(r'$y~[\si{\angstrom}]$')
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect(abs(x1-x0)/abs(y1-y0))
        
        fig.savefig(self.plotdir + 'latticeCheck.pdf',bbox_inches='tight')
        if not self.notebook:
            plt.close(fig)

    def saveLattice(self, x = 4, y = 2, z=None, dim=None,fname='carbonPositions.txt'):
        """Save the lattice"""
        Lx = x*self.lattice.a1x
        Ly = y*self.lattice.a1y

        xpoints, ypoints = self.hexagonal_lattice(x,y)
        

        if dim:
            ind = np.where(xpoints >= dim[0])
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]
            
            ind = np.where(xpoints <= dim[1])
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]
            
            ind = np.where(ypoints >= dim[2])
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]
            
            ind = np.where(ypoints <= dim[3])
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]

        else:
            ind = np.where(np.abs(xpoints) <= Lx/2.)
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]

            ind = np.where(np.abs(ypoints) <= Ly/2.)
            xpoints = xpoints[ind]
            ypoints = ypoints[ind]
        
        if z:
            r = np.squeeze(np.dstack((xpoints.T,ypoints.T,z*np.ones_like(xpoints).T)))
        else:
            r = np.squeeze(np.dstack((xpoints.T,ypoints.T)))
        np.savetxt(fname,r)

    
    def V_0(self,z,A,epsilon,sigma):
        v = (4.*np.pi/A)*epsilon*sigma*sigma*( ((2./5.)*((sigma/z)**10)) - ((sigma/z)**4) )
        return v

    def V_g(self,z,A,epsilon,sigma,g):
        k5term = ((g*sigma*sigma/2./z)**5.)*kn(5,g*z)/30.
        k2term = 2.*((g*sigma*sigma/2./z)**2.)*kn(2,g*z)
        prefactor = epsilon*sigma*sigma*2.*np.pi/A

        v = prefactor*(k5term-k2term)
        return v

    def C_g(self,z,sigma, m, n):
        g = np.sqrt(((m*self.lattice.g1x + n*self.lattice.g2x)**2) + ((m*self.lattice.g1y + n*self.lattice.g2y)**2));
        print(g,"\n")
        k5term = ((g*sigma*sigma/2./z)**5.)*kn(5,g*z)/30.
        k2term = 2.*((g*sigma*sigma/2./z)**2.)*kn(2,g*z)
        #prefactor = epsilon*sigma*sigma*2.*np.pi/A

        v = (k5term-k2term)
        return v

    def V_LJ(self,r,epsilon,sigma):
        t12=(sigma/r)**12
        t6=(sigma/r)**6
        return 4*epsilon*(t12-t6)

    def V(self,z,x=0.0,y=0.0,epsilon=16.24642,sigma=2.74,gnum=2):
        v = self.V_0(z,self.lattice.A,epsilon,sigma)

        for j,zz in enumerate(z):
            for m in np.arange(-gnum,gnum+1):
                for n in np.arange(-gnum,gnum+1):
                    if ((m != 0) or (n !=0)):
                        g = np.sqrt(((m*self.lattice.g1x + n*self.lattice.g2x)**2) + ((m*self.lattice.g1y + n*self.lattice.g2y)**2));

                        gdotb1 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b1x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b1y+y));
                        gdotb2 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b2x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b2y+y));

                        v_g = self.V_g(zz,self.lattice.A,epsilon,sigma,g)                    

                        v[j] += (np.cos(gdotb1)+np.cos(gdotb2))*v_g
        return v

    def V_test(self,z,k1,k2,x=0.0,y=0.0,epsilon=16.24642,sigma=2.74,gnum=2):
        v = 0*self.V_0(z,self.lattice.A,epsilon,sigma)

        for j,zz in enumerate(z):
            for m in [k1]:#np.arange(-gnum,gnum+1):
                for n in [k2]:#np.arange(-gnum,gnum+1):
                    if ((m != 0) or (n !=0)):
                        g = np.sqrt(((m*self.lattice.g1x + n*self.lattice.g2x)**2) + ((m*self.lattice.g1y + n*self.lattice.g2y)**2));

                        gdotb1 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b1x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b1y+y));
                        gdotb2 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b2x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b2y+y));

                        v_g = self.V_g(zz,self.lattice.A,epsilon,sigma,g)                    

                        v[j] += (np.cos(gdotb1)+np.cos(gdotb2))*v_g
        return v
		
    def V2(self,z,mm, nn, x=0.0,y=0.0,epsilon=16.24642,sigma=2.74,gnum=2):
        v = 0*self.V_0(z,self.lattice.A,epsilon,sigma)

        for j,zz in enumerate(z):
            for m in np.arange(-gnum,gnum+1):
                for n in np.arange(-gnum,gnum+1):
                    if ((m == mm) and (n == nn)):
                        g = np.sqrt(((m*self.lattice.g1x + n*self.lattice.g2x)**2) + ((m*self.lattice.g1y + n*self.lattice.g2y)**2));

                        gdotb1 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b1x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b1y+y));
                        gdotb2 = ((m*self.lattice.g1x + n*self.lattice.g2x)*(self.lattice.basis.b2x+x)) + ((m*self.lattice.g1y + n*self.lattice.g2y)*(self.lattice.basis.b2y+y));

                        v_g = self.V_g(zz,self.lattice.A,epsilon,sigma,g)                    

                        v[j] += (np.cos(gdotb1)+np.cos(gdotb2))*v_g
                        expigb = (np.cos(gdotb1)+np.cos(gdotb2))
        return v
    
    def junk(self):
        """Needs fixing"""

        xpoints = []
        ypoints = []
        for i in range(10):
            for j in range(10):
                x1 = a1x*i + a2x*j + b1x
                x2 = a1x*i + a2x*j + b2x
                y1 = a1y*i + a2y*j + b1y
                y2 = a1y*i + a2y*j + b2y
                xpoints += [x1]
                xpoints += [x2]
                ypoints += [y1]
                ypoints += [y2]

        xpoints = np.array(xpoints)
        ypoints = np.array(ypoints)
        plt.figure()
        plt.plot(xpoints,ypoints,'ks')
        plt.ylim((0,10))
        plt.xlim((-10,10))
        boxdimx = np.array([49.19, 51.65, 49.19])
        boxdimy = np.array([42.60, 46.86, 51.12])
        boxdimx = np.array([40., 42, 40., 24., 10., 8., 4.])*a1x
        boxdimy = np.array([20., 22., 24., 12., 6., 4., 2.])*a1y

        '''
        b2y = (1. + STRAIN) * a0

        xpoints = []
        ypoints = []
        for i in range(10):
            for j in range(10):
                x1 = a1x*i + a2x*j + b1x
                x2 = a1x*i + a2x*j + b2x
                y1 = a1y*i + a2y*j + b1y
                y2 = a1y*i + a2y*j + b2y
                xpoints += [x1]
                xpoints += [x2]
                ypoints += [y1]
                ypoints += [y2]

        xpoints = np.array(xpoints)
        ypoints = np.array(ypoints)
        '''

        rx = 0. + (2. * a1x) + (3. * a2x) + (b1x)
        ry = 0. + (2. * a1y) + (3. * a2y) + (b1y)
        plt.plot(rx,ry,'bo')
        plt.savefig(self.plotdir + 'junk.pdf')

if __name__ == '__main__':
    """Graphene tools."""
