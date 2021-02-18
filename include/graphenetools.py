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
            conventional_lj_parameters = np.array([ 2.74, 16.2463 ])
            optimized_lj_parameters = np.array(
                [[ 2.63023681, 17.44768835],
                 [ 2.66077881, 17.08534452],
                 [ 2.69725342, 16.66725476],
                 [ 2.74310303, 16.1616591 ],
                 [ 2.75562744, 16.01401717],
                 [ 2.7677124 , 15.87153265],
                 [ 2.78001709, 15.73263853],
                 [ 2.79349365, 15.58078085],
                 [ 2.80623779, 15.43983747],
                 [ 2.82220459, 15.26679227],
                 [ 2.8380249 , 15.09702875],
                 [ 2.85479736, 14.92298762],
                 [ 2.87230224, 14.74216293],
                 [ 2.88739014, 14.58550123],
                 [ 2.91156006, 14.35380853],
                 [ 2.93309326, 14.1474397 ],
                 [ 2.95718994, 13.92186853],
                 [ 2.98238525, 13.69287949],
                 [ 3.13267822, 12.54164905]
                ]
            )
            optimized_strain_array = np.array(
                [ 0.00,
                  0.05,
                  0.10,
                  0.15,
                  0.16,
                  0.17,
                  0.18,
                  0.19,
                  0.20,
                  0.21,
                  0.22,
                  0.23,
                  0.24,
                  0.25,
                  0.26,
                  0.27,
                  0.28,
                  0.29,
                  0.34]
            )


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
                    verbose = False,
                    use_conventional = False):
        """Useful tools"""

        self.datadir = datadir
        self.plotdir = plotdir
        self.notebook = notebook
        self.verbose = verbose

        if not os.path.exists(datadir):
            os.makedirs(datadir)

        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

        self.const.use_conventional = use_conventional
        self.const.strain = 0.00 if strain is None else strain
        self.const.a0 = 1.42 if a0 is None else a0
        self.const.poisson = 0.165 if poisson is None else poisson
        self.const.conventional_lj_parameters = np.array([ 2.74, 16.2463 ])
        self.const.optimized_lj_parameters = np.array(
            [[ 2.63023681, 17.44768835],
             [ 2.66077881, 17.08534452],
             [ 2.69725342, 16.66725476],
             [ 2.74310303, 16.1616591 ],
             [ 2.75562744, 16.01401717],
             [ 2.7677124 , 15.87153265],
             [ 2.78001709, 15.73263853],
             [ 2.79349365, 15.58078085],
             [ 2.80623779, 15.43983747],
             [ 2.82220459, 15.26679227],
             [ 2.8380249 , 15.09702875],
             [ 2.85479736, 14.92298762],
             [ 2.87230224, 14.74216293],
             [ 2.88739014, 14.58550123],
             [ 2.91156006, 14.35380853],
             [ 2.93309326, 14.1474397 ],
             [ 2.95718994, 13.92186853],
             [ 2.98238525, 13.69287949],
             [ 3.13267822, 12.54164905]
            ]
        )
        self.const.optimized_strain_array = np.array(
            [ 0.00,
              0.05,
              0.10,
              0.15,
              0.16,
              0.17,
              0.18,
              0.19,
              0.20,
              0.21,
              0.22,
              0.23,
              0.24,
              0.25,
              0.26,
              0.27,
              0.28,
              0.29,
              0.34]
        )

        self.setStrain()


    def setStrain(self,strain = None, a0 = None, poisson = None):
        """Change strain and lattice vectors"""
        self.const.strain = self.const.strain if strain is None else strain
        self.const.a0 = self.const.a0 if a0 is None else a0
        self.const.poisson = self.const.poisson if poisson is None else poisson
        if self.const.use_conventional:
            self.const.lj_sigma, self.const.lj_epsilon = self.const.conventional_lj_parameters
        else:
            if self.const.strain > self.const.optimized_strain_array[-1]:
                print("WARNING: strain larger than {0}, using Lennard-Jones parameters for {0}".format(self.const.optimized_strain_array[-1]))
                self.const.lj_sigma, self.const.lj_epsilon = self.const.optimized_lj_parameters[-1]
            elif self.const.strain in self.const.optimized_strain_array:
                lj_parameters_index = self.const.optimized_strain_array == self.const.strain
                self.const.lj_sigma, self.const.lj_epsilon = self.const.optimized_lj_parameters[lj_parameters_index,:].squeeze()
            else:
                lj_parameters_index_gt = self.const.optimized_strain_array > self.const.strain
                lj_parameters_index_lt = self.const.optimized_strain_array < self.const.strain

                lj_parameters_gt = self.const.optimized_lj_parameters[lj_parameters_index_gt,:][0,:]
                lj_parameters_lt = self.const.optimized_lj_parameters[lj_parameters_index_gt,:][-1,:]

                strain_gt = self.const.optimized_strain_array[lj_parameters_index_gt][0]
                strain_lt = self.const.optimized_strain_array[lj_parameters_index_lt][-1]
                print("WARNING: strain {2} not in table, linear interpolating Lennard-Jones parameters for strain between {0} and {1}".format(strain_lt,strain_gt,self.const.strain))

                slope = (lj_parameters_gt - lj_parameters_lt)/(strain_gt - strain_lt)
                intercept = lj_parameters_gt - lj_parameters_lt*slope
                lj_parameters_linear = slope*self.const.strain + intercept

                self.const.lj_sigma, self.const.lj_epsilon = lj_parameters_linear


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
        self.lattice.g1x = 8.*np.pi*np.sqrt(3.)*(4. + (3.*self.const.strain) - (self.const.strain*self.const.poisson))/(3.*self.const.a0*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))*(4.+(3*self.const.strain)-(self.const.strain*self.const.poisson)));
        self.lattice.g1y = 8.*np.pi*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))/(3.*self.const.a0*(4. + self.const.strain - (3.*self.const.strain*self.const.poisson))*(4.+(3*self.const.strain)-(self.const.strain*self.const.poisson)));
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
        # FIXME this code should be rewritten
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
                        u = -135.0,
                        z = 40.0,
                        W = 29,
                        N = None,
                        M = 8,
                        epsilon = None,
                        sigma = None,
                        primitive = False,
                        internal = 'szalewicz',
                        external = 'graphenelut'):
        """ Happacher: x = {40,42,40}, y = {20,22,24}"""
        N = x * y // 6 if N is None else N
        epsilon = self.const.lj_epsilon if epsilon is None else epsilon
        sigma = self.const.lj_sigma if sigma is None else sigma
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

    def plotLattice(self, x = 4, y = 2,dim=None,):
        """Check the lattice"""
        Lx = x*self.lattice.a1x
        Ly = y*self.lattice.a1y

        fig,ax = plt.subplots()
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
        ax.set_xlabel(r'$x\ \mathrm{[\AA]}$')
        ax.set_ylabel(r'$y\ \mathrm{[\AA]}$')
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect(abs(x1-x0)/abs(y1-y0))
        
        return fig,ax


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

    def V_LJ(self,r,epsilon,sigma):
        t12=(sigma/r)**12
        t6=(sigma/r)**6
        return 4*epsilon*(t12-t6)

    def V(self,z,x=0.0,y=0.0,epsilon=16.24642,sigma=2.74,gnum=4):
        #FIXME uncertain if the is the latest potential equation
        print("FIXME: need to look at this equation")
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

if __name__ == '__main__':
    """Graphene tools."""
    #FIXME need to implement this part
    print("FIXME")
