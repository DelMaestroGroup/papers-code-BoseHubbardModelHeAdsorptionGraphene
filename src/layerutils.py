''' layerutils.py

Helper utililities for analzing the layering of 4He on graphene.
Author: Adrian Del Maestro
Date: 2020-07-11

'''

from fractions import Fraction

# get the graphenetools
import sys
sys.path.append("..")
import graphenetools as gt
import numpy as np
import matplotlib.pyplot as plt
import os

# ----------------------------------------------------------------------
def lab(T=None,Lz=None,N=None,n=None,β=None):
    '''A common simulation label based on named parameters.'''
    
    lvars = dict(vars().items())
    label = ''
    unit = {'T':' K', 'Lz':' Å', 'N':'','n':'','β':' K⁻¹'}
    lformat = {'T':'.2f', 'Lz':'.2f', 'N':'03d','n':'','β':'.4f'}
    
    # turn the density into a filling fraction
    if n:
        lvars['n'] = Fraction(n).limit_denominator()
        
    # create the descriptive key
    for key,val in lvars.items():
        if val is not None:
            label += f'{key} = {val:{lformat[key]}}{unit[key]}, ' 
            
    label = label.rstrip(', ')
    return label

# ----------------------------------------------------------------------
def vals(sim):
    'Return dictionary of values based on a simulation label.'
    lvals = {}
    lconvert = {'Lz':float, 'N':int, 'T':float, 'n':Fraction, 'β':float}
    
    for kvp in sim.split(','):
        key,val = kvp.split('=')
        key = key.strip()
    
        # remove any units from the string
        for unit in ['Å','K']:
            val = val.replace(unit,'')
        lvals[key] = lconvert[key](val)
    
    return lvals
    
# ----------------------------------------------------------------------
def texformat(sim,split_lines=False):
    '''Wrap a key label in some tex for pretty printing.'''
    
    labs = [f"${l.strip()}$" for l in sim.split(',')]
    labs = [l.replace("Å",r"\,\mathrm{\AA}") for l in labs]
    labs = [l.replace("Lz",r"L_z") for l in labs]
    labs = [l.replace("K",r"\,\mathrm{K}") for l in labs]
    labe = [l.replace("⁻¹",r"^{-1}") for l in labs]
    
    if split_lines:
        return "\n".join(labs)
    else:
        return ", ".join(labs)
    
# ----------------------------------------------------------------------
def get_graphene_lattice(L):
    graphene = gt(notebook=True, strain=0.0)
    aₒ = 1.42
    a1 = np.array([graphene.lattice.basis.b1x,graphene.lattice.basis.b1y])
    a2 = np.array([graphene.lattice.basis.b2x,graphene.lattice.basis.b2y])
    lattice_x,lattice_y = graphene.hexagonal_lattice()

    # Determine the number of hexagonal plaquettes
    num_x = int(round(L[0]/graphene.lattice.a1x))
    num_y = int(round(L[1]/graphene.lattice.a1y))

    #number of adsorption sites = (x/2)*y
    print('number of adsorption sites:\t{}'.format(num_x*num_y/2))
    #N_sites = int(num_x*num_y/2)
    print('number in C1/3 phase:\t\t{}'.format(num_x*num_y/2/3))
    print(f'Lattice Dimensions: {L[0]:.3f} Å x {L[1]:.3f} Å')

    # This is a subset of the lattice
    lattice_x_idx = np.where(np.abs(lattice_x)<1.5*L[0]/2)[0]
    lattice_y_idx = np.where(np.abs(lattice_y)<1.5*L[1]/2)[0]
    lattice_idx = np.intersect1d(lattice_x_idx,lattice_y_idx)

    fig,ax = plt.subplots(figsize=(4,4))

    ax.hlines(-L[1]/2,-L[0]/2.,L[0]/2., color='r', ls=':')
    ax.hlines(L[1]/2.,-L[0]/2.,L[0]/2., color='r', ls=':')
    ax.vlines(-L[0]/2.,-L[1]/2.,L[1]/2., color='r', ls=':')
    ax.vlines(L[0]/2.,-L[1]/2.,L[1]/2., color='r', ls=':')

    plt.scatter(lattice_x[lattice_idx],lattice_y[lattice_idx], s=5, c='k',zorder=10)

    plt.axis('equal')
    plt.xlim(-1.2*L[0]/2,1.2*L[0]/2)
    plt.ylim(-1.2*L[1]/2,1.2*L[1]/2)
    plt.xlabel('x (Å)')
    plt.ylabel('y (Å)');
    
    return aₒ,a1,a2,lattice_x,lattice_y,lattice_idx,fig

# ----------------------------------------------------------------------
def get_base_dir(N,T=None):
    base_dir = f'{os.environ["HeGrapheneData"]}/N_eq_{N:03d}/'
    if T is not None and T <= 0.0:
        base_dir += 'PIGS/'
    return base_dir