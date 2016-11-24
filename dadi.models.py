#! ~/py33/bin
# -*- coding: utf-8 -*-

"""modeldemo module, contains 7 different divergence models, can be modified to include more models"""

import numpy
import dadi


def SI(params, (n1,n2), pts):
    nu1, nu2, Ts = params
    """
    Strict Isolation model, no gene flow during divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts = params
    """
    Isolation-with-Migration model.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def AM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts = params
    """
    Ancient Migration model.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def SC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc = params
    """
    Secondary Contact model.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between secondary contact and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def IM2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, P = params
    """
    Isolation-with-Migration model with two categories of loci experiencing different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P
    fs = P*fsN+(1-P)*fsI
    return fs

def AM2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Tam, Ts, P = params
    """
    Ancient Migration model with two categories of loci experiencing different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P
    fs = P*fsN+(1-P)*fsI
    return fs

def SC2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P = params
    """
    Secondary Contact model with two categories of loci experiencing different migration rates. 

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P
    fs = P*fsN+(1-P)*fsI
    return fs

def DB((nu1B, nu1F, nu2B, nu2F, nuA, m12, m21, Tp, T), (n1,n2), pts):
    """
    Model with growth, split, bottleneck in pop2 , exp recovery, asymetrical migration
    nuA : The ancestral population size after growth. (Initial size is 1.)
    nu1B: The bottleneck size for pop1
    nu1F: The final size for pop1 (The ancestral population size after growth is nu1F)
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m12: migration from pop 1 to pop2
    m21: migration from pop 2 to pop 1
    Tp: The scaled time between ancestral population growth and the split.
    T: The time between the split and present
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We need to define a function to describe the non-constant population
    # size. lambda is a convenient way to do so.
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/T)
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

