# -*- coding: utf-8 -*-
"""
pw_stix.py

Solver for cold plasma waves using Stix parameter method

Yasuhito Narita, Uwe Motschmann, and Tohru Hada
y.narita@tu-braunschweig.de
March 2026
License: MIT License

Usage:
    Refer to the README for detailed instructions.
    If you use this code, please cite it as follows:
    DOI: TBD.
"""

# MIT License (Short form)
# Copyright (c) 2026 Y. Narita, U. Motschmann, and T. Hada
# This software is released under the MIT License.
# http://opensource.org


import numpy as np
import matplotlib.pyplot as plt

def dispcold():

    #----------------------------#
    # ion-to-electron mass ratio #
    #----------------------------#
    mass_ratio = 1836.15

    #-------------------------------------------------------------------#
    # electron cyclotron frequency to electron plasma frequency wce/wpe #
    #-------------------------------------------------------------------#
    wce = 0.01

    #----------------------------------------#
    # propagation angle (degrees in integer) #
    #----------------------------------------#
    th = 00

    #----------------------------#
    # output file (numpy binary) #
    #----------------------------#
    outfile = 'stix_elc_mass_1836_00deg'
    pdffile = outfile + '.pdf'

    #-----------------------------------#
    # frequency array (in units of wpe) #
    #-----------------------------------#
    wmin = 0.000000001
    wmax = wce + wmin
    nomega = 20000

    #---------------------------------#
    # avoid singularity at 90 degrees #
    #---------------------------------#
    if th == 90:
        th = 89.9999

    #------------------------#
    # vectorized omega array #
    #------------------------#
    om = np.linspace(wmin, wmax, nomega)
    rm = 1 / mass_ratio
    rd = np.deg2rad(th)
    sin2, cos2 = np.sin(rd)**2, np.cos(rd)**2

    #----------------------------#
    # vectorized Stix parameters #
    #----------------------------#
    rr = 1. - rm / (om * (om + rm * wce)) - 1 / (om * (om - wce))
    ll = 1. - rm / (om * (om - rm * wce)) - 1 / (om * (om + wce))
    pp = 1. - (rm + 1) / (om**2)
    ss = 0.5 * (rr + ll)
    dd = 0.5 * (rr - ll)

    #--------------------------------------#
    # coefficients for dispersion relation #
    #--------------------------------------#
    aa = ss * sin2 + pp * cos2
    bb = (ss**2 - dd**2) * sin2 + pp * ss * (1 + cos2)
    cc = pp * rr * ll
    ff = bb**2 - 4 * aa * cc

    #-----------------------------------#
    # filter for valid physics (ff > 0) #
    #-----------------------------------#
    mask_ff = ff > 0
    results = []

    #----------------------------------------------#
    # solve for both roots (+ and -) using masking #
    #----------------------------------------------#
    for sign in [1, -1]:
        n_sq = (bb[mask_ff] + sign * np.sqrt(ff[mask_ff])) / (2 * aa[mask_ff])

        #----------------------------------#
        # valid refraction index (n^2 > 0) #
        #----------------------------------#
        valid = n_sq > 0
        curr_n_sq = n_sq[valid]
        curr_om = om[mask_ff][valid]
        curr_ss = ss[mask_ff][valid]
        curr_dd = dd[mask_ff][valid]
        curr_pp = pp[mask_ff][valid]

        #-----------------------------#
        # wavenumber through n-square #
        #-----------------------------#
        ak = np.sqrt(curr_n_sq) * curr_om

        #-------------------------------------------------#
        # polarization and wave electric field components #
        #-------------------------------------------------#
        pol = (curr_n_sq - curr_ss) / np.where(np.abs(curr_dd) > 1e-12, curr_dd, 1e-12)
        ez = -curr_n_sq * np.sqrt(cos2 * sin2) / (curr_pp - curr_n_sq * sin2 + 1e-12)
        ey2 = 1 / (pol**2 + 1e-12)

        #----------------#
        # ES-to-EM ratio #
        #----------------#
        denom = np.sqrt(1 + ey2 + ez**2)
        esem = (np.sqrt(sin2) + np.sqrt(cos2) * ez) / denom

        #--------------------------------------------------------#
        # log10 ratio (safe-guarded against non-positive values) #
        #--------------------------------------------------------#
        ratio_arg = np.maximum(1e-12, 1 - esem**2)
        val4 = np.log10(np.maximum(1e-12, esem / np.sqrt(ratio_arg)))

        #---------------------------------------------#
        # Save: [k, omega, polarization, log10_ratio] #
        #---------------------------------------------#
        results.append(np.column_stack((ak, curr_om, pol, val4)))

    #------------------------#
    # convert to NumPy array #
    #------------------------#
    if not results:
        print("No valid physical solutions found.")
        return

    #-----------------------------#
    # combine roots and transpose #
    #-----------------------------#
    omsave = np.vstack(results).T

    #-------------------------------#
    # ion-scale conversion of units #
    #-------------------------------#
    freq = omsave[1, :] * mass_ratio / wce
    wavenum = omsave[0, :] * np.sqrt(mass_ratio)
    wmax_plot = wmax * mass_ratio / wce
    kmax_plot = 15 * np.sqrt(mass_ratio)

    #---------------------------------#
    # numpy save arrays in npz format #
    #---------------------------------#
    np.savez(outfile, wavenum=wavenum, freq=freq)

    #----------#
    # plotting #
    #----------#
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    axes[0].plot(wavenum, freq, '.', markersize=1)
    axes[0].set_xlim(0, kmax_plot)
    axes[0].set_ylim(0, wmax_plot)
    axes[0].set_title('Dispersion diagram')
    axes[0].set_xlabel('$kc/\omega_\mathrm{pi}$')
    axes[0].set_ylabel('$\omega/\Omega_\mathrm{i}$')
 

    axes[1].plot(omsave[2, :], freq, '.', markersize=1)
    axes[1].set_xlim(-10, 10)
    axes[1].set_ylim(0, wmax_plot)
    axes[1].set_title('Wave polarization')
    axes[1].set_xlabel('polarization')

    axes[2].plot(omsave[3, :], freq, '.', markersize=1)
    axes[2].set_xlim(-5, 5)
    axes[2].set_ylim(0, wmax_plot)
    axes[2].set_title('Es/Em')
    axes[2].set_xlabel('log10(ES/EM)')

    plt.tight_layout()
    plt.savefig(pdffile)
    plt.show()

if __name__ == "__main__":
    dispcold()

