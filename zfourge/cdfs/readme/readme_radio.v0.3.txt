**************************
 Radio Sources in ZFOURGE
**************************

RELEASE: V0.3
DATE: JULY 28TH, 2016
SOURCE: http://zfourge.tamu.edu/
EMAIL: michael.cowley@students.mq.edu.au

**************************
         DETAILS
**************************

ZFOURGE v3.4 crossmatched with published radio sources, based on VLA data:

cdfs.v1.6.9  -------------  VLA 1.4 GHz Survey of the Extended Chandra Deep Field South: 2nd DR (Miller et al. 2013, M13)
cosmos.v1.3.6  -----------  VLA-COSMOS Survey IV Deep Data and Joint Catalogue (Schinnerer et al. 2010, S10)
uds.v1.5.8  --------------  Subaru/XMM-Newton Deep Field-I 100 μJy Catalogue (Simpson et al. 2006, S06)

**************************
        REVISIONS
**************************

V0.1 ---------------------  Initial Release
V0.2 ---------------------  Recalculated Luminosities using updated photo-zs
V0.3 ---------------------  Formatting
v0.4 ---------------------  Recalculated Luminosities using updated photo-zs

**************************
    CATALOGUE ENTRIES
**************************

id  ----------------------  ZFOURGE ID
f_14  --------------------  Peak or integrated flux density at 1.4GHz (uJy/beam)
e_14  --------------------  Error in f_14 (uJy/beam)
l_14  --------------------  Rest-frame 1.4GHz luminosity (W/Hz)

NOTES:

(1)

A total of 370 radio sources (218 in CDFS, 122 in COSMOS, 30 in UDS) were cross-matched with the primary ZFOURGE catalogues

(2)

Radio sources were matched to ZFOURGE to within 1” and positionally corrected by calculating the median positional uncertainty in each field. In arcsecs, the corrections are CDFS: RA = -0.19 DEC = -0.23, COSMOS: RA = +0.08 DEC = +0.07, UDS: RA = +0.03 DEC = +0.37. 

(3)  

For M13, f_14 includes a combination of peak and integrated fluxes, with the final value chosen based on the author’s recommendation (see Section 3.4 for further details). For S10 and S06, the integrated flux is used. 

(4)

The rest-frame radio luminosity is calculated using:

	l_14 = (4*pi*ld**2)*f_14*(1+z_peak)**(α-2)

where ld is the luminosity distance (cm) and α is the radio spectral index (fixed to α = -0.3)
