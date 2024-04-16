**************************
 X-ray Sources in ZFOURGE
**************************

RELEASE: V0.4
DATE: JULY 28TH, 2016
SOURCE: http://zfourge.tamu.edu/
EMAIL: michael.cowley@students.mq.edu.au

**************************
         DETAILS
**************************

ZFOURGE v3.4 crossmatched with published X-ray sources, based on Chandra and XMM-Newton data:

cdfs.v1.6.9  -------------  The CDF-S Survey 4ms Source Catalogs (Xue et al. 2011, X11)
cosmos.v1.3.6  -----------  The C-COSMOS Survey I (Elvis et al. 2009, E09)
uds.v1.5.8  --------------  The SXDS III X-ray Data (Ueda et al. 2008, U08)

**************************
        REVISIONS
**************************

V0.1 ---------------------  Initial Release
V0.2 ---------------------  Formatting
V0.3 ---------------------  Recalculated Luminosities using updated photo-zs
v0.4 ---------------------  Recalculated Luminosities using updated photo-zs

**************************
    CATALOGUE ENTRIES
**************************

id  ----------------------  ZFOURGE ID
f_xray  ------------------  Full band X-ray flux (erg/cm2/s)
l_xray  ------------------  Rest-frame X-ray luminosity (erg/s)
HR  ----------------------  Hardness ratio

NOTES:

(1)

A total of 736 X-ray sources (596 in CDFS, 109 in COSMOS, 31 in UDS) were cross-matched with the primary ZFOURGE catalogues

(2)

X-ray sources were matched to ZFOURGE to within 1” and positionally corrected by calculating the median positional uncertainty in each field. In arcsecs, the corrections are CDFS: RA = -0.15 DEC = -0.26, COSMOS: RA = +0.13 DEC = +0.14, UDS: RA = +0.02 DEC = -0.16. 

(3)  

In the absence of a full-band flux value, f_xray is calculated considering a counts-to-flux conversion factor with photon index Γ = 1.4

(4)

For sources from X11, the intrinsic flux is derived from counts in the 0.5-8 keV full band, while from E09 and U08 it is derived from the sum of the counts in the relevant bands over 0.5-10 keV. Flux values were adjusted for sources from E09 and U08 to align with the full bandpass values of X11 by assuming a power-law model with photon index Γ = 1.4 (i.e. E09 and U08 fluxes are multiplied by a factor of 0.95)

(5)  

The rest-frame X-ray luminosity is calculated using:

	l_xray = (4*pi*ld**2)*f_x*(1+z_peak)**(Γ-2)

where ld is the luminosity distance (cm) and Γ is the photon index (fixed to Γ = 1.4)

(6)  

The hardness ratio (HR) is calculated using:

	HR = (H-S)/(H+S)
 
where H and S are the count rates in the hard (2-8 or 2-10 keV) and soft (0.5-2 keV) bands, respectively

(7)  

The following AGN classification criteria are adopted in Cowley at al. (2015):

	(l_xray >= 10**41) & (HR > −0.2)
	(l_xray >= 10**42) & (HR <= −0.2)

