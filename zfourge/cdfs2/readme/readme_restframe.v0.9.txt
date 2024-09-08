**************************
ZFOURGE REST FRAME CATALOGS
**************************
Compiled by Ben Forrest ----- bforrest@physics.tamu.edu

ZFOURGE v3.4 catalogs run through EAZY for various filters
  cdfs.v1.6.9
  cosmos.v1.3.6 
  uds.v1.5.8 
EAZY downloaded from github on 30 October, 2015.

Each rest frame flux is fit individually.  This yields slightly different results compared to fitting multiple fluxes simultaneously.
Reported errors are 1-sigma spreads on the flux determined from the width of the p(z) as reported by EAZY 1/2*((p84-L_#) + (L_# - p16)).
UV filters are simple tophat filters. Their response curves are unity to 175* angstroms on either side of the central wavelength for which it is named and zero elsewhere, i.e.:

r = / 1, if |lambda_c - lambda| < 175
    \ 0, if |lambda_c - lambda| > 175

All filter curves can be found in the FILTER_RES_UV file.
*The uv13 filter is only 110A wide to avoid contamination from Lya.


**************************
        REVISIONS
**************************

v0.1 ----------------------  Initial Release
v0.2 ----------------------  Added units, filter central wavelengths
v0.3 ----------------------  Author info added to README; naming scheme altered
v0.4 ----------------------  README description improved
v0.5 ----------------------  Previous versions utilized earlier use flags. 
                             This run is updated to the catalog versions listed above.
v0.6 ----------------------  Error calculation updated.
v0.7 ----------------------  Units changed to mJy 
v0.8 ----------------------  Units changed to uJy, EAZY runs with zbin file
v0.9 ----------------------  Updated to data release v3.4.
                             z_spec utilized, uv13 changed to 110A wide.

**************************
    CATALOG ENTRIES
**************************

id        ----  id number from ZFOURGE v3.4 catalog for a given field
f_* [uJy] ----  flux in the rest frame filter *
e_* [uJy] ----  error in the flux of rest frame filter *

FILTERS: (filter number in FILTER.RES_UV)
u  -  SDSS u filter (73)
g  -  SDSS g filter (74)
r  -  SDSS r filter (75)
i  -  SDSS i filter (76)
z  -  SDSS z filter (77)
U  -  Johnson U filter (153)
B  -  Johnson B filter (154)
V  -  Johnson V filter (155)
J  -  2MASS J filter (161)
H  -  2MASS H filter (162)
K  -  2MASS K filter (163)
uv13  -  tophat filter 110 angstroms in width, centered at 1300 angstroms (293)
uv15  -  tophat filter 350 angstroms in width, centered at 1500 angstroms (294)
uv19  -  tophat filter 350 angstroms in width, centered at 1900 angstroms (295)
uv22  -  tophat filter 350 angstroms in width, centered at 2200 angstroms (296)
uv28  -  tophat filter 350 angstroms in width, centered at 2800 angstroms (297)


**************************
    CHANGING UNITS
**************************
These catalogs are in units of uJy.
To convert back to catalog fluxes (AB flux densities with zeropoint 25), divide by 10**((23.9-25)/2.5).
To convert to solar luminosities:
  1. Convert to ergs/cm**2/s/A -- flux density [uJy] *10**-29 *3*10**18 /lambda_c**2,
     where lambda_c is the central wavelength of the filter bandpass.
  2. Convert from flux density to flux -- multiply by lambda_c
  3. Convert from flux to luminosity using the distance modulus
     ( distance = 10**((DM+5)/5)*3.086*10**18 ) --
     L = flux *4pi *distance**2 /(1+z)
  4. Convert to solar luminosities -- L /(3.828*10**33)



