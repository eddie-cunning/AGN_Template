
####################################################
#########   Star-Formation Rate Catalogs   #########
####################################################

AUTHOR: Adam R. Tomczak
EMAIL: artomczak@ucdavis.edu
RELEASE: v0.4
DATE: July, 27th 2016
SOURCE: http://zfourge.tamu.edu/



###########################################
#########   Far-IR Imaging Data   #########
###########################################

cdfs.v1.6.9
  - 24um: Spitzer MIPS (GOODS: Dickinson et al. 2007)
  - 100um: Herschel PACS (Elbaz et al. 2011)
  - 160um: Herschel PACS (Elbaz et al. 2011)

cosmos.v1.3.6
  - 24um: Spitzer MIPS (S-COSMOS: Sanders et  al. 2007)

uds.v1.5.8
  - 24um: Spitzer MIPS (SpUDS: PI Dunlop)



#######################################
#########   Catalog Details   #########
#######################################

Star-formation rates are derived from adding the contribution of UV and IR emission. All far-IR fluxes are in units of milli-Janskys (AB zeropoint 16.4) and luminosities are in units of solar luminosities.

Far-IR fluxes were measured using a custom source deblending code written by I. Labbe. For each source, neighboring objects detected in the Ks band imaging are modeled at the resolution of each far-IR image and subtracted before measuring photometry. For more details regarding this procedure see Wuyts et al. (2007).

Bolometric IR luminosities (LIR) are estimated by redshifting the IR spectral template used by Wuyts et al. (2008) to a galaxy's photometric redshift (z_peak), scaling to the measured 24um photometry, and integrating between rest-frame 8-1000um. For cdfs.v1.6.2.herschel.cat, the IR template is instead scaled based on a weighted least-squares fit to the 24, 100, and 160um photometry.

Rest-frame 2800 Angstrom luminosities (L2800) are derived using EAZY to fit templates near rest-frame 2800 Angstroms. We use top-hat filter centered on 2800 Angstroms of width 350 Angstroms.

UV+IR star-formation rates are derived using the calibration of Bell et al. (2005) scaled to a Chabrier (2003) IMF:

    SFR = 1.09e-10 * (LIR + 2.2*LUV)
    LUV = 1.5 * nu * L2800



#############################
#########   Files   #########
#############################

readme_sfr.txt
cdfs.v1.6.9.herschel.cat
cdfs.v1.6.9.sfr.cat
cosmos.v1.3.6.sfr.cat
uds.v1.5.8.sfr.cat

./data/MipsFilter_24um.txt     Transmission curve of MIPS 24um filter
./data/PacsFilter_100um.txt    Transmission curve of PACS 100um filter
./data/PacsFilter_160um.txt    Transmission curve of PACS 160um filter
./data/template_W08.txt        Infrared template of Wuyts et al. (2008)



############################################
#########   Fields for *.sfr.cat   #########
############################################

id                           ZFOURGE id
z_peak                       photometric redshift
f24        [ mJy ]           total MIPS 24um flux
e24        [ mJy ]           total MIPS 24um flux error
LIR24      [ Lsol ]          integrated 8-1000um luminosity
L2800      [ Lsol ]          rest-frame 2800A luminosity
SFR_UVIR   [ Msol / yr ]     UV+IR star-formation rate

Note: A value of -99 is reserved for objects with
      missing data (i.e. objects that fall off the
      footprint of the far-IR images).



#################################################
#########   Fields for *.herschel.cat   #########
#################################################

id                           ZFOURGE id
z_peak                       photometric redshift
f24        [ mJy ]           total MIPS 24um flux
e24        [ mJy ]           total MIPS 24um flux error
f100       [ mJy ]           total PACS 100um flux
e100       [ mJy ]           total PACS 100um flux error
f160       [ mJy ]           total PACS 160um flux
e160       [ mJy ]           total PACS 160um flux error
LIR        [ Lsol ]          integrated 8-1000um luminosity
L2800      [ Lsol ]          rest-frame 2800A luminosity
SFR_UVIR   [ Msol / yr ]     UV+IR star-formation rate

Note: A value of -99 is reserved for objects with
      missing data (i.e. objects that fall off the
      footprint of the far-IR images).



##################################
########    References    ########
##################################

Bell et al. (2005), ApJ, 625, 23
Dickinson & FIDEL Team (2007), BAAS, 39, 822
Elbaz et al. (2011), A&A, 533, A119
Sanders et al. (2007), ApJS, 172, 86
Wuyts et al. (2007), ApJ, 655, 51
Wuyts et al. (2008), ApJ, 682, 985

