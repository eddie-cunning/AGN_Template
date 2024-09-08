***************************
 AGN Identified in ZFOURGE
***************************

RELEASE: V0.5
DATE: JULY 28th, 2016
SOURCE: http://zfourge.tamu.edu/
EMAIL: michael.cowley@students.mq.edu.au

**************************
         DETAILS
**************************

AGN flags for the ZFOURGE v3.4 catalogues

**************************
        REVISIONS
**************************

V0.1 ---------------------  Initial Release
V0.2 ---------------------  Formatting
V0.3 ---------------------  Recut using updated photo-zs
V0.4 ---------------------  Formatting
V0.5 ---------------------  Recut using updated photo-zs and use flags

**************************
    CATALOGUE ENTRIES
**************************

id  ----------------------  ZFOURGE ID
ir_agn  ------------------  -1 = no data, 0 = nil detection, 1 = positive detection
radio_agn  ---------------  -1 = no data, 0 = nil detection, 1 = positive detection
xray_agn  ----------------  -1 = no data, 0 = nil detection, 1 = positive detection

NOTES:

(1)

Infrared AGN (ir_agn) are selected using the colour-colour classifications of Messias et al. (2012):
       
z < 1.8: Ks - [4.5] > 0, [4.5] - [8.0] > 0
z > 1.8: [8.0] - [24] > 2.9 x ([4.5] - [8.0]) + 2.8, [8.0] - [24] > 0.5

(2)

Radio AGN (radio_agn) are selected using the Rees et al. (2015) Radio AGN activity index:

SFR_{Radio} / SFR_{IR+UV} > 3.0

The minimum root-mean-square (RMS) sensitivity for the 1.4GHz data in CDFS, COSMOS and UDS fields are 6, 10 and 100 μJy/beam, respectively.

(3)

X-ray AGN (xray_agn) are selected using the Szokoly et al. (2004) classification:

L_x >= 1e41 erg/s & HR > -0.2
L_x >= 1e42 erg/s & HR <= -0.2

The on-axis limiting X-ray flux in the soft and hard bands for the CDFS, COSMOS and UDS fields are 9.1e−18 and 5.5e−17 erg/cm^2/s, 1.9e−16 and 7.3e−16 erg/cm^2/s, and 6.0e−16 and 3.0e−15 erg/cm^2/s, respectively. 

(4)

AGN selection requires use_flag==1 as the AGN selection techniques require a reliable redshift. 

(5)

More details can be found in Cowley at al. (2016).