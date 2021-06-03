      PROGRAM CTSFG5              

***** versione 2 marzo 2005

* ATTENZIONE A X10 e X23 SE USI MISFIT
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* GAS SYSTEM IS H-O-S. Carbon is not accounted for.....                   *
*                                                                         *
* computes PO2 ,            given T,Ptot and Fe2+/Fe3+                    *
*          PO2 ,            given T,Ptot and Fe2+/Fe3+                    *
* or....   Fe2+/Fe3+,       given PO2, Ptot, T                            *
* or....   T,               given PO2 and Fe2+/Fe3+                       *
* or the sulphide capacity, given PO2/PS2 or FeII/FeIII, PS2  and T       *
* or the sulphate capacity, given pO2/pS2 or FeII/FeIII, pS2  and T       *
* or dissolved sulphur,     given  pO2 or FeII/FeIII, pS2 and T           *
* or PS2, PSO2, PSO3, PSO,  given dissolved sulphide and PO2 and T        *
*    """""""""""""""""""",  given dissolved sulphide and FeII/FeIII       *
*    """""""""""""""""""",  given dissolved sulphate and PO2 and T        *
*    """""""""""""""""""",  given dissolved sulphate and FeII/FeIII       *
*                                                                         *
* Il rapporto Fe''/Fe''' da considerare Š REDOZ !                         *
*                                                                         *
* N.B. N.B. N.B. N.B. N.B. N.B. N.B. N.B.:                                *
*                                                                         *
* PUT logfO2 = 0. TO GET logfO2 as OUTPUT                                 *
*                                                                         *
* N.B. N.B. N.B. N.B. N.B. N.B. N.B. N.B.:                                *
*                                                                         *
*                                                                         *
* It adopts Moretti & Ottonello model for water speciation                *           
*_________________________________________________________________________*
* based on the simplified polymeric approach of
* Ottonello et al. (2000) for the evaluation of redox conditions
* and on the Moretti & Ottonello (2002) S= solubility model
* and on the Moretti & Ottonello (2002) total S solubility model
*_________________________________________________________________________
*                                                                         *
* ancillary informations                                                  *
*                                                                         *
* Z: basicity moderating parameters (eq. 29-II in Ottonello et al. (2000))*
*     exp. data of Young et al. (1992) with modifications given in        *
*     Moretti & Ottonello (2002) (e.g. water set to 2.56)                 *
*                                                                         *
*  Y: oxides gfw                                                          *
*                                                                         *
*  Zcat: n. cations per formula unit                                      *
*                                                                         *
*  X: a,b constants in a lnK=a+b/T equation of type:                      *
*     AO + 0.5S2 <==> AS + 0.5O2                                          *
*     on a 1-cation basis for sulphide capacity                           *
*                                                                         *
* XX: AO + 0.5S2 + 1.5O2 <==> ASO4                                        *
*     on a 1-cation basis for sulphate capacity                           *
*                                                                         *
* volx: molar volumes of oxides in the melt phase (cc/mol)                *
*                                                                         *
* eoxm: expansivity of oxides in melt phase (cc/k mol)                    *
*                                                                         *
* coxm: compressibility of oxides in melt phase                           *
*                                                                         *
* Reference standard state is 1 bar @ T of interest                       *
*                                                                         *
* Fugacities are employed. Activity terms in the liquid phase are         *
* corrected for pressure greater than 1 bar.                              *
*                                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     DECLARATIONS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      common /setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),  XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,     A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         regg(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),focz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi), sobs(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ,ISYSMF,isysvo,
     &         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
*     UNIT NUMBERS
 
      ISYSRD = 17   ! Read unit
      IVIDEO =  6
      ISYSWR =  8
      ISYSIN =  9
      ISYSCK = 12
      ISYSQZ = 13
      ISYSWZ = 14
      ISYSPR = 15
      ISYSPU = 16
      ISYSMF = 19
      ISYSTT=  21
      ISYSVO = 18
      ISYSHO = 11
      ISYSFE = 10
 
c     open(File='1.IN'        ,unit=ISYSRD,status='UNKNOWN')
c     open(File='tange2.IN'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='tanghero.in' ,unit=ISYSRD,status='UNKNOWN')
c     open(File='tangeman.IN' ,unit=ISYSRD,status='UNKNOWN')
c     open(File='kd.IN'       ,unit=ISYSRD,status='UNKNOWN')
c     open(File='metrich.IN'  ,unit=ISYSRD,status='UNKNOWN')
c     open(File='metrich2.IN' ,unit=ISYSRD,status='UNKNOWN')
c     open(File='inclus.IN'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='mutual.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='empty.in'    ,unit=ISYSRD,status='UNKNOWN')
c     open(File='wempty.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='1wempty.in'  ,unit=ISYSRD,status='UNKNOWN')
c     open(File='qempty.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='kempty.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='redox.in'    ,unit=ISYSRD,status='UNKNOWN')
c     open(File='kfs.in'      ,unit=ISYSRD,status='UNKNOWN')
c     open(File='nafs.in'     ,unit=ISYSRD,status='UNKNOWN')
c     open(File='cafs.in'     ,unit=ISYSRD,status='UNKNOWN')
c     open(File='casio2.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='138gut.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='K2Oadd.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='CaOadd.in'   ,unit=ISYSRD,status='UNKNOWN')
c     open(File='Na2Oadd.in'  ,unit=ISYSRD,status='UNKNOWN')
c     open(File='visco.in'    ,unit=ISYSRD,status='UNKNOWN')
      open(File='INPUT.txt'   ,unit=ISYSRD,status='OLD')
      open(File='ctsfg6.dat'  ,unit=ISYSQZ,status='UNKNOWN')
      open(File='ctsfg6.jet'  ,unit=ISYSPU,status='UNKNOWN')
      open(File='ctsfg6.s2'   ,unit=ISYSWR,status='UNKNOWN')
      open(File='ctsfg6.tt'   ,unit=ISYSTT,status='UNKNOWN')
      open(File='ctsfg6.inp'  ,unit=ISYSIN,status='UNKNOWN')
      open(File='ctsfg6.log'  ,unit=ISYSPR,status='UNKNOWN')
      open(File='ctsfg6.so4'  ,unit=ISYSWZ,status='UNKNOWN')
      open(File='ctsfg6.fit'  ,unit=ISYSMF,status='UNKNOWN')
      open(File='ctsfg6.vol'  ,unit=Isysvo,status='UNKNOWN')
      open(File='ctsfg6.chk'  ,unit=Isysck,status='UNKNOWN')
      open(File='ctsfg6.h2o'  ,unit=Isysho,status='UNKNOWN')
      open(File='ctsfg6.fer'  ,unit=Isysfe,status='UNKNOWN')
  
      write (*,*) ' Do you want to fix fSO2 and fO2 ? (y=1) '
      read   (*,*)  answ
 
      CALL NEWSULF
 
c     write(ivideo,*)
c     write(IVIDEO,*) 'data have been correctly read'
c     write(ivideo,*)
c     pause
 
      call SULPHIDE
      call SULPHATE
 
      END
 
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo,
     &         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      DATA Z/2.09,1.54,2.5,1.67,2.09,1.82,1.354,1.40,1.,1.28,.87,
     & .71,2.50/
c N.B gamma per Mn Š 1.69 nel lavoro precedente e per H2O 2.56!
      DATA Y/60.0848,79.899,141.94,101.96,159.69,151.99,71.846,
     & 70.937,56.079,40.311,61.979,94.203,18.015/
      DATA Zcat/1.,1.,2.,2.,2.,2.,1.,1.,1.,1.,2.,2.,2./
         
***********************************************************************
cc Best water parameter would be -2.607 if we use Burnham model...
cc  Old regression, valable with old volumes. Papale value inserted
cc    DATA X/12.752,0.27953,1.7132,7.70870,8.279,0.,6.1184,-0.917,
cc   $ -4.0239,-4.1698,14.5300,11.266,-0.8281,-36418.,
cc   $-23460.,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
cc   $-15181.,-16396.,-16432.,-18764./
cc
cc New values 21 October 2001 from nohrs2.out in mawnew
cc    DATA X/12.752,2.1962,1.7132,7.70870,14.278,0.,5.6,-0.917,
cc   $ -4.0239,-4.1698,14.5300,11.266,-1.457,-36418.,
cc   $-28698.63,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
cc   $-15181.,-16396.,-16432.,-18764./
cc
cc Last values 18 February 2002. con S of annealing "P dependent".
cc determined with a Clausius - Clapeyron approximation.
cc Water value is referred to Papale H2O model (-148.98233 for all data
cc Burnham water best values would be: -126.556167 for all data.
cc    DATA XX/-26.436,-33.625,0.,-141.784,-55.963,0.,-35.4375,
cc   $-42.0084,-0.601,-36.708,-22.089,13.7,
cc below Papale model derived value
cc   $-148.98233,77291.1,77291.1,0.,
cc below Burnham model derived value
cc   $-126.556167,77291.1,77291.1,0.,
cc   $293557.81,132676.85,0.,85972.16,85972.16,27384.29,95831.,
cc   $100934.,100741.,67292./
cc
cc  New values (10 aprile 2002)
cc    DATA XX/-26.436,-25.4631,0.,-30.8846,-31.411,0.,-42.4136,
cc   $-49.7195,-41.8793,-36.8460,-22.089,13.7,
cc below Papale model derived value
cc   $-148.98233,77291.1,98230.0,0.,
cc below Burnham model derived value
cc   $-126.556167,77291.100,98230.000,0.,
cc   $85879.800,85879.800,0.,99540.800,106414.10,104207.80,95831.,
cc   $100934.,100741.,67292./
 
***********************************************************************
c    NO MORE GOOD VALUES ARE BELOW (3 june 2002)
***********************************************************************
cc   VERY IMPORTANT: Verify x27 and x227 (i.e. the annealing term -alfaV)
 
cc  OXIDE-SULFIDE disproportions: new values 14 April 2002
c      DATA X/12.752,-4.7934,1.7132,7.70870,15.0328,0.,5.4414,-0.917,
cc  Below Papale value without P effect
cc    $ -4.0239,-4.1698,14.5300,11.266,-2.7750,-36418.,
cc  Below Papale value with P effect (S annealing)
c*    $ -4.0239,-4.1698,14.5300,11.266,-19.9782,-36418.,
cc  Below Burnham value without P effect
c     $ -4.0239,-4.1698,14.5300,11.266,-6.06220,-36418.,
cc  Below Burnham value with P effect (S annealing)
c*    $ -4.0239,-4.1698,14.5300,11.266,-15.0592532,-36418.,
c     $-16575.9,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
c     $-15181.,-16396.,-16432.,-18764./
 
cc   OXIDE_SULFATE disproportions: new values (14 aprile 2002)
c      DATA XX/-26.436,-27.0906,0.,-32.98960,-33.89291,0.,-40.958569,
c     $-60.128854,-42.38865,-36.408443,-22.089,13.70000,
cc  below Papale model derived value without P effect
cc    $-65.985,77291.1,98230.0,0.,
cc  below Papale model derived value with P effect (S annealing)
c*    $-87.044140,77291.1,98230.0,0.,
cc  below Burnham model derived value without P effect
c     $-60.852000,77291.100,98230.000,0.,
cc  below Burnham model derived value with  P effect (S annealing)
c*    $-87.044140,77291.100,98230.000,0.,
c     $85879.800,85879.800,0.,99540.800,106414.10,107754.894,95831.,
c     $100934.,100741.,67292./
**********************************************************************
c    GOOD VALUES ARE BELOW (27 november 2003)
***********************************************************************
c   VERY IMPORTANT: Verify 3 and x227 (i.e. the annealing term -alfaV)
 
c  OXIDE-SULFIDE disproportions: new values 23 March 2004 for H2O
cf      DATA X/12.752,-4.7934,1.7132,7.70870,15.0328,0.,5.4414,-0.917,
ccf  Below Papale value without P effect
ccf    $ -4.0239,-4.1698,14.5300,11.266,0.0,-36418.,
ccf  Below Papale value with P effect (S annealing)
cf     $ -4.0239,-4.1698,14.5300,11.266,1.34446741,-36418.,
ccf  Below Burnham value without P effect
ccf    $ -4.0239,-4.1698,14.5300,11.266,-4.069,-36418.,
ccf  Below Burnham value with P effect (S annealing)
*cf    $ -4.0239,-4.1698,14.5300,11.266,-4.069,-36418.,
cf     $-16575.9,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
cf     $-15181.,-16396.,-16432.,-18764./
 
c   OXIDE_SULFATE disproportions: new values (23 March 2004)
cf      DATA XX/-26.436,-27.0906,0.,-32.98960,-33.89291,0.,-40.958569,
cf     $-60.128854,-42.38865,-36.408443,-22.089,13.70000,
ccf  below Papale model derived value without P effect
ccf    $-65.4839395,77291.1,98230.0,0.,
ccf  below Papale model derived value with P effect (S annealing)
cf     $-120.683,77291.1,98230.0,0.,
ccf  below Burnham model derived value without P effect
ccf    $-60.8531,77291.100,98230.000,0.,
ccf  below Burnham model derived value with  P effect (S annealing)
*cf    $-105.488,77291.100,98230.000,0.,
cf     $85879.800,85879.800,0.,99540.800,106414.10,107754.894,95831.,
cf     $100934.,100741.,67292./
***********************************************************************************

c   VERY IMPORTANT: Verify x27 and x227 (i.e. the annealing term -alfaV)
 
c  OXIDE-SULFIDE disproportions: new values 13 August 2004 for all oxide
c  GCA paper revision followoin reviewers
c REVISED ON 7 APRIL 2005
      DATA X/12.9576967,-4.8017622,2.15696154,7.81050803,6.66128037,0.,
     $8.41324988,-1.0201951,-4.2289568,-5.0761932,11.253287,10.6606042,
*  Below Papale value (only without P effect) 7 APRIL 2005
     $ 0.781537236943571,
*  Below Burnham value (only without P effect) 9 Aug 2004
*    $ 11.7851076580539,
     $-36418.,-16575.9,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,
     $-6804.2,-15181.,-16396.,-16432.,-18764./
 
c   OXIDE_SULFATE disproportions: new values 13 August 2004 for all ox except H2O
c   GCA paper revision followoin reviewers
c REVISED ON 7 APRIL 2005
      DATA XX/-26.883627,-22.661437,0.,-31.095548,-31.927338,
     $ 0.,-42.051629,-51.400798,-43.517305,
     $ -36.371917,-21.612439,14.8139739,
*    BELOW PAPALE VALUE NO ANNEALING 13 Aug 2004
*     $-60.3993819146732,
*    BELOW PAPALE VALUE ANNEALING INCLUDED 7 April 2005
     $-119.74,   
*    BELOW BURNHAM VALUE NO ANNEALING 13 Aug 2004
*     $-56.3439946571902,
*    BELOW BURNHAM VALUE ANNEALING INCLUDED 13 Aug 2004
*     $-101.959,
     $77291.1,98230.0,0.,85879.800,85879.800,0.,99540.800,106414.10,
     $107754.894,95831.,100934.,100741.,67292./
***********************************************************************************

* MODIFIED TO SEE ERROR PROGRESSION (Parameters are modified in the same direction,
* both positive and negative)  NOW POSITIVE
* A + B/T former   = logK: if they increase log Cs and Cso4 increase
* A + B/T modifier = logK: if they increase log Cs and Cso4 increase
 
c   VERY IMPORTANT: Verify x27 and x227 (i.e. the annealing term -alfaV)
c   SiO2 B' term was locked at 77291.1. If introduced goes to 83000
 
*c  OXIDE-SULFIDE disproportions: new values 14 April 2002
*      DATA X/12.822,-1.1307,2.04933,7.89874,16.2314,0.,5.70743,-0.84616,
*c  Below Papale value without P effect
*     $-3.92242,-4.0347,14.822,11.71141,-1.3772,-36418.,
*     $-10134,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
*     $-15181.,-16396.,-16432.,-18764./
 
*c   OXIDE_SULFATE disproportions: new values (14 aprile 2002)
*      DATA XX/-23.2668,-21.8977,0.,-32.1198,-33.15832,0.,-40.2435,
*     $-56.0869,-36.2243,-35.3041,-21.28064,14.9927,
*c below Papale value without P effect
*     $-65.11081,77291.1,98230.0,0.,
*     $85879.800,85879.800,0.,99540.800,106414.10,107754.894,95831.,
*     $100934.,100741.,67292./
 
***********************************************************************************
* MODIFIED TO SEE ERROR PROGRESSION (Parameters are modified in the same direction,
* both positive and negative)  NOW NEGATIVE
* A + B/T former   = logK: if they decrease log Cs and Cso4 decrease
* A + B/T modifier = logK: if they decrease log Cs and Cso4 decrease
 
c   VERY IMPORTANT: Verify x27 and x227 (i.e. the annealing term -alfaV)
c   SiO2 B' term was locked at 77291.1. If introduced goes to 72000
 
*c  OXIDE-SULFIDE disproportions: new values 14 April 2002
*      DATA X/12.682,-8.4561,1.37707,7.51866,13.8346,0.,5.17537,-0.98784,
*c  Below Papale value without P effect
*     $-4.12538,-4.30473,14.238,10.8206,-4.173,-36418.,
*     $-23017.9,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
*     $-15181.,-16396.,-16432.,-18764./
 
*c   OXIDE_SULFATE disproportions: new values (14 aprile 2002)
*      DATA XX/-29.6082,-32.2843,0.,-33.8598,-34.62768,0.,-41.67453,
*     $-64.1711,-48.5537,-37.5119,-22.88836,12.4073,
*c  below Papale value without P effect
*     $-66.85919,77291.1,98230.0,0.,
*     $85879.800,85879.800,0.,99540.800,106414.10,107754.894,95831.,
*     $100934.,100741.,67292./
 
****************************************************************************
c NEW VOLUME REGRESSION (file Volumes_melts.xls) Expansivities for SiO2 and Al2O3
c have been adjusted on the basis of the OptBas correlation. Water Expansivity obtained
c by considering H2O estimated V at 1673K and Richet Volume at room T (12 cc/mol).
c About expansivities:
c SiO2 : use 6.24 (recalculated) or 0.10 (Lange)
c Al2O3: use 8.10 (recalculated) or 2.62 (Lange)
c P2O5:  use 18.3 (recalculated) or 2.62 (P2O5 as Al2O3 of Lange)
      data voxm/26.9,23.16,82.16,37.11,42.13,36.36,13.65,11.62,16.57,
     $          11.45,28.78,45.84,16.79/
      data eoxm/0.10,7.24,2.62,2.62,9.09,8.22,2.92,2.73,2.92,2.62,7.41,
     $          11.91,3.48/
      data coxm/-1.89,-2.31,-8.93,-2.26,-2.53,-1.96,-0.45,-0.37,-1.34,
     $          -0.40,-2.4,-6.75,-1.88/
***********************************************************************
c  segue per calcolo porosit… ionica
 
      data noxy/2,2,5,3,3,3,1,1,1,1,1,1,1/
      data radi/0.26,0.605,0.17,0.39,0.49,0.615,0.78,0.67,1.,0.72,1.02,
     $          1.37,0./
 
      END
 
      SUBROUTINE NEWSULF
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* this is the main program, disguised as a subroutine for *
* reasons of compatibility between systems.               *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      common/setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      Read (isysrd,*) ncomp
	  Read (isysrd,*) ! empty line to skip header
      write (isysin,*) 'SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO MnO MgO CaO Na2O
     & K2O P2O5 H2O'
C     rapmed=0.0
      do 111 i=1,ncomp
C     jflag(i)=0
      Read (isysrd,*) poss(i,1),poss(i,2),poss(i,4),poss(i,5),poss(i,6),
     &poss(i,7),poss(i,8),poss(i,10),
     &poss(i,9),poss(i,11),poss(i,12),poss(i,3), poss(i,13),
     &stot(i),so2(I),tcent(i),pbar(i),xossi(i),fs2(i),iaut(I),kflag(i),
     &wmol(i),kir(i)

        if (wmol(i).eq.0.) then
         wmol(i)=0.1*poss(i,13)
        endif
	wmol(i)=wmol(i)/18.015 ! this is molecular water in moles
c        write (*,*) wmol(i)
c        pause
c      pbar(i)=1.
c       poss(i,13)=0.
c Nota Bene:
c kflag=0 il programma itera per Fe2+/Fe3+
c kflag=2 il programma accetta Fe2+/Fe3+ di input e log fO2
c      kflag(i)=2
 
      if (answ.eq.1.) then
       write (*,*)
       write (*,*) '-------------------------------------------------'
       write  (*,*) 'set fSO2 (bar) @ T = ',tcent(i),' P = ',pbar(i),
     & ' for composition n. ',i
       write (*,*)
       read (*,*) fso2
c      fso2=0.021
       fso2=dlog10(fso2)
798    write  (*,*) 'will you set dNNO or dFMQ or logfo2 (1/2/3)?'
       write  (*,*)
       read (*,*) kans
       if (kans.eq.1) then
        write  (*,*) 'set deltaNNO '
        read (*,*) dnno
        xossi(i)=8.951-24556/(tcent(i)+273.15)
     $           +0.046*(pbar(i)-1)/tkelv(i)+dnno
       endif
       if (kans.eq.2) then
        write  (*,*) 'set deltaFMQ '
        read (*,*) dfmq
        xossi(i)=10.50-26913/(tcent(i)+273.15)+dfmq
       endif
       if (kans.eq.3) then
        write(*,*)'set log fO2 '
        read (*,*) xossi(i)
       endif
       if (kans.ne.1.and.kans.ne.2.and.kans.ne.3) go to 798
       cost=18896.92/(tcent(i)+273.15)-3.82
       fS2(i)=2.*fso2-2.*xossi(i)-2.*cost
       write (*,*) fsO2,fS2(i),xossi(i)
c      pause
      endif
 
      write(6,*)'ho letto ',i
      sulf(i)=stot(i)
      stotT(i)=stot(i)/(10000*32.064)
      tkelv(i)=tcent(i)+273.15
      FSS2=10.**FS2(I)
      FOO2=10.**XOSSI(I)
      CSTART(I)=STOTT(I)*32.064*DSQRT(FOO2/FSS2)
      cstar(I)=STOTT(I)*32.064/(FOO2**1.5*FSS2**0.5)
      WRITE (*,*) '  '
      WRITE (*,*)'---------------- COMPOSITION ',I,'-------------------'
      WRITE (ISYSPR,*) '  '
      WRITE (ISYSPR,*)'------------- COMPOSITION ',I,'-----------------'
 
      if (poss(i,5).ne.0.or.poss(i,7).ne.0) then
      write (*,*) 'Composition n. ',i,' holds Iron !!!'
      write (*,*) 'Fe2O3 =',poss(i,5)
      write (*,*) 'FeO  =',poss(i,7)
      write (ISYSPR,*) 'Composition n. ',i,' holds Iron !!!'
      write (ISYSPR,*) 'Fe2O3 =',poss(i,5)
      write (ISYSPR,*) 'FeO  =',poss(i,7)
      endif
      jflag=0.
 
c  calcola le proporzioni molari a 100% catione e ossido
c  poi si assume acqua misurata indipendentemente (tipo FT-IR)
c  per cui non entra nella normalizzazione.....
      summa = 0.
      summb = 0.
      do 961 lk=1,12
      summa=poss(i,lk)+summa
961   continue
      do 871 lk=1,12
      poss(i,lk)=poss(i,lk)*100./summa
871   continue
      wt(i)=summa
      summc=0.
      do 9623 lh = 1,12
      poss(i,lh)= poss(i,lh)*(1.-poss(i,13)/100.)
      write (*,*) poss(i,lh)
      summc=summc+poss(i,lh)
9623  continue
c     write (*,*) poss(i,13),summc+poss(i,13)
c     pause
 
      sumfe(i)=poss(i,5)+poss(i,7)
      NIT = 50
      do 686 kj=1,NIT
 
      write (*,*) 'Poss 5 = ',poss(i,5)
      write (*,*) 'Poss 7 = ',poss(i,7)
      write (ISYSPR,*) 'Poss 5 = ',poss(i,5)
      write (ISYSPR,*) 'Poss 7 = ',poss(i,7)
 
      if (poss(i,5).ne.0.or.poss(i,7).ne.0) then
      write (*,*) 'Iteration N. ',kj,' @ composition ',i
      write (ISYSPR,*) 'Iteration N. ',kj,' @ composition ',i
      endif
 
c  Procedura di normalizzazione su FeOtotale (somma degli ossidi wt%)
      if (sumfe(i).ne.0) then
      summb=poss(i,5)+poss(i,7)
      poss(i,5)=poss(i,5)*sumfe(i)/summb
      poss(i,7)=poss(i,7)*sumfe(i)/summb
      endif
 
      write (*,878) poss(i,1),poss(i,2),poss(i,3),poss(i,4),poss(i,5),
     &poss(i,6),poss(i,7),poss(i,8),poss(i,9),poss(i,10),poss(i,11),
     &poss(i,12),poss(i,13),summa,summb
      write (ISYSPR,878) poss(i,1),poss(i,2),poss(i,3),poss(i,4),
     &POSS(I,5),poss(i,6),poss(i,7),poss(i,8),poss(i,9),poss(i,10),
     &POSS(I,11),poss(i,12),poss(i,13),summa,summb
878   format (x,13(f5.2,x),/,x,f6.2,x,f6.2)
c
c     write (*,*)'jwflag = ',jwflag
c     pause
c     if (jwflag.gt.1) go to 519
      d(i)=0
      dd(i)=0
      do 100 j=1,13
c     write(*,*) j,poss(i,j),zcat(j),y(j)
      xcat(j)=poss(i,j)*zcat(j)/y(j)
c     pause
      d(i) = d(i)+ xcat(j)
      dd(i)=dd(i)+poss(i,j)/y(j)
 100  continue
      som(i)=d(i)
      dd(i)=dd(i)+stotT(i)
      xcat0(i)=xcat(13)
c        write (*,*) wmol(i),xcat0(i)
c        wmol(i)=(xcat0(i)/2.-wmol(i))*2. ! moli di OH IR-like
c        write (*,*) wmol(i),xcat0(i)
c        pause
c considero tutto H2O come OH (o H!)
      dd(i)=dd(i)-poss(i,13)/y(13)
      xw(i)=2.*poss(i,13)/y(13)
      dd(i)=dd(i)+xw(i)
      somm(i)=dd(i)
      do 200 j=1,13
      xcat(j)=xcat(j)/d(i)
 200  continue
      xfetot(i)=xcat(5)+xcat(7)
      if (kj.eq.1) xftot=xfetot(i)
      if (xftot.ne.0.) then
       write(*,*) 'xfetot-xftot = ',xfetot(i)-xftot
       write(15,*) 'xfetot-xftot = ',xfetot(i)-xftot
      endif
 
      stot(i)=stott(i)/dd(i)
 
CCC   if (xcat(5).ne.0.and.xcat(7).ne.0.and.kj.eq.1) jflag=2
      IF (KFLAG(I).EQ.2) JFLAG=2
 
      if (xcat(5).ne.0.or.xcat(7).ne.0) then
      write (*,*) 'Fe3 = ',xcat(5),' Fe2 = ',xcat(7)
      write (*,*) 'Fetot = ',xftot,' newFetot = ',xfetot(i)
      write (ISYSPR,*) 'Fe3 = ',xcat(5),' Fe2 = ',xcat(7)
      write (ISYSPR,*) 'Fetot = ',xftot,' newFetot = ',xfetot(i)
      endif
       if (xcat(5).ne.0.and.xcat(7).ne.0.) then
        write (*,*) 'redoz is = ',dlog10(xcat(7)/xcat(5))
        write (ISYSPR,*) 'redoz is = ',dlog10(xcat(7)/xcat(5))
        redoz(i)=dlog10(xcat(7)/xcat(5))
       endif
       if (xcat(5).ne.0.and.xcat(7).eq.0.) then
        write (*,*) 'redoz is UNDEFINED'
        write (ISYSPR,*) 'redoz is UNDEFINED'
        redoz(i)=1.0
       endif
       if (xcat(5).eq.0.and.xcat(7).ne.0.) then
        write (*,*) 'redoz is UNDEFINED'
        write (ISYSPR,*) 'redoz is UNDEFINED'
        redoz(i)=1.00
       endif
       if (xfetot(i).eq.0) then
        write (*,*) 'redoz is UNDEFINED'
        write (ISYSPR,*) 'redoz is UNDEFINED'
        xcat(5)=0.0000000001
        xcat(7)=0.0000000001
       endif
 
       if (kj.eq.1) then
        refirst(i)=redoz(i)
       endif
 
519   continue
c
c  Verifica la consistenza.
      FRIT=0.
      KWFLAG=0
      jwflag=0
      ITER=50
      if (poss(I,13).eq.0.) then
       iter=1
       kwflag=3
      endif
      do 789 jk=1,iter
      jwflag=jwflag+1
 
      CALL CALCOMP
      CALL BASOPT
      CALL TOOPSAMIS
      CALL ACTION
 
      frit=dabs(aossi(i)-aossiz(i))
        if (poss(i,13).ne.0.and.frit.lt.0.00001.and.jk.gt.1) then
c        write (13,*) i,jk,xoh(i),xohz(i),aossi(i),ah(i)
         kwflag = 3
        endif
 
      if (kwflag.eq.3) goto 790
      if (jk.eq.iter) then
      write (isyspr,*) ' STRONZO...SU', i
      pause
      goto 790
 
      endif
 
      write (*,*) 'ESCO DA ',i,jk,' nOH-= ',xoh(i),' nH+= ',xh(i),
     $'ITERATION VALUE = ',frit
c     pause
 
789   continue
 
790   CALL TOOPSAMIS2
      CALL ACTION2
      CALL FEREDOX
 
      if (xfetot(i).eq.0.or.jflag.eq.1.) then
      goto 696
      endif
 
      write (*,*) 'redox last (redoX) = ',redox(i)
      write (*,*) 'redox before (redoZ) = ',redoz(i)
      write (15,*) 'redox last (redoX) = ',redox(i)
      write (15,*) 'redox before (redoZ) = ',redoz(i)
      write (15,*) 'Fe tot = ',xfetot(i)
 
      zxf3=(xfetot(i)/(10**redoz(i)+1))
      zxf2=(xfetot(i)/(10**redoz(i)+1))*(10**redoz(i))
 
      xf3=(xfetot(i)/(10**redox(i)+1))
      xf2=(xfetot(i)/(10**redox(i)+1))*(10**redox(i))
 
c      xff3=xf3*d(i)
c      xff2=xf2*d(i)
      xff3=xf3*d(i)
      xff2=xf2*d(i)
      pfe3=xff3*y(5)/zcat(5)
      pfe2=xff2*y(7)
 
      if (poss(i,5).eq.0.) then
       poss(i,5)=0.01
       poss(i,7)=poss(i,7)-0.01
      endif
 
      if (poss(i,7).eq.0.) then
       poss(i,7)=0.01
       poss(i,5)=poss(i,5)-0.01
      endif
c    l'ho spostato sotto!!!
c     poss(i,5)=dsqrt(pfe3*poss(i,5))
c     poss(i,7)=dsqrt(pfe2*poss(i,7))
c
      write (*,*) 'Fe2/Fetot calc  = ', xf2/(xf2+xf3)
      write (*,*) 'Fe2/Fetot before= ', zxf2/(zxf2+zxf3)
      write (15,*) 'Fe2/Fetot calc  = ', xf2/(xf2+xf3)
      write (15,*) 'Fe2/Fetot before= ', zxf2/(zxf2+zxf3)
 
      chi=dabs(bossi(i)-xossi(i))
      chio=dabs(xf2/(xf2+xf3)-zxf2/(zxf2+zxf3))
      write (*,*) 'bossi = ',bossi(i), 'xossi = ',xossi(i)
      write (*,*) 'REMINDER: bossi Š ricalcolato da redox!'
      write (*,*) 'REMINDER: bossi = xossi se INPUT logfO2 = 0'
      write (*,*) ' JFLAG is ',jflag
      write (15,*) 'bossi = ',bossi(i), 'xossi = ',xossi(i)
      write (15,*) 'REMINDER: bossi Š ricalcolato da redox!'
      write (15,*) 'REMINDER: bossi = xossi se INPUT logfO2 = 0'
      write (15,*) ' JFLAG is ',jflag
c     write (15,*) 'NOW Fe2+ is '
c     write (15,*) 'NOW Fe3+ is '
c     pause
      if (jflag.eq.2) then
      write (*,*)' !!! I HAVE USED DEFAULT FeO/Fe2O3 AND INPUT fO2 !!!'
      write (*,*)' !!!           I MADE NO ITERATIONS              !!!'
c     pause
      go to 696
      endif
 
c L'ho spostato qua ma stava sopra!!!!!
 
      poss(i,5)=dsqrt(pfe3*poss(i,5))
      poss(i,7)=dsqrt(pfe2*poss(i,7))
 
      if (chi.lt.0.0001.and.chio.lt.0.0001) then
      jflag=1
      write (*,*) 'JFLAG= ',jflag
      write (15,*) 'JFLAG= ',jflag
c     pause
      endif
 
      if (kj.eq.nit) then
      write (*,*) 'This was the last iteration. '
      write (15,*) 'This was the last iteration. '
c     pause
      endif
c     if (i.lt.ncomp) pause
686   continue
 
696   continue
 
C     rapmed=dlog10(cstart(i))+rapmed
 
*********************************************************************
      write(ISYSIN,888) poss(i,1),poss(i,2),poss(i,4),poss(i,5),
     &POSS(I,6),poss(i,7),poss(i,8),poss(i,10),
     &poss(i,9),poss(i,11),poss(i,12),poss(i,3), poss(i,13)
888   format(13(F9.5,1x))
 
      CALL IONPOR
 
c  calcola le proporzioni molari a 100% catione e ossido
      summit = 0.
      do 991 lk=1,13
       summit=poss(i,lk)/y(lk)+summit
991   continue
      xh2o(i)=(poss(i,13)/y(13))/summit
***************************************************************************
 
      CALL CALCPROP
      CALL DATAOUT
111   CONTINUE
      END
 
      SUBROUTINE CALCOMP
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes stoichiometric amounts                             *
*     charge balance for melt complexes                           *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      common/setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
c     if (poss(i,13).ne.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NOTA BENE............                                                      c
c INIZIA DA QUI IL RICALCOLO DI H2O SPECIATA SECONDO MORETTI & OTTONELLO     c
c (2003) IN PREP. TRA H+ E OH-. OH- Š ANIONE LIBERO SULLA MATRICE ANIONICA   c
c MENTRE H+ CONTRIBUISCE ALLA MATRICE W CON QUANTO NE CONSEGUE.              c
c QUINDI H+ ENTRA NEL COMPUTO DELLA DEPOLIMERIZZAZIONE SECONDO QUANTO EVINTO c
c DA FRASER.I RUOLI STRUTTURALI NON CI INTERESSANO (CHIUSURA TERMINAZIONI    c
c POLIMERICHE O QUANT'ALTRO....)                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      xcat(1)=poss(i,1)*zcat(1)/y(1)
      xcat(3)=poss(i,3)*zcat(3)/y(3)
      xcat(4)=poss(i,4)*zcat(4)/y(4)
      xcat(5)=poss(i,5)*zcat(5)/y(5)
      xcat(6)=poss(i,6)*zcat(6)/y(6)
      xcat(7)=poss(i,7)*zcat(7)/y(7)
      xcat(2)=poss(i,2)*zcat(2)/y(2)
      xcat(13)=poss(i,13)*zcat(13)/y(13)
      xcat(8)=poss(i,8)*zcat(8)/y(8)
      xcat(9)=poss(i,9)*zcat(9)/y(9)
      xcat(10)=poss(i,10)*zcat(10)/y(10)
      xcat(11)=poss(i,11)*zcat(11)/y(11)
      xcat(12)=poss(i,12)*zcat(12)/y(12)
      xcat(1)=xcat(1)/som(i)
      xcat(3)=xcat(3)/som(i)
      xcat(4)=xcat(4)/som(i)
      xcat(5)=xcat(5)/som(i)
      xcat(6)=xcat(6)/som(i)
      xcat(2)=xcat(2)/som(i)
      xcat(7)=xcat(7)/som(i)
      xcat(8)=xcat(8)/som(i)
      xcat(9)=xcat(9)/som(i)
      xcat(10)=xcat(10)/som(i)
      xcat(11)=xcat(11)/som(i)
      xcat(12)=xcat(12)/som(i)
      xcat(13)=xcat(13)/som(i)
      aossiz(i)=aossi(i)
      if (jk.eq.1.and.iter.gt.1.) then
      aossiz(i)=1.
      xoh(i)=xcat0(i)
      endif
 
      xh(i)=xcat0(i)-xoh(i)
      xohz(i)=xoh(i)
      xcat(13)=xh(i)/som(i)
      xwd(i)=xoh(i)/som(i)
 
c     endif
 
       w=3.0*xcat(6)+2.0*(xcat(7)+xcat(8)+xcat(9)+xcat(10))+xcat(11)+
     & xcat(12)+xcat(13)+xcat(3)+4.0*xcat(2)
 
      WRITE(*,*) ' W (charges due to modifiers) = ',W
      WRITE(ISYSPR,*) ' W (charges due to modifiers) = ',W
 
      IF (Xcat(4).GT.W) GO TO 540
      xfor=xcat(1)+xcat(4)/2.0
      w=w-xcat(4)
      xallu=xcat(4)/2.0
      xcat(4)=0.
      GO TO 590
 540  u=xcat(4)-w
      xfor=xcat(1)+w/2.0
      xcat(4)=u
      XALLU=W/2.0
      w=0.
      GO TO 680
 590  IF (xcat(5).GT.W) GO TO 640
      xfor=xfor+xcat(5)/2.0
      w=w-xcat(5)
      xfe2o3=xcat(5)/2.0
c      feo2m=xcat(5)
      afeo2(i)=xcat(5)/2.
      afeo22(i)=xcat(5)/2.0
c      fe3p=0
      xcat(5)=0.
      GO TO 680
 640  v=xcat(5)-w
      xfor=xfor+w/2.0
      afeo2(i)=w/2.
      afeo22(i)=w/2.0
      xcat(5)=v
      xfe2o3=w/2.0
c      feo2m=w
c      fe3p=v
      w=0.
c     pause
 680  xfor=xfor+xcat(3)/2.0
c      afeo2(i)=afeo2(i)/xfor
 1691 CONTINUE
 
      END
 
      SUBROUTINE BASOPT
 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes optical basicity of network formers and network          *
*     modifiers;                                                        *
*     defines oxides stoichiometric proportions for oxide-sulfide       *
*     reactions                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 k,mnom,mgom,nao05m,ko05m
      PARAMETER (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
C
      basfor=xcat(1)*z(1)/xfor+xcat(3)*z(3)/xfor+xallu*z(4)/xfor+
     & xfe2o3*z(5)/xfor
c
      scat(1)=xcat(1)
      scat(3)=xcat(3)
      scat(2)=0.0
      scat(4)=xallu
      scat(5)=xfe2o3
      write (isysck,*)
      write (isysck,*) 'BASOPT Composition: ',i
      write (isysck,*) 'BASOPT Iteration; ',kj
      write (isysck,*) 'BASOPT XFe2O3:', xfe2o3
      write (isysck,*) 'BASOPT Xfor: ', xfor
      write (isysck,*) 'BASOPT Scat(5): ',scat(5)
      write (isysck,*) 'BASOPT Redox: ' ,redox(i)
      write (isysck,*) 'BASOPT Redoz: ',redoz(i)
      write (isysck,*)
      write (isyspr,*) ' + + + INCIDENTAL S CAPACITIES LOG BLOCK + + +'
      write (isyspr,*) 'F&G elect. equiv. frac. per Fe3+ = ',scat(5)
      xmod=0.0
      do 6969 ij=2,13
      xmod=xmod+xcat(ij)/zcat(ij)
6969  continue
      basmod=0.0
      do 6996 ij=2,13
      xcat(ij)=xcat(ij)/(zcat(ij)*xmod)
      basmod=basmod+xcat(ij)*z(ij)
6996  continue
      write (isyspr,*) 'F&G elect. equiv. frac. per Fe2+ = ',xcat(7)
      write (isyspr,*) ' + + + INCIDENTAL S CAPACITIES LOG BLOCK + + +'
      afe2(i)=xcat(7)
      afe22(i)=xcat(7)
c     afe2(i)=xcat(7)/xcat(10) ! usato solo per fornire [Fe]/[Mg]
      afe3(i)=xcat(5)
      afe32(i)=xcat(5)
      ana(i)=xcat(11)
      ak(i)=xcat(12)
      amg(i)=xcat(10)
      amn(i)=xcat(8)
      aca(i)=xcat(9)
      acr(i)=xcat(6)
      ati(i)=xcat(2)
      aprot(i)=xcat(13)
      ap(i)=xcat(3)
      k=dexp(((basmod-basfor)/0.2145)-1.1445)
      polcos(i)=k
      basdif(i)=basmod-basfor
      a=1.0-4.0*k
      xcat(1)=xfor/(xfor+xmod)
      acidic(i)=xcat(1)
      e=1.0-xcat(1)
      totcat(i)=e
      afe22(i)=afe22(i)/(xfor+xmod)
      afe32(i)=afe32(i)/(xfor+xmod)
      afe22(i)=afe22(i)/totcat(i)*xmod
      afe32(i)=afe32(i)/totcat(i)*xmod
      afeo2(i)=afeo2(i)/(xfor+xmod)
      afeo22(i)=afeo22(i)/(xfor+xmod)
      END
 
      SUBROUTINE TOOPSAMIS
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  polymeric aproach to melt structure (see Ottonello et al., 2000; *
*  Ottonello, 2001; for details)                                    *
*  computes singly bonded, doubly bonded and free oxygen ion on a   *
*  1-mole of melt basis and furnishes structural details such as    *
*  mean extension of polymeric units, acidity etc. etc.             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      atoop=-a
      btoop=2.0+2.0*xcat(1)
      ctoop=8.0*xcat(1)*(xcat(1)-1.0)
      o(i)=(-btoop+dsqrt(btoop**2.-4.0*atoop*ctoop))/(2.0*atoop)
 
      o1(i)=1.0-xcat(1)-o(i)/2.0
      o3(i)=(4.0*xcat(1)-o(i))/2.0
      o4(i)=o(i)/(o(i)+o3(i)+xcat(1))
c      write (*,*) ' stong ca!!!!!!!!!!!!!!!!'
c      write (*,*) o1(i)+o3(i)+o(i)
c      pause
      s=dexp(-1.7165*dlog(o4(i))+2.8776)
      s1=acidic(i)/s
      o2=o1(i)/(o1(i)+s1+stot(i)+xwd(i)) ! o2 Š attivit… di O=!
c     aossi(i)=o2
      totani(i)=o1(i)+s1+stot(i)+xwd(i)
c      write (*,*)
c      write (*,*) 'structons = ',s1
c      write (*,*) 'O2meno = ',o1(i)
c      write (*,*) 'S anions = ',stot(i)
c      write (*,*) 'OH- = ',xwd(i)
c      write (*,*) 'total anions = ',totani(i)
c      pause
      totali(i)=totani(i)
      sobs(i)=stot(i)/totani(i)
c     afeo2(i)=afeo2(I)/s
c     aFeO2(i)=afeo2(i)/totani(i)
c     afeo22(i)=afeo22(i)/s
c     afeo22(i)=afeo22(i)/totani(i)
      ah(i)=totcat(i)/totani(i)
      END
 
      SUBROUTINE ACTION
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes the activity of sulphide and of iron ionic species   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
C
c     sobs(i)=stot(i)/totani(i)
      totold(i)=o1(i)+s1
      struct(i)=s1
      aossi(i)=o2
      aozzi(i)=aossi(i)

c  I volumi sotto sono in joule/bar ...fattore 0.1...
c  Raggi ionici di Shannon
c     voh=(8.28528253/3.)*3.14159*(1.6066-0.38+0./10000*
c    $ (tkelv(i)-298.15))**3
c      voh=(4./3.)*3.14159*(1.6066+0./10000*(tkelv(i)-298.15))**3
       voh=(4./3.)*3.14159*(1.40+0./10000*(tkelv(i)-298.15))**3
      voh=voh*0.6022045
c      vo2=(4./3.)*3.14159*(1.6066+0./10000*(Tkelv(i)-298.15))**3
      vo2=(4./3.)*3.14159*(1.40+0./10000*(Tkelv(i)-298.15))**3
      vo2=vo2*0.6022045
c     vh=(4./3.)*3.14159*(-0.38+0./10000*(tkelv(i)-298.15))**3
      vh=(4./3.)*3.14159*(0.+0./10000*(tkelv(i)-298.15))**3
      vh=vh*0.6022045
      if (VH.LT.0) vh=0.
c Calcolo delv assumendo espansione termica = 0.
      delv0=vh+vo2-voh+0./1000.*(tkelv(i)-298.15)
      delv(i)=delv0
      delv0=delv0*0.1/(8.3147*2.303)
      delv0=delv0*(pbar(i)-1)/tkelv(i)
c INSERISCO I LAST VALUES
c      partw=10.**(-6.2274+7947.29/tkelv(i)-delv0)
c      partw=1.18e38
c       partw=10**(3001.6/tkelv(i)-1.7419-delv0)
c      partw=10**(3696.9/tkelv(i)-3.147-delv0)
c     partw=10**(4679.34/tkelv(i)-3.81613-delv0)
c      partw=10**(7727.019164/tkelv(i)-5.969157)
c      partw=10**(607.514123515489/tkelv(i)-1.05834853820627)
       partwa=-2.89964871774198                    
       partwb=1371.16508083715
       sann=0.218216405605124d-3
       partw=10.**(partwa+partwb/tkelv(i)+0.5*sann*(pbar(i)-1.))
c      partw=10.**(38.)
      ratius=(aossi(i)/ah(i))/partw
      xoh(i)=xcat0(i)*ratius/(ratius+1)
      xh(I)=XCAT0(I)-XOH(I)
      xoh(i)=dsqrt(xohz(i)*xoh(i))
 
      END
 
      SUBROUTINE TOOPSAMIS2
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  polymeric aproach to melt structure (see Ottonello et al., 2000; *
*  Ottonello, 2001; for details)                                    *
*  computes singly bonded, doubly bonded and free oxygen ion on a   *
*  1-mole of melt basis and furnishes structural details such as    *
*  mean extension of polymeric units, acidity etc. etc.             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
c speculazione su H+ e O-
c     voh=(4./3.)*3.14159*((1.6066-0.38+0./1.d6*(tkelv(i)-298.15))**3)
c      voh=(4./3.)*3.14159*((1.6066+0./1.d6*(tkelv(i)-298.15))**3)
      voh=(4./3.)*3.14159*((1.40+0./1.d6*(tkelv(i)-298.15))**3)
      voh=voh*0.6022045
c  sarebbe x(7) il coeff lineare di espansione
c      vo=(4./3.)*3.14159*((1.6066+0./1.d6*(Tkelv(i)-298.15))**3)
      vo=(4./3.)*3.14159*((1.40+0./1.d6*(Tkelv(i)-298.15))**3)
      vo=vo*0.6022045
c     vh=(4./3.)*3.14159*((-0.38+0./1.d6*(tkelv(i)-298.15))**3)
      vh=(4./3.)*3.14159*((0.+0./1.d6*(tkelv(i)-298.15))**3)
      vh=vh*0.6022045
      if (VH.LT.0) vh=0.
      dv1=-vh-vo+voh
      dv(i)=dv1
      dv1=dv1*0.1/(8.3147*2.303)
      dv1=dv1*(pbar(i)-1)/tkelv(i)
c      costx=10.**(-0.952+0./tkelv(i)-dv1)
c      costx=1.18e-38
c     costx=10**(-0.99855)
c      costx=10**(-0.9857)
c      costx=10**(-0.99917503)
c      costx=10**(-132)
      aoin=o(i)
      ahin=xh(i)/som(i)
c      aspec=costx
c      bspec=-(costx*ahin+costx*aoin+totcat(i))
c      cspec=costx*ahin*aoin
c      write (*,*) som(i), 'MERDA'
c      root1=(-bspec-dsqrt(bspec**2-4*aspec*cspec))/(2*aspec)
c     root2=(-bspec+dsqrt(bspec**2-4*aspec*cspec))/(2*aspec)
c ATTENZIONE EFFETTUO UN CAMBIO IMPORTANTE:
      root1=0.
      o(i)=o(i)-root1
 
      o1(i)=1.0-xcat(1)-o(i)/2.0
      o3(i)=(4.0*xcat(1)-o(i))/2.0
      o4(i)=o(i)/(o(i)+o3(i)+xcat(1)+root1)
 
      root(i)=root1
 
      s=dexp(-1.7165*dlog(o4(i))+2.8776)
      sii(i)=s
      s1=acidic(i)/s
      o2=o1(i)/(o1(i)+s1+stot(i)+xwd(i))
c     aossi(i)=o2
      totani(i)=o1(i)+s1+stot(i)+xwd(i)
c      afeo2(i)=afeo2(I)/s
      aFeO2(i)=afeo2(i)/totani(i)/xcat(1)
      afeo22(i)=afeo22(i)/s
      afeo22(i)=afeo22(i)/totani(i)
c     write (*,*) 's = ',s,' afeo22 = ',afeo22(i)
      ah(i)=totcat(i)/totani(i)
      END
 
      SUBROUTINE ACTION2
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes the activity of sulphide and of iron ionic species   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
C
c      sobs(i)=stot(i)/totani(i)
      totold(i)=o1(i)+s1
      struct(i)=s1
      aossi(i)=o2
 
      END
 
      SUBROUTINE FEREDOX
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* computes the redox state of iron in terms of oxides proportions *
* and/or the ensuying oxygen partial pressure at equilibrium      *
* (utilized whenever the data are unsufficient or not precise)    *
c the adopted constants are from Ottonello et al. (2000)          *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
c Pressure effect on cost 54 (cost 21 in published work) accounting
c We assume that for the other constants -equilibrium among components
c in the liquid phase- the Dvolume of reaction is 0.
c calcola la fo2 in base al rapporto di ossidazione del ferro
 
c  cost54 => FeO + 0.25O2 <=> FeO1.5
c  cost57 => Fe2O3 + O2- <=> 2FeO2-
c  cost58 => Fe2O3 <=> 2Fe3+ + 3O2-
c  cost59 => FeO <=> Fe2+ + O2-
 
********* OLD VALUES (Chemical Geology) **************
 
      cost54=-2.8792+6364.8/Tkelv(i)
      cost57=2.913128*2.-1723.12*2./Tkelv(i)
      cost58=0.856683*2.-2033.56*2./Tkelv(i)
      cost59=1.1529-1622.4/Tkelv(i)
 
c      cos54=10.**cost54
c      cos57=10.**cost57
c      cos58=10.**cost58
c      cos59=10.**cost59

********* NEW VALUES Moretti and Ottonello (2003 in prep.) *******
 
c     cost54=-2.8792+6364.8/Tkelv(i)
c     cost57=5.2588-3023.2/Tkelv(i)
c     cost58=1.9269-4203.3/Tkelv(i)
c     cost59=1.1529-1622.4/Tkelv(i)
 
********* NEWNEW VALUES 04/08/03 Moretti and Ottonello (2003 in prep.) *******
 
c     cost54=-2.8792+6364.8/Tkelv(i)
c     cost57=5.2878-3358.9/Tkelv(i)
c     cost58=1.9221-4103.9/Tkelv(i)
c     cost59=1.1529-1622.4/Tkelv(i)
 
ccc   Now we use cost54 at pressures greater than 1 bar. FeO and FeO1.5
ccc   activities will be corrected for volume terms.
ccc   We use the values given by Lange in Carroll and Webster, 1994
ccc   They have also been recomputed with correlations with optical basicity.
ccc   First, we must recalculate for every T at 1 bar pressure
ccc   Values are in cc/mole.
ccc   These are parameters recommended for metaluminous and peralkaline
ccc   liquid !!!!!!  But our extrapolation makes them good for every melt.
 
c   VFeO2 usando i CR (assumo struttura vincolata e non free anions)
c   0.63 Š HS CR di Shannon per Fe3+ (IR sarebbe 0.49) in coord IV
c   0.55 Š LS IR di Shannon per Fe3+ (HS sarebbe 0.645) in coord VI
c   0.61 Š LS IR di Shannon per Fe2+ (HS sarebbe 0.780) in coord VI
c   1.24 Š CR di Shannon per O2- (IR sarebbe 1.38) in coord IV
c   1.26 Š CR di Shannon per O2- (IR Š 1.40...) in coord VI
c     vfeo2m=4./3.*3.14159*((0.63**3)+2*(1.24**3))
c     vfeo2m=4./3.*3.14159*((2*1.24+0.63)**3)
C     vfeo2m=4./3.*3.14159*(3.25+0./10000*(Tkelv(i)-298.15))**3
C SOTTO LAST VAlues
c     vfeo2m=4./3.*3.14159*(0.49+2*1.6066+0./10000*(Tkelv(i)-298.15))**3
      vfeo2m=4./3.*3.14159*(0.49+2*1.40+
     $ 0./10000*(Tkelv(i)-298.15))**3
c      vo2m=4./3.*3.14159*((1.6066+0./10000*(Tkelv(i)-298.15))**3)
      vo2m=4./3.*3.14159*((1.40+0./10000*(Tkelv(i)-298.15))**3)
      vfe3p=4./3.*3.14159*((0.645+0./10000*
     $ (Tkelv(i)-298.15))**3)
      vfe2p=4./3.*3.14159*((0.78+0./10000*(Tkelv(i)-298.15))**3)
      conv=0.6022045
      vfeo2m=vfeo2m*conv
      vfeo2m=vfeo2m+0.*0.001*(Tkelv(i)-298.15)
      vo2m=vo2m*conv
      vo2m=vo2m+0.*0.001*(Tkelv(i)-298.15)
c Di seguito equivale ad assumere che DV di FeO2- + 2O2- = FeO45- Š 0.
      vfeo2m=vfeo2m-2*vo2m
      vfe3p=vfe3p*conv
      vfe3p=vfe3p+0.*0.001*(Tkelv(i)-298.15)
      vfe2p=vfe2p*conv
      vfe2p=vfe2p+0.*0.001*(Tkelv(i)-298.15)
c     write (*,*) 'vfeo2m   vo2m    vfe3p    vfe2p'
c     write (*,*) vfeo2m,vo2m,vfe3p,vfe2p
c     pause
c sopra ho ripetuto: x(8)=x(9). CioŠ le espansioni di Fe3 e Fe2 sono uguali
c espansione libera: sui volumi molari espressione intera
c     vfeO15=21.065+4.545*0.001*(Tkelv(i)-1673)
c     vfeO=13.65+2.92*0.001*(Tkelv(i)-1673)
c espansione fissata: T = 298.15
      vfeO15=21.065+4.545*0.001*(298.15-1673)
      vfeO=13.65+2.92*0.001*(298.15-1673)
 
      dv57=2*(vfeo2m-vfeo15-0.5*vo2m)
      dv58=2*(vfe3p+1.5*vo2m-vfeo15)
      dv59=vfe2p+vo2m-vfeo
c da scommentare per verificare consistenza con gaslast
c     dv57=0
c     dv58=0
c     dv59=0

c     da Fe3/Fe2 GCA
c      sann=0.303739519291106d-3     ! good value
      sann=0.218216405605124d-3
c      sann=0.
      cos57a=cost57
      cost57=cost57-(dv57*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
      cost57=cost57-sann*(pbar(i)-1.)
      cos58a=cost58
      cost58=cost58-(dv58*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
      cost58=cost58-1.5*sann*(pbar(i)-1.)
      cos59a=cost59
      cost59=cost59-(dv59*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
      cost59=cost59-sann*(pbar(i)-1.)
c Therefore:
      cost54=10.**cost54
      cost57=10.**cost57
      cost58=10.**cost58
      cost59=10.**cost59
 
      cos57a=10.**cos57a
      cos58a=10.**cos58a
      cos59a=10.**cos59a

*********************************************************
 
      V0II1 =voxm(7)+eoxm(7)*0.001*(tkelv(i)-1673.)
      V0III1=voxm(5)+eoxm(5)*0.001*(tkelv(i)-1673.)
c     v0iip =v0ii1 +(pbar(i)-1.)*0.0001*(-1.8)
c     v0iiip=(v0iii1+(pbar(i)-1.)*0.0001*3.1)*0.5
 
c     terms above are converted in joule/bar
      v0ii1 =0.1*v0ii1
      v0iii1=0.1*v0iii1
 
      cFeii =coxm(7)*0.0001*0.1
      cfeiii=coxm(5)*0.0001*0.1
 
      dvs = (0.5*v0iii1-v0ii1)*(pbar(i)-1.)+(0.5*cfeiii-cfeii)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   Dvs becomes dvs/RT
      dvs=dvs/(tkelv(i)*8.3147)
c   Correction term dvs entered below...
 
c serve solo quando non sono dati o sono imprecisi
c qui fissiamo FeII/FeIII e calcoliamo fO2
c partiamo da redoZ
c     rappox(i)=dlog10(xcat(7)/xcat(5))
      raptru=10.0**redoz(i)
      den1=cost54*raptru
      den1=den1/(dexp(dvs))
      den2=cost57**(0.5)*aossi(i)**2*totani(i)+cost58**(0.5)*totcat(i)
      razio=aossi(i)**(0.5)*cost59*totcat(i)/(den1*den2)
      bossi(i)=razio**4.
      bossi(i)=dlog10(bossi(i))
 
***  Vale nel caso che fO2 non sia data in input e si ha Fe2+/Fe3+ ***
      if (xossi(i).eq.0.and.xftot.gt.0.) xossi(i)=bossi(i)
**************************************************************************
c     Per calcolo sulfuro
 
      beppe(i)=(10.**(xossi(i)-fs2(i)))
c     sppm(i)=aossi(i)/(beppe(i)**0.5)
      sppm(i)=aozzi(i)/(beppe(i)**0.5)
 
c     Per calcolo solfato
 
      beppi(i)=((10.**xossi(i))**1.5)*((10.**fs2(i))**0.5)
      sppb(i)=aossi(i)*beppi(i)
      sppb(i)=aozzi(i)*beppi(i)
 
      XOSS =10.0**XOSSI(I)
      den1=cost54*xoss**0.25
      den1a=cost54*(10**bossi(i))**0.25
      den11=den1

c sotto abbiamo fFeO1.5/fFeO se FeII/FeIII sono fissati costanti
      fa3fa2(i)=den11
c sotto abbiamo fFeO1.5/fFeO se fO2 Š fissata costante
      fa3f2a(i)=den1a
      fcox(i)=0.
      fcoz(i)=0.

      write (isysck,*) cost54,xoss
      den1=den1/(dexp(dvs))
      den1a=den1a/(dexp(dvs))
      den2=cost57**(0.5)*aossi(i)**2*totani(i)+cost58**(0.5)*totcat(i)
      ratio=aossi(i)**(0.5)*cost59*totcat(i)/(den1*den2)
c     write (isysck,*) aossi(i),totcat(i),cost59,den1,den2
      a3a2(i)=den1
      a3a2z(i)=den1a
      ca3ca2(i)=aossi(i)**(0.5)*cost59*totcat(i)/den2
      ca3c2a(i)=aossi(i)**(0.5)*cost59*totcat(i)/den2
      afeo2m(i)=0.
      ffo2m(i)=0.
      ffeo2m(i)=0.
      afe3m(i)=0.
      ffe3m(i)=0.
      ffe2(i)=0.

c  questa sotto Š una mezza intoppata....
      pippa=fa3f2a(i)*(10**redoz(i))

      if (xfetot(i).ne.0.) then

c  sotto abbiamo il rapporto [FeO2-]/[Fe2+] come log
      afeo2m(i)=((cost57**0.5)/cost59)*(aossi(i)**1.5)*den1 !!! OK lineare!!!
      afeo2m(i)=dlog10(afeo2m(i))                           !!! log !!!

c  iniziamo con correggereaO2-  in fO2-
      ffo2m(i)=dlog10(aossi(i))+(1./(8.3147*2.303*tkelv(i)))
     $ *(pbar(i)-1)*vo2m*0.1   !log!
      ffo2m(i)=10.**ffo2m(i)   !lineare!

c  sotto abbiamo il rapporto f(FeO2-)/f(Fe2+) come log
      ffeo2m(i)=den11*(ffo2m(i)**1.5)*((cost57**0.5)/cost59) !lineare!
      ffeo2m(i)=dlog10(ffeo2m(i))                             !log!

c  sotto abbiamo il rapporto [Fe3+/Fe2+] come log
      afe3m(i)=((cost58**0.5)/cost59)*den1/(aossi(i)**0.5) !lineare!
      afe3m(i)=dlog10(afe3m(i))                            !log!

c  sotto abbiamo il rapporto [Fe3+/Fe2+] come log
      ffe3m(i)=((cost58**0.5)/cost59)*den11/(ffo2m(i)**0.5) !lineare!
      ffe3m(i)=dlog10(ffe3m(i))                           !log!

c sotto calcolo fFe2+
      ffe2(i)=dlog10(afe2(i))+(1./(8.3147*2.303*tkelv(i)))
     $ *(pbar(i)-1)*vfe2p*0.1
      ffe2(i)=10.**ffe2(i)

      endif

c N.B.: fcox e fcoz non coincidono!!!! In realt… non possono per via dei
c volumi!

      den2a=cos57a**(0.5)*ffo2m(i)**2*totani(i)+cos58a**(0.5)*totcat(i)
      fcox(i)=ffo2m(i)**(0.5)*cos59a*totcat(i)/den2a
      fcoz(i)=pippa

c Sotto Š l'equazione da risolvere per calcolare la fug. di O2-
c      pippa*(cos57a**(0.5)*ffo2m(i)**2*totani(i)+cos58a**(0.5)*totcat(i))
c      =ffo2m(i)**(0.5)*cos59a*totcat(i)/den2a
c------------------------------------------------------
      redox(i)=dlog10(ratio)
c      write (*,*) 'superstronzo!!!!!'
c      pause
 
      write (isysck,*)
      write (isysck,*) 'FEREDOX Composition: ',i
      write (isysck,*) 'FEREDOX Iteration; ',kj
      write (isysck,*) 'FEREDOX XFe2O3:', xfe2o3
      write (isysck,*) 'FEREDOX Xfor: ', xfor
      write (isysck,*) 'FEREDOX Scat(5): ',scat(5)
      write (isysck,*) 'FEREDOX Redox: ' ,redox(i)
      write (isysck,*) 'FEREDOX Redoz: ',redoz(i)
      write (isysck,*)
 
      END
 
      SUBROUTINE IONPOR
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes ionic porosity                                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      common/setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      dq=0
      ddq=0
      xoan=0.
      vsio2=0
      vtio2=0
      vp2o5=0
      val2o3=0
      vfe2O3=0
      vCr2o3=0
      vfeo=0
      vmno=0
      vmgo=0
      vcao=0
      vna2o=0
      vk2o=0
      vh2o=0
      vcsi=0
      vcti=0
      vcp=0
      vcal=0
      vcfe=0
      vcCr=0
      vcfe=0
      vcmn=0
      vcmg=0
      vcca=0
      vcna=0
      vck=0
      vch=0
      vao=0
      sumol=0
 
      do 110 j=1,13
      xxcat(i,j)=poss(i,j)*zcat(j)/y(j)
      xmol(i,j)=poss(i,j)/y(j)
      dq=dq+xxcat(i,j)
      ddq=ddq+poss(i,j)/y(j)
      xoan=xmol(i,j)*noxy(j)+xoan
 110  continue
c     xxcat(i,13)=xoan              !!! okkio a non confondersi...
c     xoan=xoan-xmol(i,13)*noxy(13)
      sumol=dq+xoan
      do 210 j=1,13
      xmol(i,j)=xmol(i,j)/ddq
      xxcat(i,j)=xxcat(i,j)/sumol
 210  continue
c      xmol0(i)=xmol(i,13)
 
c     ora calcola i volumi molari e ionico per il calcolo di IP
 
      vsio2=voxm(1)+eoxm(1)*0.001*(tkelv(i)-1673.)
      vtio2=voxm(2)+eoxm(2)*0.001*(tkelv(i)-1673.)
      vp2o5=voxm(3)+eoxm(3)*0.001*(tkelv(i)-1673.)
      val2o3=voxm(4)+eoxm(4)*0.001*(tkelv(i)-1673.)
      vfe2o3=voxm(5)+eoxm(5)*0.001*(tkelv(i)-1673.)
      vcr2o3=voxm(6)+eoxm(6)*0.001*(tkelv(i)-1673.)
      vfeO=voxm(7)+eoxm(7)*0.001*(tkelv(i)-1673.)
      vmnO=voxm(8)+eoxm(8)*0.001*(tkelv(i)-1673.)
      vcaO=voxm(9)+eoxm(9)*0.001*(tkelv(i)-1673.)
      vmgO=voxm(10)+eoxm(10)*0.001*(tkelv(i)-1673.)
      vna2O=voxm(11)+eoxm(11)*0.001*(tkelv(i)-1673.)
      vk2O=voxm(12)+eoxm(12)*0.001*(tkelv(i)-1673.)
      vh2O=voxm(13)+eoxm(13)*0.001*(tkelv(i)-1673.)
 
      vsio2=vsio2*xmol(i,1)
      vtio2=vtio2*xmol(i,2)
      vp2o5=vp2o5*xmol(i,3)
      val2o3=val2o3*xmol(i,4)
      vfe2o3=vfe2o3*xmol(i,5)
      vcr2o3=vcr2o3*xmol(i,6)
      vfeo=vfeo*xmol(i,7)
      vmno=vmno*xmol(i,8)
      vcao=vcao*xmol(i,9)
      vmgo=vmgo*xmol(i,10)
      vna2o=vna2o*xmol(i,11)
      vk2o=vk2o*xmol(i,12)
      vh2o=vh2o*xmol(i,13)
 
      volm(i)=vsio2+vtio2+vp2o5+val2o3+vfe2o3+vcr2o3+vfeo+vmno+
     $vcao+vmgo+vna2o+vk2o+vh2o
 
      fact=4./3.*3.14159265358979*0.6022045
      radio=1.40
 
      vcsi=fact*radi(1)**3
      vcti=fact*radi(2)**3
      vcp=fact*radi(3)**3
      vcal=fact*radi(4)**3
      vcfe=fact*radi(5)**3
      vccr=fact*radi(6)**3
      vcfe=fact*radi(7)**3
      vcmn=fact*radi(8)**3
      vcca=fact*radi(9)**3
      vcmg=fact*radi(10)**3
      vcna=fact*radi(11)**3
      vck=fact*radi(12)**3
      vch=fact*radi(13)**3
      vao=fact*radio**3
 
      vcsi=vcsi*xxcat(i,1)
      vcti=vcti*xxcat(i,2)
      vccp=vccp*xxcat(i,3)
      vcal=vcal*xxcat(i,4)
      vcfe=vcfe*xxcat(i,5)
      vccr=vccr*xxcat(i,6)
      vcfe=vcfe*xxcat(i,7)
      vcmn=vcmn*xxcat(i,8)
      vcca=vcca*xxcat(i,9)
      vcmg=vcmg*xxcat(i,10)
      vcna=vcna*xxcat(i,11)
      vck=vck*xxcat(i,12)
      vch=vch*xxcat(i,13)
      vao=vao*xxcat(1,13)
 
      volion(i)=vcsi+vcti+vcp+vcal+vcfe+vccr+vcfe+vcmn+vcca+vcmg
     $+vcna+vck+vch+vao
 
      por(i)=100*(1-volion(i)/volm(i))
c     write (*,*) i,volion(i),volm(i),por(i)
c     pause
      END
 
 
 
      SUBROUTINE CALCPROP
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* assigns network-former and network-modifier fractions for     *
* sulfidation reactions of type                                 *
*   SiO2 + 0.5S2(g)  <=>  SiOS + 0.5O2(g)                       *
*   TiO2 +0.5S2(g)   <=>  TiOS + 0.5O2(g)                       *
*   P2O5 +0.5S2(g)   <=>  P2O4S + 0.5O2(g)                      *
*   Al2O3 +0.5S2(g)  <=>  Al2O2S + 0.5O2(g)                     *
*   Fe2O3 +0.5S2(g)  <=>  Fe2O2S + 0.5O2(g)                     *
*   Cr2O3 +0.5S2(g)  <=>  Cr2O2S + 0.5O2(g)                     *
*   FeO + 0.5S2(g)   <=>  FeS + 0.5O2(g)                        *
*   MnO + 0.5S2(g)   <=>  MnS + 0.5O2(g)                        *
*   CaO + 0.5S2(g)   <=>  CaS + 0.5O2(g)                        *
*   MgO + 0.5S2(g)   <=>  MgS + 0.5O2(g)                        *
*   Na2O +0.5S2(g)   <=>  Na2S + 0.5O2(g)                       *
*   K2O + 0.5S2(g)   <=>  K2S +  0.5O2(g)                       *
*   H2O + 0.5S2(g)   <=>  H2S +  0.5O2(g)                       *
*   values are normalized to the total number of cations in one *
*   mole of melt                                                *
*                                                               *
*   We HAVE ALSO ANALOGOUS SULPHATATION REACTIONS               *
* * * * * * * * ** * * * * * * * * * ** * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,  ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      write (isyspr,*)
      write (isyspr,*)'ENTERING CALCPROP: BEGINS THE SULPHUR DOMAIN'
      write (isyspr,*)'Composition treated is ',i
      write (*,*)'ENTERING CALCPROP: BEGINS THE SULPHUR DOMAIN'
      write (*,*)'Composition treated is ',i
c     pause
 
      write (isysck,*)
      write (isysck,*) 'CALCPROP Composition: ',i
      write (isysck,*) 'CALCPROP XFe2O3:', xfe2o3
      write (isysck,*) 'CALCPROP Xfor: ', xfor
      write (isysck,*) 'CALCPROP Scat(5): ',scat(5)
      write (isysck,*) 'CALCPROP Redox: ' ,redox(i)
      write (isysck,*) 'CALCPROP Redoz: ',redoz(i)
      write (isysck,*)
 
      si=2.0*scat(1)
      p=2.5*scat(3)
      al=1.50*scat(4)
      fe=1.50*scat(5)
      wat=0.5*xwd(i)
      totale=si+p+al+fe+wat
      scat(1)=si/totale
      scat(2)=wat/totale
      scat(3)=p/totale
      scat(4)=al/totale
      scat(5)=fe/totale
      sio2f(i)=(scat(1)*(1.0-totcat(i)))
      po25f(i)=(scat(3)*(1.0-totcat(i)))
      alo15f(i)=(scat(4)*(1.0-totcat(i)))
      feo15f(i)=(scat(5)*(1.0-totcat(i)))
      ho05f(i)=(scat(2)*(1.0-totcat(i)))
C     totfor=(scat(1)+scat(3)+scat(4)+scat(5))*(1.0-totcat(i))
      ti=2.0*xcat(2)
      p=2.5*xcat(3)
      al=1.5*xcat(4)
      fe3=1.5*xcat(5)
      cr=1.5*xcat(6)
      fe2=xcat(7)*1.
      xmn=xcat(8)*1.
      ca=xcat(9)*1.
      xmg=xcat(10)*1.
      xna=xcat(11)*0.5
      xk=xcat(12)*0.5
      xhh=xcat(13)*0.5
      totale=ti+p+al+fe3+cr+fe2+xmn+ca+xmg+xna+xk+xhh
      xcat(2)=ti/totale
      xcat(3)=p/totale
      xcat(4)=al/totale
      xcat(5)=fe3/totale
      xcat(6)=cr/totale
      xcat(7)=fe2/totale
      xcat(8)=xmn/totale
      xcat(9)=ca/totale
      xcat(10)=xmg/totale
      xcat(11)=xna/totale
      xcat(12)=xk/totale
      xcat(13)=xhh/totale
      tio2m(i)=(xcat(2)*totcat(i))
      po25m(i)=(xcat(3)*totcat(i))
      alo15m(i)=(xcat(4)*totcat(i))
      feo15m(i)=(xcat(5)*totcat(i))
      cro15m(i)=(xcat(6)*totcat(i))
      feom(i)=xcat(7)*totcat(i)
      mnom(i)=xcat(8)*totcat(i)
      caom(i)=xcat(9)*totcat(i)
      mgom(i)=xcat(10)*totcat(i)
      nao05m(i)=(xcat(11)*totcat(i))
      ko05m(i)=(xcat(12)*totcat(i))
      ho05m(i)=(xcat(13)*totcat(i))
C     totmod=(xcat(2)+xcat(3)+xcat(4)+xcat(5)+xcat(6)+xcat(7)+
C    $ xcat(8)+xcat(9)+xcat(10)+xcat(11)+xcat(12)+xcat(13))*totcat(i)
 
      write (isyspr,*) 'FeO15f =', feo15f(i)
      write (isyspr,*) 'FeO15m =', feo15m(i)
      write (isyspr,*) 'FeOm =', feom(i)
      if (i.eq.ncomp) then
       write (isyspr,*) '-------------------------FINAL NEWSULF-CALLED S
     &TEP -----------'
       write (isyspr,*)
      endif
 
      END
 
      SUBROUTINE DATAOUT
 
* prints informations on melt composition and structure
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
c
      write(13,*)'--------------------'
      write(13,*)'composition n. ',i
      write(13,*)'N acidic=',acidic(i)
      write(13,*)'polymerization constant =',polcos(i)
      write(13,*)'polianions =',struct(i)
      write(13,*)'total anions (including S=  ) =',totani(i)
      write(13,*)'total anions (excluding S=  ) =',totold(i)
      write(13,*)'total cations =',totcat(i)
      write(13,*)'parameter p of the polymerization equation =',o4(i)
      write(13,*)'mean number of cations per polianion =',s
      write(13,*)'moles of non-bridging oxygen (NBO) (O-) =',o(i)
      write(13,*)'moles of free oxygen ion (O=) =',o1(i)
      write(13,*)'moles of bridging oxygen (Oø) =',o3(i)
      write(13,*)'NBO/T =',o(i)/ACIDIC(I)
      write(13,*)'total sulphur anions = ',stot(i)                                      '
      write(13,*)'activity of free oxygen ion aO= =',aossi(i)
      write(13,*)'activity of polyanions (struct) =',struct(i)/totani(i)
      write(13,*)'activity of indifferentiated sulphur  =',sobs(i)
      write(13,*)'trivalent ANION iron activity aFeO2-  =',afeo2(i)
      write(13,*)'New trivalent ANION iron activity aFeO2-  =',afeo22(i)
      write(13,*)'trivalent CATION iron activity aFe3+ =',afe3(i)
      write(13,*)'New trivalent CATION iron activity aFe3+ =',afe32(i)
      write(13,*)'CATION ferrous iron activity aFe2+  =',afe2(i)
      write(13,*)'CATION titanium activity aTi4+ =',ati(i)
      write(13,*)'CATION cromium activity aCr3+  =',acr(i)
      write(13,*)'CATION phosphorous activity aP5+  =',ap(i)
      write(13,*)'CATION manganese activity aMn2+  =',amn(i)
      write(13,*)'CATION magnesium activity aMg2+  =',amg(i)
      write(13,*)'CATION calcium activity aCa2+  =',aca(i)
      write(13,*)'CATION sodium activity aNa+  =',ana(i)
      write(13,*)'CATION potassium activity aK+  =',ak(i)
      write(13,*)'CATION hydrogen activity aH+  =',aprot(i)
      totti=afe2(i)+afe3(i)+ati(i)+acr(i)+ap(i)+amn(i)+amg(i)+aca(i)+
     $ana(i)+ak(i)+aprot(i)
      write(13,*)'SUM of CATION activities =',totti
      write(13,*)'                                      '
 
C     pause
      END
 
      SUBROUTINE SULPHIDE
****************************************************************
*  computes sulphide capacities and writes the results on      *
*  gasjob.out                                                  *
****************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
c     rapmed=rapmed/ncomp
 
      write (ISYSMF,*) '    x10        x23        misfit FCN'
 
      write (ISYSWR,*) '   TK    Pbars  logfO2  logfS2  logCs(o)  logCs
     &(c)  logKS2'
 
      write (16,*) '     aO=          Totani       totcat         dd
     &    kpol    Basdif  Stotobs      nO=       nOø     nO-   struct    
     $  xOH-   Nacidic  si     NBO/T      S2(wt%)    s6(wt%)    Stot(wt%)
     $   S6+/Stot  log(KSo4/KS2)   logfO2in    logfO2ric   XH2O   nHtot
     $ nH+    nOH-   nOH     NHres  som   IP   DVH+O-toOH   DVH+O=toOH-
     $ aSO42-   aS2-'
c     do 1235 ik= -10,3
      do 1235 ik = 3,3
 
        x10 = ik
        if (ik.eq.3) x10 = x(10)
c       x10 = x(10)
 
c      do 1236 il = 1,21
       do 1236 il = 21,21
 
        x23 = -25000.+il*1000.
        if (il.eq.21) x23 = x(23)
c       x23 = x(23)
 
      xxfit = 0.
 
      do 1234 i=1,ncomp
 
* * * * * * * * * DeltaG di volume per coppie ox-sulf * * * * * * * * * *
c     NEW VALUES
      x28=57.71
      x29=12.88
      x30=-7.39
c     X31 è una differenza tra solfato e sulfuro, not dep on P e T
c      x31= 40.   !!! Guess Value, change it if you like!!!
c      x31= 45.1937964902720 ! Papale model adopted 24 March 2004
c      x31= 40.6234 ! Burnham model adopted 24 March 2004
      x31= 41 ! Papale model adopted 13 August 2004

c     OLD VALUES
c     x28=34.405
c     x29=0.374
c     x30=-2.44
c     x31= 20.   !!! Guess Value, change it if you like!!!
 
c                       SiO2-SiS2
      vsi1o =voxm(1)+eoxm(1)*0.001*(tkelv(i)-1673.)
      vsi1s=(voxm(1)+2.*x28)+(eoxm(1)+2.*x29)*0.001*(tkelv(i)-1673)
     &     -2.*x31
c   terms above are converted in joule/bar
      vsi1o =0.1*vsi1o
      vsi1s =0.1*vsi1s
 
      csiox =coxm(1)*0.0001*0.1
      csisu =(coxm(1)+2.*x30)*0.0001*0.1
 
      vsi=(0.5*vsi1s-0.5*vsi1o)*(pbar(i)-1.)+(0.5*csisu-0.5*csiox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
c
c   term becomes term/RT
      vsi =vsi/(tkelv(i)*8.3147)
 
c                       TiO2-TiS2
      vti1o =voxm(2)+eoxm(2)*0.001*(tkelv(i)-1673.)
      vti1s=(voxm(2)+2.*x28)+(eoxm(2)+2.*x29)*0.001*(tkelv(i)-1673)
     &     -2.*x31
c   terms above are converted in joule/bar
      vti1o =0.1*vti1o
      vti1s =0.1*vti1s
 
      ctiox =coxm(2)*0.0001*0.1
      ctisu =(coxm(2)+2.*x30)*0.0001*0.1
 
      vti = (0.5*vti1s-0.5*vti1o)*(pbar(i)-1.)+(0.5*ctisu-0.5*ctiox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vti =vti/(tkelv(i)*8.3147)
 
c                       P2O5-P2S5
      vp1o =voxm(3)+eoxm(3)*0.001*(tkelv(i)-1673.)
      vp1s=(voxm(3)+5.*x28)+(eoxm(3)+5.*x29)*0.001*(tkelv(i)-1673)
     &     -2.*x31
c   terms above are converted in joule/bar
      vp1o =0.1*vp1o
      vp1s =0.1*vp1s
 
      cpox =coxm(3)*0.0001*0.1
      cpsu =(coxm(3)+5.*x30)*0.0001*0.1
 
      vp=(0.2*vp1s-0.2*vp1o)*(pbar(i)-1.)+(0.2*cpsu-0.2*cpox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vp =vp/(tkelv(i)*8.3147)
 
c                       Al2O3-Al2S3
      val1o =voxm(4)+eoxm(4)*0.001*(tkelv(i)-1673.)
      val1s=(voxm(4)+3.*x28)+(eoxm(4)+3.*x29)*0.001*(tkelv(i)-1673)
     &     -3.*x31
c   terms above are converted in joule/bar
      val1o =0.1*val1o
      val1s =0.1*val1s
 
      calox =coxm(4)*0.0001*0.1
      calsu =(coxm(4)+3.*x30)*0.0001*0.1
 
      val=(1./3.*val1s-1./3.*val1o)*(pbar(i)-1.)+(1./3.*calsu-1./3.
     $*calox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      val =val/(tkelv(i)*8.3147)
 
c                       Fe2O3-Fe2S3
      vfe31o =voxm(5)+eoxm(5)*0.001*(tkelv(i)-1673.)
      vfe31s=(voxm(5)+3.*x28)+(eoxm(5)+3.*x29)*0.001*(tkelv(i)-1673)
     &     -3.*x31
c   terms above are converted in joule/bar
      vfe31o =0.1*vfe31o
      vfe31s =0.1*vfe31s
 
      cfe3ox =coxm(5)*0.0001*0.1
      cfe3su =(coxm(5)+3.*x30)*0.0001*0.1
 
      vfe3=(1./3.*vfe31s-1./3.*vfe31o)*(pbar(i)-1.)+(1./3.*cfe3su-1./3.
     $*cfe3ox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vfe3 =vfe3/(tkelv(i)*8.3147)
 
c                       Cr2O3-Cr2S3
      vcr1o =voxm(6)+eoxm(6)*0.001*(tkelv(i)-1673.)
      vcr1s=(voxm(6)+3.*x28)+(eoxm(6)+3.*x29)*0.001*(tkelv(i)-1673)
     &     -3.*x31
c   terms above are converted in joule/bar
      vcr1o =0.1*vcr1o
      vcr1s =0.1*vcr1s
 
      ccrox =coxm(6)*0.0001*0.1
      ccrsu =(coxm(6)+3.*x30)*0.0001*0.1
 
      vcr=(1./3.*vcr1s-1./3.*vcr1o)*(pbar(i)-1.)+(1./3.*ccrsu-1./3.*
     $ccrox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vcr =vcr/(tkelv(i)*8.3147)
 
c                       FeO-FeS
      vfe21o =voxm(7)+eoxm(7)*0.001*(tkelv(i)-1673.)
      vfe21s=(voxm(7)+x28)+(eoxm(7)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vfe21o =0.1*vfe21o
      vfe21s =0.1*vfe21s
 
      cfe2ox =coxm(7)*0.0001*0.1
      cfe2su =(coxm(7)+x30)*0.0001*0.1
 
      vfe2=(vfe21s-vfe21o)*(pbar(i)-1.)+(cfe2su-
     $cfe2ox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vfe2 =vfe2/(tkelv(i)*8.3147)
 
c                       MnO-MnS
      vmn1o =voxm(8)+eoxm(8)*0.001*(tkelv(i)-1673.)
      vmn1s=(voxm(8)+x28)+(eoxm(8)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vmn1o =0.1*vmn1o
      vmn1s =0.1*vmn1s
 
      cmnox =coxm(8)*0.0001*0.1
      cmnsu=(coxm(8)+x30)*0.0001*0.1
 
      vmn=(vmn1s-vmn1o)*(pbar(i)-1.)+(cmnsu-
     $cmnox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vmn =vmn/(tkelv(i)*8.3147)
 
c                       CaO-CaS
      vca1o =voxm(9)+eoxm(9)*0.001*(tkelv(i)-1673.)
      vca1s=(voxm(9)+x28)+(eoxm(9)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vca1o =0.1*vca1o
      vca1s =0.1*vca1s
 
      ccaox =coxm(9)*0.0001*0.1
      ccasu=(coxm(9)+x30)*0.0001*0.1
 
      vca=(vca1s-vca1o)*(pbar(i)-1.)+(ccasu-
     $ccaox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vca =vca/(tkelv(i)*8.3147)
 
c                       MgO-MgS
      vmg1o =voxm(10)+eoxm(10)*0.001*(tkelv(i)-1673.)
      vmg1s=(voxm(10)+x28)+(eoxm(10)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vmg1o =0.1*vmg1o
      vmg1s =0.1*vmg1s
 
      cmgox =coxm(10)*0.0001*0.1
      cmgsu=(coxm(10)+x30)*0.0001*0.1
 
      vmg=(vmg1s-vmg1o)*(pbar(i)-1.)+(cmgsu-
     $cmgox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vmg =vmg/(tkelv(i)*8.3147)
 
c                       Na2O-Na2S
      vna1o =voxm(11)+eoxm(11)*0.001*(tkelv(i)-1673.)
      vna1s=(voxm(11)+x28)+(eoxm(11)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vna1o =0.1*vna1o
      vna1s =0.1*vna1s
 
      cnaox =coxm(11)*0.0001*0.1
      cnasu=(coxm(11)+x30)*0.0001*0.1
 
      vna=(vna1s-vna1o)*(pbar(i)-1.)+(cnasu-
     $cnaox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vna =vna/(tkelv(i)*8.3147)
 
c                       K2O-K2S
      vk1o =voxm(12)+eoxm(12)*0.001*(tkelv(i)-1673.)
      vk1s=(voxm(12)+x28)+(eoxm(12)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vk1o =0.1*vk1o
      vk1s =0.1*vk1s
 
      ckox =coxm(12)*0.0001*0.1
      cksu=(coxm(12)+x30)*0.0001*0.1
 
      vk=(vk1s-vk1o)*(pbar(i)-1.)+(cksu-
     $ckox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vk =vk/(tkelv(i)*8.3147)
 
c                       H2O-H2S
      vh1o =voxm(13)+eoxm(13)*0.001*(tkelv(i)-1673.)
      vh1s=(voxm(13)+x28)+(eoxm(13)+x29)*0.001*(tkelv(i)-1673)
     &     -x31
c   terms above are converted in joule/bar
      vh1o =0.1*vh1o
      vh1s =0.1*vh1s
 
      chox =coxm(13)*0.0001*0.1
      chsu=(coxm(13)+x30)*0.0001*0.1
 
      vh=(vh1s-vh1o)*(pbar(i)-1.)+(chsu-
     $chox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vh=vh/(tkelv(i)*8.3147)
c     vh=0.
 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
      x227=0.
c  Below annealing value to adopt with Papale  water solubility model
c     x227=0.
c  Below annealing value to adopt with Burnham water solubility model
c     x227=0.
 
      reg(i)=sio2f(i)*(x(1)+x(14)/tkelv(i)-0.5*x227*(pbar(i)-1)-vsi)
     $   +tio2m(i)*(x(2)+x(15)/tkelv(i)-0.5*x227*(pbar(i)-1)-vti)
     $   +po25f(i)*(x(3)+x(16)/tkelv(i)-0.2*x227*(pbar(i)-1)-vp)
     $   +po25m(i)*(x(3)+x(16)/tkelv(i)-0.2*x227*(pbar(i)-1)-vp)
     $   +alo15f(i)*(x(4)+x(17)/tkelv(i)-1./3.*x227*(pbar(i)-1)-val)
     $   +alo15m(i)*(x(4)+x(17)/tkelv(i)-1./3.*x227*(pbar(i)-1)-val)
     $   +feo15f(i)*(x(5)+x(18)/tkelv(i)-1./3.*x227*(pbar(i)-1)-vfe3)
     $   +feo15m(i)*(x(5)+x(18)/tkelv(i)-1./3.*x227*(pbar(i)-1)-vfe3)
     $   +cro15m(i)*(x(6)+x(19)/tkelv(i)-1./3.*x227*(pbar(i)-1)-vcr)
     $   +feom(i)*(x(7)+x(20)/tkelv(i)-x227*(pbar(i)-1)-vfe2)
     $   +mnom(i)*(x(8)+x(21)/tkelv(i)-x227*(pbar(i)-1)-vmn)
     $   +caom(i)*(x(9)+x(22)/tkelv(i)-x227*(pbar(i)-1)-vca)
     $   +mgom(i)*(x(10)+x(23)/tkelv(i)-x227*(pbar(i)-1)-vmg)
     $   +nao05m(i)*(x(11)+x(24)/tkelv(i)-x227*(pbar(i)-1)-vna)
     $   +ko05m(i)*(x(12)+x(25)/tkelv(i)-x227*(pbar(i)-1)-vk)
     $   +ho05m(i)*(x(13)+x(26)/tkelv(i)-x227*(pbar(i)-1)-vh)
     $   +ho05f(i)*(x(13)+x(26)/tkelv(i)-x227*(pbar(i)-1)-vh)

c      write (*,*) 'mgom x10  x23 reg',mgom(i),x101,x231,reg(i)
c      write (*,*) '------',x(10),x(23)
c      pause
      write (15,*) '+ + + BEGINNING S CAPACITIES BLOCK + + +'
      write (15,*) 'COMP ',i,' feo15f =', feo15f(i)
      write (15,*) 'COMP ',i,' feo15m=', feo15m(i)
      write (15,*) 'COMP ',i,' feom=', feom(i)
      write (15,*) ' + + + END S CAPACITIES LOG BLOCK + + +'
c     pause
      REG22=0.
      REG2(I)=REG22
      REG332=0.
      REG3(I)=REG332
c     reg3(i)=sio2f(i)*vsi+tio2m(i)*vti+(po25m(i)+po25f(i))*vp
c    $        +(alo15f(i)+alo15m(i))*val+(feo15f(i)+feo15m(i))*vfe3
c    $        +cro15m(i)*vcr+feom(i)*vfe2+mnom(i)*vmn+caom(i)*vca
c    $        +mgom(i)*vmg+nao05m(i)*vna+ko05m(i)*vk+ho05m(i)*vh
 
      reg(i)=dexp(reg(i))
 
c     Cscal=reg(i)*sppm(i)*totani(i)*dd(i)*32.064*beppe(i)**0.5
      Cscal=reg(i)*aozzi(i)*totali(i)*dd(i)*32.064
       write (*,*) reg(i),aozzi(i),totali(i),dd(i)
c       pause

c      write (*,*) 'reg aozzi totali dd foga pbar'
c      write (*,*) reg(i),aozzi(i),totali(i),dd(i),foga(i),pbar(i)
      write (isysmf,*) 'SULFURI @ COMP ',i
      write (isysmf,*) 'reg = ',reg(i)
      write (isysmf,*) 'aozzi = ',aozzi(i)
      write (isysmf,*) 'totali = ',totali(i)
      write (isysmf,*) 'dd = ',dd(i)

c     Cscal=cscal/(dexp(reg3(i)))
      S2MENO(I)=CSCAL/(BEPPE(I)**0.5)
      s2k(i)=cscal/(aozzi(i)*totali(i)*dd(i)*32.064)
      s2k(i)=dlog10(s2k(i))
      Csobs=sobs(i)*totali(i)*dd(i)*32.064*beppe(i)**0.5
      if (csobs.eq.0.) csobs=1.
      csc=dlog10(cscal)
      cso=dlog10(csobs)
      cs(i)=csc
 
c   Calculation of MISFIT function for sulphide solubility
 
      if (cso.eq.0.) cso=1.
      xfit=(csc-cso)**2./dabs(cso)
c     write (19,*) x10,x23,xfit
      WRITE(ISYSWR,123) tkelv(I),pbar(I),xossi(I),fs2(I),
     $dlog10(Csobs),dlog10(Cscal),s2k(i)
123   format (f8.2,f7.0,5f9.3)
      xxfit = xxfit + xfit
 
1234  CONTINUE
      write (isysmf,788)  x10, x23, xxfit
788   format (x,f8.3,x,f12.3,x,d12.5)
1236  continue
1235  continue
 
16    CONTINUE
      END
 
      SUBROUTINE SULPHATE
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  computes sulphate capacities and writes the results on     *
*  gasjob.out                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=300)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   xx(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13),
     $         ffeo2m(ndi),ffe2(ndi),ffe3m(ndi),ffo2m(ndi),fa3fa2(ndi),
     $         fa3f2a(ndi),a3a2z(ndi),ca3c2a(ndi),fcox(ndi),fcoz(ndi)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),wmol(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,  ISYSWZ, ISYSMF,isysvo,
     $         isysho,isysfe
     6/COUNT/  I,kj,jk,iter
     7/FLAG/   kir(ndi),jwflag
     8/WASSER/ partwa,partwb
 
      dimension s6stot(ndi)
C     rapmed=rapmed/ncomp
 
      write (ISYSWZ,*) 'TK    Pbars   LogCs(o)   LogCs(c) LogKSO4
     $ logfO2in    logfO2cal'
      write (isysfe,*)'TK    Pbars  xH2O  afeo2  afeo22  l[FeO2]/[Fe2]
     & l[Fe3]/[Fe] lfFeO2/fFe2 lfFe3/Fe2 fFe fO2-  afe2   afe22   afe3
     & afe32   redox   redoz refirst   act3/2X act3/2Z fug3/2X fug3/2Z
     $ gam3/2X gam3/2Z  fcof3/2X fcof3/3Z logfO2in   logfO2cal  KFLAG'
c
      write (isysho,*) 'REMINDER: '
      write (isysho,*) '1. (OH-/H)<(nH2Otot+1/2nOHIR)/(nH2Otot-1/2nOHIR)
     $    cioŠ <dx  =>nOHveri<nOHIR'
      write (isysho,*) '2. (OH-/H)>(nH2Otot-1/2nOHIR)/(nH2Otot+1/2nOHIR)
     $    cioŠ >left =>nOHveri>0'
      write (isysho,*)
      write (isysho,*) 'H2Ototwt% H2Omolwt% nH2Otot nH2Omol nOH-IR nHtot
     $ nOHveri nHveri Totcat  nO2-  nO-  nOø TOHexpected O-residual
     $    1/T  TK  >left   <dx    OH-/H+   nH+   nOH-   cost22  >LEFT
     $  <DX  TøC  Pbar  xOH-'

      write (isysvo,*)
 
      DO 1234 I=1,ncomp
 
c     afeo22(i)=afeo22(i)/afe22(i)
c     afe32(i)=afe32(i)/afe22(i)
c     afeo2m(i)=afeo2m(i)/afe22(i)
c     afe3m(i)=afe3m(i)/afe22(i)
 
****************** DeltaG di volume per coppie ox-sulf ****************
c     NEW VALUES
      x28=57.71
      x29=12.88
      x30=-7.39
c     OLD VALUES
c     x28=34.405
c     x29=0.374
c     x30=-2.44
 
c                       SiO2-Si(SO4)2
      vsi1o =voxm(1)+eoxm(1)*0.001*(tkelv(i)-1673.)
      vsi1s=(voxm(1)+2.*x28)+(eoxm(1)+2.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vsi1o =0.1*vsi1o
      vsi1s =0.1*vsi1s
 
      csiox =coxm(1)*0.0001*0.1
      csisu =(coxm(1)+2.*x30)*0.0001*0.1
 
      vsi=(0.5*vsi1s-0.5*vsi1o)*(pbar(i)-1.)+(0.5*csisu-0.5*csiox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
c
c   term becomes term/RT
      vsi =vsi/(tkelv(i)*8.3147)
 
c                       TiO2-Ti(SO4)2
      vti1o =voxm(2)+eoxm(2)*0.001*(tkelv(i)-1673.)
      vti1s=(voxm(2)+2.*x28)+(eoxm(2)+2.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vti1o =0.1*vti1o
      vti1s =0.1*vti1s
 
      ctiox =coxm(2)*0.0001*0.1
      ctisu =(coxm(2)+2.*x30)*0.0001*0.1
 
      vti = (0.5*vti1s-0.5*vti1o)*(pbar(i)-1.)+(0.5*ctisu-0.5*ctiox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vti =vti/(tkelv(i)*8.3147)
 
c                       P2O5-P2(SO4)5
      vp1o =voxm(3)+eoxm(3)*0.001*(tkelv(i)-1673.)
      vp1s=(voxm(3)+5.*x28)+(eoxm(3)+5.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vp1o =0.1*vp1o
      vp1s =0.1*vp1s
 
      cpox =coxm(3)*0.0001*0.1
      cpsu =(coxm(3)+5.*x30)*0.0001*0.1
 
      vp=(0.2*vp1s-0.2*vp1o)*(pbar(i)-1.)+(0.2*cpsu-0.2*cpox)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vp =vp/(tkelv(i)*8.3147)
 
c                       Al2O3-Al2(SO4)3
      val1o =voxm(4)+eoxm(4)*0.001*(tkelv(i)-1673.)
      val1s=(voxm(4)+3.*x28)+(eoxm(4)+3.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      val1o =0.1*val1o
      val1s =0.1*val1s
 
      calox =coxm(4)*0.0001*0.1
      calsu =(coxm(4)+3.*x30)*0.0001*0.1
 
      val=(1./3.*val1s-1./3.*val1o)*(pbar(i)-1.)+(1./3.*calsu-1./3.
     $*calox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      val =val/(tkelv(i)*8.3147)
 
c                       Fe2O3-Fe2(SO4)3
      vfe31o =voxm(5)+eoxm(5)*0.001*(tkelv(i)-1673.)
      vfe31s=(voxm(5)+3.*x28)+(eoxm(5)+3.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vfe31o =0.1*vfe31o
      vfe31s =0.1*vfe31s
 
      cfe3ox =coxm(5)*0.0001*0.1
      cfe3su =(coxm(5)+3.*x30)*0.0001*0.1
 
      vfe3=(1./3.*vfe31s-1./3.*vfe31o)*(pbar(i)-1.)+(1./3.*cfe3su-1./3.
     $*cfe3ox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vfe3 =vfe3/(tkelv(i)*8.3147)
 
c                       Cr2O3-Cr2(SO4)3
      vcr1o =voxm(6)+eoxm(6)*0.001*(tkelv(i)-1673.)
      vcr1s=(voxm(6)+3.*x28)+(eoxm(6)+3.*x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vcr1o =0.1*vcr1o
      vcr1s =0.1*vcr1s
 
      ccrox =coxm(6)*0.0001*0.1
      ccrsu =(coxm(6)+3.*x30)*0.0001*0.1
 
      vcr=(1./3.*vcr1s-1./3.*vcr1o)*(pbar(i)-1.)+(1./3.*ccrsu-1./3.*
     $ccrox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vcr =vcr/(tkelv(i)*8.3147)
 
c                       FeO-FeSO4
      vfe21o =voxm(7)+eoxm(7)*0.001*(tkelv(i)-1673.)
      vfe21s=(voxm(7)+x28)+(eoxm(7)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vfe21o =0.1*vfe21o
      vfe21s =0.1*vfe21s
 
      cfe2ox =coxm(7)*0.0001*0.1
      cfe2su =(coxm(7)+x30)*0.0001*0.1
 
      vfe2=(vfe21s-vfe21o)*(pbar(i)-1.)+(cfe2su-
     $cfe2ox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vfe2 =vfe2/(tkelv(i)*8.3147)
 
c                       MnO-MnSO4
      vmn1o =voxm(8)+eoxm(8)*0.001*(tkelv(i)-1673.)
      vmn1s=(voxm(8)+x28)+(eoxm(8)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vmn1o =0.1*vmn1o
      vmn1s =0.1*vmn1s
 
      cmnox =coxm(8)*0.0001*0.1
      cmnsu=(coxm(8)+x30)*0.0001*0.1
 
      vmn=(vmn1s-vmn1o)*(pbar(i)-1.)+(cmnsu-
     $cmnox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vmn =vmn/(tkelv(i)*8.3147)
 
c                       CaO-CaSO4
      vca1o =voxm(9)+eoxm(9)*0.001*(tkelv(i)-1673.)
      vca1s=(voxm(9)+x28)+(eoxm(9)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vca1o =0.1*vca1o
      vca1s =0.1*vca1s
 
      ccaox =coxm(9)*0.0001*0.1
      ccasu=(coxm(9)+x30)*0.0001*0.1
 
      vca=(vca1s-vca1o)*(pbar(i)-1.)+(ccasu-
     $ccaox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vca =vca/(tkelv(i)*8.3147)
 
c                       MgO-MgSO4
      vmg1o =voxm(10)+eoxm(10)*0.001*(tkelv(i)-1673.)
      vmg1s=(voxm(10)+x28)+(eoxm(10)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vmg1o =0.1*vmg1o
      vmg1s =0.1*vmg1s
 
      cmgox =coxm(10)*0.0001*0.1
      cmgsu=(coxm(10)+x30)*0.0001*0.1
 
      vmg=(vmg1s-vmg1o)*(pbar(i)-1.)+(cmgsu-
     $cmgox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vmg =vmg/(tkelv(i)*8.3147)
 
c                       Na2O-Na2SO4
      vna1o =voxm(11)+eoxm(11)*0.001*(tkelv(i)-1673.)
      vna1s=(voxm(11)+x28)+(eoxm(11)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vna1o =0.1*vna1o
      vna1s =0.1*vna1s
 
      cnaox =coxm(11)*0.0001*0.1
      cnasu=(coxm(11)+x30)*0.0001*0.1
 
      vna=(vna1s-vna1o)*(pbar(i)-1.)+(cnasu-
     $cnaox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vna =vna/(tkelv(i)*8.3147)
 
c                       K2O-K2SO4
      vk1o =voxm(12)+eoxm(12)*0.001*(tkelv(i)-1673.)
      vk1s=(voxm(12)+x28)+(eoxm(12)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vk1o =0.1*vk1o
      vk1s =0.1*vk1s
 
      ckox =coxm(12)*0.0001*0.1
      cksu=(coxm(12)+x30)*0.0001*0.1
 
      vk=(vk1s-vk1o)*(pbar(i)-1.)+(cksu-
     $ckox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vk =vk/(tkelv(i)*8.3147)
 
c                       H2O-H2SO4
      vh1o =voxm(13)+eoxm(13)*0.001*(tkelv(i)-1673.)
      vh1s=(voxm(13)+x28)+(eoxm(13)+x29)*0.001*(tkelv(i)-1673)
c   terms above are converted in joule/bar
      vh1o =0.1*vh1o
      vh1s =0.1*vh1s
 
      chox =coxm(13)*0.0001*0.1
      chsu=(coxm(13)+x30)*0.0001*0.1
 
      vh=(vh1s-vh1o)*(pbar(i)-1.)+(chsu-
     $chox)*(pbar(i)**2./2.-pbar(i)+0.5)
 
c   term becomes term/RT
      vh=vh/(tkelv(i)*8.3147)
c     vh=0.
 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
      x27=0.
c     Below SO4 annealing value to adopt with Papale  water solubility model
      x27=-0.00407081493980624   !LAST: 7 April 2005
c     Below SO4 annealing value to adopt with Burnham water solubility model
c      x27=-0.00344957395889646   !LAST: 13 August  2004
 
      regg(i)=sio2f(i)*(xx(1)+xx(14)/tkelv(i)-0.5*x27*(pbar(i)-1)-vsi)
     $ +tio2m(i)*(xx(2)+xx(15)/tkelv(i)-0.5*x27*(pbar(i)-1)-vti)
     $ +po25f(i)*(xx(3)+xx(16)/tkelv(i)-0.2*x27*(pbar(i)-1)-vp)
     $ +po25m(i)*(xx(3)+xx(16)/tkelv(i)-0.2*x27*(pbar(i)-1)-vp)
     $ +alo15f(i)*(xx(4)+xx(17)/tkelv(i)-1./3.*x27*(pbar(i)-1)-val)
     $ +alo15m(i)*(xx(4)+xx(17)/tkelv(i)-1./3.*x27*(pbar(i)-1)-val)
     $ +feo15f(i)*(xx(5)+xx(18)/tkelv(i)-1./3.*x27*(pbar(i)-1)-vfe3)
     $ +feo15m(i)*(xx(5)+xx(18)/tkelv(i)-1./3.*x27*(pbar(i)-1)-vfe3)
     $ +cro15m(i)*(xx(6)+xx(19)/tkelv(i)-1./3.*x27*(pbar(i)-1)-vcr)
     $ +feom(i)*(xx(7)+xx(20)/tkelv(i)-x27*(pbar(i)-1)-vfe2)
     $ +mnom(i)*(xx(8)+xx(21)/tkelv(i)-x27*(pbar(i)-1)-vmn)
     $ +caom(i)*(xx(9)+xx(22)/tkelv(i)-x27*(pbar(i)-1)-vca)
     $ +mgom(i)*(xx(10)+xx(23)/tkelv(i)-x27*(pbar(i)-1)-vmg)
     $ +nao05m(i)*(xx(11)+xx(24)/tkelv(i)-x27*(pbar(i)-1)-vna)
     $ +ko05m(i)*(xx(12)+xx(25)/tkelv(i)-x27*(pbar(i)-1)-vk)
     $ +ho05m(i)*(xx(13)+xx(26)/tkelv(i)-x27*(pbar(i)-1)-vh)
     $ +ho05f(i)*(xx(13)+xx(26)/tkelv(i)-x27*(pbar(i)-1)-vh)
 
      REG32=0.
      REGG3(I)=REG32
      reg43=0.
      reg4(i)=reg43
 
c      reg4(i)=sio2f(i)*vsi+tio2m(i)*vti+(po25m(i)+po25f(i))*vp
c    $        +(alo15f(i)+alo15m(i))*val+(feo15f(i)+feo15m(i))*vfe3
c     $        +cro15m(i)*vcr+feom(i)*vfe2+mnom(i)*vmn+caom(i)*vca
c     $        +mgom(i)*vmg+nao05m(i)*vna+ko05m(i)*vk+ho05m(i)*vh
 
      regg(i)=dexp(regg(i))
      Cscal2=regg(i)*sppb(i)*totali(i)*dd(i)*32.064/beppi(i)
c      write (*,*) 'regg sppb totali dd beppi foga pbar'
c      write (*,*) regg(i),sppb(i),totali(i),dd(i),beppi(i),foga(i),
c     $ pbar(i)
c      write (21,145) dlog(regg(i)),dlog(sppb(i)),totali(i),dd(i),
c     $ dlog(beppi(I)),foga(i), pbar(i)
145   format(2(f8.4,x),f9.6,x,f9.6,x,f8.4,x,f9.3,x,f9.3)
      write (isysmf,*) 'SOLFATI @ COMP ',i
      write (isysmf,*) 'regg = ',regg(i)
      write (isysmf,*) 'sppb = ',sppb(i)
      write (isysmf,*) 'totali = ',totali(i)
      write (isysmf,*) 'dd = ',dd(i)
      write (isysmf,*) 'beppi = ',beppi(i)

c     Cscal2=cscal2/(dexp(reg4(i)))
C     Cscal2=regg(i)*aozzi(i)*totali(i)*dd(i)*32.064
      Csobs2=sobs(i)*totali(i)*dd(i)*32.064/beppi(i)
      so4k(i)=cscal2/(sppb(i)*totali(i)*dd(i)*32.064/beppi(i))
      so4k(i)=dlog10(so4k(i))
      wec=so4k(i)-s2k(i)
      if (csobs2.eq.0.) csobs2=1.
      SO4(I)=cscal2*BEPPI(I)
      TOTS(I)=SO4(I)+S2MENO(I)
      s6stot(i)=so4(i)/tots(i)
      OXSOLF(I)=(SO4(I)/32.064)*4*15.9994
      csc2=dlog10(cscal2)
      cso2=dlog10(csobs2)
      write (*,*)
c      pause

      write (*,*)
      write(*,*)' COMP    S=      SO4=       STOT  S6/Sot    TK     P'
      write(*,901)i,S2MENO(I),SO4(I),TOTS(I),S6STOT(i),Tkelv(i),pbar(i)
      write (*,*)
      write(*,*) ' logfo2     logfO2rec     logfS2 '
      write(*,902) xossi(i),bossi(i),fs2(i)
902   format (3(x,f9.3))
      write (*,*)
      write (ivideo,*) '    logCs  logCso4  logKS2-  logKSO4   wec'
      write (ivideo,*) cs(i),csc2,s2k(i),so4k(i),wec
      write (ivideo,*) 'S6+/Stot'
      write (ivideo,*) s6stot(i)
      write (ivideo,*) 'Scalc = ',tots(i)*10000.
      write (ivideo,*) 'sobs =',sulf(i)
c      if (i.eq.ncomp) pause
       pause
      WRITE(ISYSWZ,123) tkelv(I),pbar(I),dlog10(Csobs2),
     $dlog10(Cscal2),so4k(i),xossi(i),bossi(i)
      write (isysfe,497) tkelv(i),pbar(i),xh2o(i),afeo2(i),afeo22(i)
     $,afeo2m(i),afe3m(i),ffeo2m(i),ffe3m(i),ffe2(i),ffo2m(i),afe2(i),
     $afe22(i),afe3(i),afe32(i),redox(i),redoz(i),refirst(i),a3a2(i),
     $a3a2z(i),fa3fa2(i),fa3f2a(i),ca3ca2(i),ca3c2a(i),fcox(i),fcoz(i),
     $xossi(i),bossi(i),kflag(i)
      write (ISYSPU,965) aOSSI(I),Totani(I),Totcat(I),DD(I),POLCOS(I),
     &basdif(i),Stott(i)*32.064,o1(I),o3(I),o(i),struct(i),xwd(i),
     $ACIDIC(I),sii(i),O(I)/ACIDIC(I),
     &s2meno(i),so4(i),tots(i),s6stot(i),wec,xossi(i),bossi(i),xh2o(i)
     &,xcat0(i),xh(i),xoh(i),root(i)*som(i),xh(i)-root(i)*som(i),som(i),
     &por(i),dv(i),delv(i),sobs(i)*s6stot(i),sobs(i)*(1-s6stot(i))


      if (poss(i,13).eq.0.or.kir(i).eq.1) goto 1234
      sann=0.218216405605124d-3
      cost22=1/10.**(partwa+partwb/tkelv(i)+0.5*sann*(pbar(i)-1.))
      freeox=aozzi(i)*totali(i)
      wmol(i)=(xcat0(i)/2.-wmol(i))*2. ! moli di OH IR-like
c     wattot=poss(i,13)/y(13)
      wattot=xcat0(i)*0.5
      watmol=wattot-0.5*(wmol(i))
      aden=-totcat(i)*wattot+cost22*wattot*freeox+0.5*totcat(i)*wmol(i)+
     $0.5*cost22*freeox*wmol(i)
      bden=cost22*freeox+totcat(i)
      ohtrue=aden/bden
c      if (ohtrue.lt.0.) ohtrue=0.
      prot=wmol(i)-ohtrue
       xleft=(wattot-0.5*wmol(i))/(wattot+0.5*wmol(i))
       dx=(wattot+0.5*wmol(i))/(wattot-0.5*wmol(i))
      write (isysho,451) poss(I,13),watmol*18.015,wattot,watmol,wmol(i),
     $xcat0(i),ohtrue,prot,totcat(i),freeox,o(i),o3(i),prot,
     $o(i)-prot,1./tkelv(i),tkelv(i),xleft,dx,cost22*freeox/totcat(i),
     $xh(i),xoh(i),dlog10(cost22),dlog10(xleft*totcat(i)/freeox),
     $dlog10(dx*totcat(i)/freeox),tkelv(i)-273.15,pbar(i)
     $,xoh(i)/xcat0(i)

451   format (8(f8.5,x),f9.6,1x,f8.3,1x,14(f12.5,x),1x,f9.3,f11.4,x
     *,f8.5)

c     pause
123   format (f8.2,f7.0,x,5(f10.7,x))
497   format (f8.2,f7.0,x,26(f11.7,x),i3)
965   format (13f12.6,x,f10.4,x,3(d11.6,x),f12.6,x,2(f12.6,x),
     &14(f10.6,x))
901   format (2x,i3,4f15.6,x,f7.2,x,f7.0)
1234  CONTINUE
      do 3333 i=1,ncomp
       write (ivideo,*)
       write (ivideo,*) 'afeo2-new  log[feo2/fe2]-K57   afe3+new   log[
     &fe3/fe2]-k58'
       write (ivideo,225) afeo22(i),afeo2m(i),afe32(i),afe3m(i)
3333  continue
c      nuovo OUT
      write (isysvo,*) 'TøC  T(K)  vFeO1.5   vFeO     vFeO2-    vFe3+
     $    vFe2+    vo2m     vhmvoh    vH+      vOH-'
      do 334 jj=1,20
      temp=700.+jj*50.
      temp2=temp+273.15
      vfeO15=21.065+4.545*0.001*(Temp2-1673)
      vfeO=13.65+2.92*0.001*(Temp2-1673)
 
      conv=0.6022045
      vfeo2m=4./3.*3.14159*(3.29**3)
      vo2m=4./3.*3.14159*(1.40**3)
      vfe3p=4./3.*3.14159*(0.645**3)
      vfe2p=4./3.*3.14159*(0.78**3)
      voh=(4./3.)*3.14159*(1.40**3)
      vh=(4./3.)*3.14159*(0.**3)
 
      vfeo2m=vfeo2m*conv
      vo2m=vo2m*conv
      vfe3p=vfe3p*conv
      vfe2p=vfe2p*conv
      voh=voh*conv
      vh=vh*conv
      vfeo2m=vfeo2m-vo2m
      vhmvoh=vh+vo2m-voh
      write (isysvo,336) temp,temp2,vfeo15,vfeo,vfeo2m,vfe3p,
     $ vfe2p,vo2m,vhmvoh,vH,voh
336   format (x,f5.0,x,f7.2,x,9(f7.2,x))
334   continue
c fine nuovo printout
225   format (4(f10.7,x))
      END
 
