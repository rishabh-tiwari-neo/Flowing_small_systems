c.....driver program for A Multi-Phase Transport model
      PROGRAM AMPT
c
      double precision xmp, xmu, alpha, rscut2, cutof2, dshadow
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER FRAME*8, PROJ*8, TARG*8
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
clin-4/2012-6/2009:
c      common/phidcy/iphidcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
clin-7/2009:
      common/cmsflag/dshadow,ishadow
clin-2/2012 allow random orientation of reaction plane:
      common /phiHJ/iphirp,phiRP
cfast: output flags, see root_init below
      common /fastio/ itxtdt, izpcdt

      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************
      OPEN (24, FILE = 'input.ampt', STATUS = 'UNKNOWN')
      OPEN (12, FILE = 'ana/version', STATUS = 'UNKNOWN')
      READ (24, *) EFRM
c     format-read characters (for ALPHA compilers):
      READ (24, 111) FRAME
      READ (24, 111) PROJ
      READ (24, 111) TARG
      READ (24, *) IAP
      READ (24, *) IZP
      READ (24, *) IAT
      READ (24, *) IZT
      READ (24, *) NEVNT
      READ (24, *) BMIN
      READ (24, *) BMAX
c     flag to select default AMPT or string melting:
      READ (24, *) isoft
c     read initialization value for hadron cascade:
      READ (24, *) NTMAX
      READ (24, *) DT
c     parj(41) and (42) are a and b parameters in Lund string fragmentation:
      READ (24, *) PARJ(41)
      READ (24, *) PARJ(42)
c     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and 
c     increase the net-baryon stopping in rapidity (value HIJING is 1):
      READ (24, *) ipop
      if(ipop.eq.1) IHPR2(11)=3
c     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
      READ (24, *) PARJ(5)
c     shadowing flag in HIJING:
      READ (24, *) IHPR2(6)
c     quenching flag in HIJING:
      READ (24, *) IHPR2(4)
c     quenching rate when quenching flag is on (=1.0 GeV/fm):
      READ (24, *) HIPR1(14)
c     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. 
      READ (24, *) HIPR1(8)
c     read initialization value for parton cascade:
      READ (24, *) xmu
      READ (24, *) izpc
      READ (24, *) alpha
c     quark coalescence radii in momentum and space for string melting:
      READ (24, *) dpcoal
      READ (24, *) drcoal
c     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
      READ (24, *) ihjsed
c     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
      READ (24, *) nseed
      READ (24, *) iseedp
      READ (24, *) iksdcy
      READ (24, *) iphidcy
      READ (24, *) ipi0dcy
c     flag for OSCAR output for final partons and hadrons:
      READ (24, *) ioscar
clin-5/2008     flag for perturbative treatment of deuterons:
      READ (24, *) idpert
      READ (24, *) npertd
      READ (24, *) idxsec
clin-6/2009 To select events that have at least 1 high-Pt minijet parton:
      READ (24, *) pttrig
      READ (24, *) maxmiss
      READ (24, *) IHPR2(2)
      READ (24, *) IHPR2(5)
clin-6/2009 To embed a back-to-back q/qbar pair into each event:
      READ (24, *) iembed
      READ (24, *) pxqembd, pyqembd
      READ (24, *) xembd, yembd
      READ (24, *) nsembd,psembd,tmaxembd
clin-7/2009 Allow modification of nuclear shadowing:
      READ (24, *) ishadow
      READ (24, *) dshadow
      READ (24, *) iphirp
c
      CLOSE (24)
 111  format(a8)
clin-6/2009 ctest off turn on jet triggering:
c      IHPR2(3)=1
c     Trigger Pt of high-pt jets in HIJING:
c      HIPR1(10)=7.
c
      if(isoft.eq.1) then
         amptvn = '1.26t9b (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.26t9b (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
      WRITE(12,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#          Version ',a25,                  '     #'/1X,
     &10X,'#               12/25/2018                       #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')
c     when ihjsed=11: use environment variable at run time for HIJING nseed:
      if(ihjsed.eq.11) then
         PRINT *,
     1 '# Read in NSEED in HIJING at run time (e.g. 20030819):'
      endif
      READ (*, *) nseedr
      if(ihjsed.eq.11) then
         nseed=nseedr
cfast: an optional 2nd number is the ZPC seed; it must be > 0 (ran1 in
c     zpc.f misbehaves for negative seeds), else iseedp of input.ampt is kept
         READ (*, *, IOSTAT=ios_seed) iseedpr
         if(ios_seed.eq.0) then
            if(iseedpr.gt.0) then
               iseedp = iseedpr
            else
               PRINT *, '# ZPC seed must be > 0, kept iseedp =', iseedp
            endif
         endif
      endif
cfast: NSEED=2*NSEED+1 below must not overflow a 32-bit integer
      if(nseed.lt.0.or.nseed.gt.1073741823) then
         PRINT *, '# NSEED out of 0..1073741823, reduced:', nseed
         nseed=mod(nseed,1073741823)
         if(nseed.lt.0) nseed=nseed+1073741823
      endif
      if(ihjsed.eq.11) then
         PRINT *, '#   read in: ', nseed, ' ZPC seed: ', iseedp
         WRITE(12,*) '# Read in NSEED in HIJING at run time:',nseed
         WRITE(12,*) '# ZPC seed (iseedp):',iseedp
      endif
      CALL root_set_seeds(nseed, iseedp)
      CLOSE(12)
clin-5/2015 an odd number is needed for the random number generator:
c      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
      NSEED=2*NSEED+1
c     9/26/03 random number generator for f77 compiler:
      CALL SRAND(NSEED)
c
c.....turn on warning messages in nohup.out when an event is repeated:
      IHPR2(10) = 1
c     string formation time:
      ARPAR1(1) = 0.7
c     smearp is the smearing halfwidth on parton z0, 
c     set to 0 for now to avoid overflow in eta.
c     smearh is the smearing halfwidth on string production point z0.
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
      nevent=NEVNT
c
c     AMPT momentum and space info at freezeout:
cfast: the ROOT tree is always written; ana/ampt.dat with AMPT_ENABLE_TXT=1
c     and ana/zpc.dat with AMPT_ENABLE_ZPC=1 only:
      itxtdt = iroot_is_txt_enabled()
      izpcdt = iroot_is_zpc_enabled()
      CALL root_init()
      if (itxtdt.eq.1) then
         OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      endif
      if (izpcdt.eq.1) then
         OPEN (14, FILE = 'ana/zpc.dat', STATUS = 'UNKNOWN')
      endif
ctest off for resonance (phi, K*) studies:
c      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
c      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
      CALL INIZPC
clin-5/2009 ctest off:
c      call flowp(0)
c      call flowh0(NEVNT,0)
c      call iniflw(NEVNT,0)
c      call frztm(NEVNT,0)
c
       DO 2000 J = 1, NEVNT
          IAEVT = J
          DO 1000 K = 1, NUM
             IARUN = K
             IF (IAEVT .EQ. NEVNT .AND. IARUN .EQ. NUM) THEN
                IOUT = 1
             END IF
             PRINT *, ' EVENT ', J, ', RUN ', K
             imiss=0
 100         CALL HIJING(FRAME, BMIN, BMAX)
             IAINT2(1) = NATT

clin-6/2009 ctest off
           if(J.eq.-2) then 
              write(98,*) HIPR1
              write(98,*) ' '
              write(98,*) IHPR2
              write(98,*) ' '
              write(98,*) (HINT1(i),i=1,20)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=21,40)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=41,60)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=61,80)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=81,100)
              write(98,*) ' '
              write(98,*) IHNT2
           endif

c     evaluate Npart (from primary NN collisions) for both proj and targ:
             call getnp
cfast: nucleons and collision numbers of this HIJING event for the ROOT tree
             CALL pass_nucleons_to_root
c     switch for final parton fragmentation:
             IF (IHPR2(20) .EQ. 0) GOTO 2000
c     In the unlikely case of no interaction (even after loop of 20 in HIJING),
c     still repeat the event to get an interaction 
c     (this may have an additional "trigger" effect):
             if(NATT.eq.0) then
                imiss=imiss+1
                if(imiss.le.20) then
                   write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
                   goto 100
                else
                   write(6,*) 'missed event: natt=0,j=',j
                   goto 2000
                endif
             endif
c.....ART initialization and run
             CALL ARINI
             CALL ARINI2(K)
 1000     CONTINUE
c
          CALL ARTAN1
clin-9/2012 Analysis is not used:
c          CALL HJANA3
          CALL ARTMN
clin-9/2012 Analysis is not used:
c          CALL HJANA4
          CALL ARTAN2
 2000  CONTINUE
c
       CALL ARTOUT(NEVNT)
clin-5/2009 ctest off:
c       call flowh0(NEVNT,2)
c       call flowp(2)
c       call iniflw(NEVNT,2)
c       call frztm(NEVNT,2)
c
       CALL root_close()
       STOP
       END
c     FYI: taken file unit numbers are 10-88, 91-93; 
c     so free file unit numbers are 1-4,7-9,89,97-99.

cfast: passes the HIJING nucleons of this event to the ROOT tree, as in
c     ana/npart-xy.dat: projectile nucleons first, then target nucleons,
c     at x,y shifted by +-b/2 along the reaction plane, Stat = NFP(i,5) or
c     NFT(i,5) of HIJING (0: no interaction; target values not negated).
c     Also passes NCOLT of HIJING and the spectator protons and neutrons
c     counted by getnp (so that Npart = A - spectators, as npart1/npart2).
      SUBROUTINE pass_nucleons_to_root
      COMMON/hjcrdn/YP(3,300),YT(3,300)
      COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      common /lastt/itimeh,bimp
      common /phiHJ/iphirp,phiRP
      common /hjncol/ ncolhj
      common /fastsp/ nsp1p, nsp1n, nsp2p, nsp2n
      INTEGER statnucl(600)
      REAL xnucl(600), ynucl(600), znucl(600)
      SAVE

      na = min(IHNT2(1), 300)
      nb = min(IHNT2(3), 300)
      nab = na + nb
      do jp = 1, na
         xnucl(jp) = YP(1,jp)+0.5*bimp*cos(phiRP)
         ynucl(jp) = YP(2,jp)+0.5*bimp*sin(phiRP)
         znucl(jp) = YP(3,jp)
         statnucl(jp) = NFP(jp,5)
      enddo
      do jt = 1, nb
         i = na + jt
         xnucl(i) = YT(1,jt)-0.5*bimp*cos(phiRP)
         ynucl(i) = YT(2,jt)-0.5*bimp*sin(phiRP)
         znucl(i) = YT(3,jt)
         statnucl(i) = NFT(jt,5)
      enddo

      CALL root_set_nucleons(na, nb, statnucl, xnucl, ynucl, znucl)
      CALL root_set_header(ncolhj, nsp1n, nsp1p, nsp2n, nsp2p)
      return
      end
