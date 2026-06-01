c     GILLESPIE SIMULATION
c	Pulsatile Stimulus Signal
c.............................................................................
c     Adipogenic Model 

      implicit real*8 (a-h,q-z)

      parameter(nu=20)
      dimension a_nu(nu)
      character(len=48) :: filename

	call cpu_time(start)
	write(*,*) 'no. of cell   adipocyte no'
c	open(1002,file='pulsatile.out',status='unknown')


c......................................................................

c     -----------------------------------------------------------------
c     STEP 0: INITIALIZATION
c     -----------------------------------------------------------------
        nk=500   !Stochastic simulation run, analogue of cell count
        l100=150 ! random number seed
c.....................................................................

    	nadipo = 0 
c......................................................................

       do 200 jj=1,nk
c      call random_seed
        l100=l100+10
        iseed=l100
        useed=dustar(iseed)

	nm = jj + 0

       write(filename,'("pulse_",I4.4,".out")') nm
       open(1001,file=filename,status='unknown')
c.....................................................................
c      Initial conditions 
c......................................................................

c       Initial values , Y0

          nmPPARG = 0
          nPPARGi = 0
          nPPARGa = 0
	   nmCEBPA = 0
	   nCEBPA = 0
	   nmFABP4 = 0
	   nFABP4 = 0

c.................................................................
c......................................................................
	ti=(-1440.0)  !-1 day started for synchronization as well as to get correct phase
        t_f=(5760.0) ! 4 day
c
        nrxn=1 ! rpd , reaction per dot
c
c     -----------------------------------------------------------------
c     STEP 1: CALCULATION OF a_nu=h_nu*c_nu AND a0
c     -----------------------------------------------------------------
c.....................................................................
c      Parameter values for the problem

c.....................................................................
  10   pulse = 1.9  ! pulse strength
	 dur = 12*60 ! duration of pulse
	 day = 24*60 ! duration of total cycle of a day	

	
	 if (ti.ge.0) then ! ensures the pulse will be given at day 0
		if (mod(ti,day).lt.dur) then
		        stim = pulse
	 	else
			stim = 0
	 	end if
	 else
		stim = 0
	 end if
	 sfpp = 1000 ! scaling factor for PPARG network
	 sf1 = 0.1 !sf1=scaling of mRNAs 
	 sfP = 1.2
	 sfC = 0.25
	 sfF = 0.2
	 sfP1 = 1.2
	 sfC1 = 4
	 sfF1 = 3
	 ay = 1
	 ak = 0.00054
	 an = 1.4
	 aksP = 0.619*sfpp*sf1*sfP*sfP1*ay
	 akP1 = 0.99*sfpp*sf1*sfP*sfP1*ay
	 akmP1 = 6.5*sfpp*sfC
	 akP2 = 1.21*sfpp*sf1*sfP*sfP1*ay
	 akmP2 = 12.3*sfpp*sfF
	 akdmP = 0.9*ay
	 aktP = (0.5*ay)/(sfP1*sf1)
	 ak0 = 140*ay
	 akact = 20*ay*exp((ak*dur)**an)
	 akinact = 4*ay
	 akdP = 0.7*ay
	 bmCEBPA = 0.016*sfpp*sf1*sfC*sfC1*ay
	 akC = 39.36*sfpp*sf1*sfC*sfC1*ay
	 akmC = 10*sfpp*sfP
	 akdmC = 0.1*ay
	 aktC = (0.1*ay)/(sfC1*sf1)
	 akdC = 0.07*ay
	 bmFABP4 = 0.01*sfpp*sf1*sfF*sfF1*ay
	 akF = 6.9*sfpp*sf1*sfF*sfF1*ay
	 akmF = 20*sfpp*sfP
	 akdmF = 0.006*ay
	 aktF = (0.15*ay)/(sfF1*sf1)
	 akdF = 0.11*ay 
	 
c.....................................................................
c      Reactions for the problem   	
	
	 a_nu(1)= aksP
	 a_nu(2)= (akP1*nCEBPA**2)/(akmP1**2+nCEBPA**2)
	 a_nu(3)= (akP2*nFABP4**2)/(akmP2**2+nFABP4**2)
	 a_nu(4)= akdmP*nmPPARG
	 a_nu(5)= aktP*nmPPARG
	 a_nu(6)= akinact*nPPARGa
	 a_nu(7)= ak0*nPPARGi
	 a_nu(8)= akact*stim*nPPARGi
	 a_nu(9)= akdP*nPPARGi
	 a_nu(10)= akdP*nPPARGa
	 a_nu(11)= bmCEBPA
	 a_nu(12)= (akC*nPPARGa**2)/(akmC**2+nPPARGa**2)
	 a_nu(13)= akdmC*nmCEBPA
	 a_nu(14)= aktC*nmCEBPA
	 a_nu(15)= akdC*nCEBPA
	 a_nu(16)= bmFABP4
	 a_nu(17)= (akF*nPPARGa**2)/(akmF**2+nPPARGa**2)
	 a_nu(18)= akdmF*nmFABP4
	 a_nu(19)= aktF*nmFABP4
	 a_nu(20)= akdF*nFABP4
	 

c.........................................................................       
      sum1=0.0
       do 85 i=1,nu
         sum1=sum1+a_nu(i)
 85    continue
       a0=sum1

c     -----------------------------------------------------------------
c     STEP 2: CALCULATION OF tau AND mu
c     -----------------------------------------------------------------

c     tau is the infinitesimal time (continuous) during which a 
c     particular reaction takes place
c     mu determines the specific reaction channel

c      r1=rand(idum)
      r1=duni()
       if(r1.eq.0.0)then
         r1=0.0001
         endif
c      call random_number(harvest=r1)
      tau=log(1.0/r1)/a0

c      r2=rand(idum)
       r2=duni()
c      call random_number(harvest=r2)
      r2a0=r2*a0

      sum2=0.0
      do 30 j=1,nu
         mu=j
         sum2=sum2+a_nu(j)

         if(sum2.ge.r2a0)then
		
		go to 100
	end if
 30   continue

c     -----------------------------------------------------------------
c     STEP 3: SAMPLING
c     -----------------------------------------------------------------
c......................................................................
c     this step updates the population of each reacting species
c......................................................................
 100	 go to(101,102,103,104,105,106,107,108,109,110,111,112,113,114,
     $115,116,117,118,119,120),mu



 101	 nmPPARG = nmPPARG+1
	 go to 50
 102  nmPPARG = nmPPARG+1
      go to 50
 103  nmPPARG = nmPPARG+1
      go to 50
 104  nmPPARG = nmPPARG-1
	 go to 50
 105	 nPPARGi = nPPARGi+1
	 go to 50
 106	 nPPARGi = nPPARGi+1
	 nPPARGa = nPPARGa-1
	 go to 50
 107	 nPPARGa = nPPARGa+1
	 nPPARGi = nPPARGi-1
	 go to 50
 108	 nPPARGi = nPPARGi-1
	 nPPARGa = nPPARGa+1
	 go to 50
 109	 nPPARGi = nPPARGi-1
	 go to 50
 110	 nPPARGa = nPPARGa-1
	 go to 50
 111	 nmCEBPA = nmCEBPA+1
	 go to 50
 112	 nmCEBPA = nmCEBPA+1
	 go to 50
 113	 nmCEBPA = nmCEBPA-1
	 go to 50
 114 	 nCEBPA = nCEBPA+1
	 go to 50
 115	 nCEBPA = nCEBPA-1
	 go to 50
 116	 nmFABP4 = nmFABP4+1
	 go to 50
 117	 nmFABP4 = nmFABP4+1
	 go to 50
 118	 nmFABP4 = nmFABP4-1
	 go to 50
 119	 nFABP4 = nFABP4+1
	 go to 50
 120	 nFABP4 = nFABP4-1
	 go to 50


c---------------------------------------------------------------
 50     nrxn=nrxn+1
        ti = ti+tau

c..................................................................
c....................................................................
c       To draw stochastic profile at any instance
c......................................................................
        if (mod(nrxn,10000).eq.0) then
                 write(1001,*) ti,nmPPARG,nPPARGi,nPPARGa,nmCEBPA,
     #nCEBPA,nmFABP4,nFABP4,stim
         end if

	
c.......................................................................
c         Time loop continues
c......................................................................



	if(ti.lt.t_f)go to 10
	if (nPPARGa.ge.1137) then
		nadipo = nadipo + 1
c		write(*,*) nadipo
		go to 201
	end if
 
 201	write(*,*) jj,nadipo,'Done'
c	write(1002,*) jj, nPPARGa ! cell_id  last PPARGa exp
 200    continue
c       end do
c........................................................................
    
       write(*,*) 'Done'

	call cpu_time(finish)
	runtime = finish-start
	write(*,*) runtime/60, 'minutes' 
       stop
       end

c......................................................

c..........................................................................
c      Subroutine for random no generator by Lane Watson
c.........................................................................
            DOUBLE PRECISION FUNCTION DUNI()
c            FUNCTION DUNI()
C***BEGIN PROLOGUE  DUNI
C***DATE WRITTEN   880714 (YYMMDD)
C***REVISION DATE  880714 (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  THIS ROUTINE GENERATES DOUBLE PRECISION UNIFORM
C             RANDOM NUMBERS ON [0,1)
C***DESCRIPTION
C        COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C
C       USAGE:
C              TO INITIALIZE THE GENERATOR
C                   USEED = DUSTAR(ISEED)
C               WHERE: ISEED IS ANY NONZERO INTEGER
C                  WILL RETURN FLOATING POINT VALUE OF ISEED.
C
C               SUBSEQUENTLY
C                       U = DUNI()
C                  WILL RETURN A REAL UNIFORM ON [0,1)
C
C                ONE INITIALIZATION IS NECESSARY, BUT ANY NUMBER OF EVALUATIONS
C                  OF DUNI IN ANY ORDER, ARE ALLOWED.
C
C           NOTE: DEPENDING UPON THE VALUE OF K (SEE BELOW), THE OUTPUT
C                       OF DUNI MAY DIFFER FROM ONE MACHINE TO ANOTHER.
C
C           TYPICAL USAGE:
C
C               DOUBLE PRECISION U,DUNI,DUSTAR,USEED
C               INTEGER ISEED
CC                 SET SEED
C               ISEED = 305
C               USEED = DUSTAR(ISEED)
C               DO 1 I = 1,1000
C                   U = DUNI()
C             1 CONTINUE
CC                 NOTE: IF K=47 (THE DEFAULT, SEE BELOW) THE OUTPUT VALUE OF
CC                           U WILL BE 0.812053811384E-01...
C               WRITE(*,*) U
C               END
C
C          NOTE ON PORTABILITY: USERS CAN CHOOSE TO RUN DUNI IN ITS DEFAULT
C               MODE (REQUIRING NO USER ACTION) WHICH WILL GENERATE THE SAME
C               SEQUENCE OF NUMBERS ON ANY COMPUTER SUPPORTING FLOATING POINT
C               NUMBERS WITH AT LEAST 47 BIT MANTISSAS, OR IN A MODE THAT
C               WILL GENERATE NUMBERS WITH A LONGER PERIOD ON COMPUTERS WITH
C               LARGER MANTISSAS.
C          TO EXERCISE THIS OPTION:  B E F O R E  INVOKING DUSTAR INSERT
C               THE INSTRUCTION        UBITS = DUNIB(K)      K >= 47
C               WHERE K IS THE NUMBER OF BITS IN THE MANTISSA OF YOUR FLOATING

C               POINT WORD (K=96 FOR CRAY, CYBER 205). DUNIB RETURNS THE
C               FLOATING POINT VALUE OF K THAT IT ACTUALLY USED.
C                    K INPUT AS .LE. 47, THEN UBITS=47.
C                    K INPUT AS .GT. 47, THEN UBITS=FLOAT(K)
C               IF K>47 THE SEQUENCE OF NUMBERS GENERATED BY DUNI MAY DIFFER
C               FROM ONE COMPUTER TO ANOTHER.
C
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE DUNI
      DOUBLE PRECISION CSAVE,CD,CM
      PARAMETER(
     *    CSAVE=0.9162596898123D+13/0.140737488355328D+15,
     *    CD=0.76543212345678D+14/0.140737488355328D+15,
     *    CM=0.140737488355213D+15/0.140737488355328D+15)
C                            2**47=0.140737488355328D+15
      DOUBLE PRECISION U(17),S,T,DUSTAR,C,DUNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED
C
      SAVE U,I,J,K,C
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH ISEED=305 AND K=96.
      DATA U/
     *0.471960981577884755837789724978D+00,
     *0.930323453205669578433639632431D+00,
     *0.110161790933730836587127944899D+00,
     *0.571501996273139518362638757010D-01,
     *0.402467554779738266237538503137D+00,
     *0.451181953427459489458279456915D+00,
     *0.296076152342721102174129954053D+00,
     *0.128202189325888116466879622359D-01,
     *0.314274693850973603980853259266D+00,
     *0.335521366752294932468163594171D-02,
     *0.488685045200439371607850367840D+00,
     *0.195470426865656758693860613516D+00,
     *0.864162706791773556901599326053D+00,
     *0.335505955815259203596381170316D+00,
     *0.377190200199058085469526470541D+00,
     *0.400780392114818314671676525916D+00,
     *0.374224214182207466262750307281D+00/
      DATA I,J,K,C/17,5,47,CSAVE/
C
C   BASIC GENERATOR IS FIBONACCI
C
      DUNI = U(I)-U(J)
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      U(I) = DUNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
C   SECOND GENERATOR IS CONGRUENTIAL
C
      C = C-CD
      IF(C.LT.0.0D0) C=C+CM
C
C   COMBINATION GENERATOR
C
      DUNI = DUNI-C
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      RETURN
C
      ENTRY DUSTAR(ISEED)
C
C          SET UP ...
C          CONVERT ISEED TO FOUR SMALLISH POSITIVE INTEGERS.
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED.
C
        DO 2 II = 1,17
          S = 0.0D0
          T = 0.5D0
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER K BITS, WHERE K IS DEFAULTED TO 47 BUT CAN
C               BE CHANGED BY USER CALL TO DUNIB(K)
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179)
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169)
                  IF(MOD(L1*M1,64).GE.32)S=S+T
                  T = 0.5D0*T
    3     continue
        U(II) = S
    2   continue
        DUSTAR = FLOAT(ISEED)
        RETURN
C
      ENTRY DUNIB(KK)
        IF(KK.LE.47)THEN
             K=47
        ELSE
             K=KK
        ENDIF
        DUNIB=FLOAT(K)
      END

