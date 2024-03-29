	PROGRAM GEN_COORDSDM
	IMPLICIT NONE
	! GENERATES RANDOM CONFIGURATION FOR A RIGID BODY SYSTEM
	!PARAMETERS: NCONFIG, NMOL, RADIUS MUST BE SET

	INTEGER, PARAMETER :: NCONFIG = 1, NMOL = 14
	DOUBLE PRECISION, PARAMETER :: RADIUS = 10.D0
	INTEGER :: I, J, K, NATOMS
	DOUBLE PRECISION :: COORDS(NCONFIG,(7*NMOL)), RANF, DUMMY, Q(4) 
	CHARACTER (LEN = 10) :: CI, FILENAME

	NATOMS = 2*NMOL !ACCOUNTS FOR NEED TO ICLUDE ORIENTATION DESCRIPTORS IN OUTPUT FILE

	DO I = 1, NCONFIG
		
		WRITE(CI,'(I2)') I
		CI = ADJUSTL(CI)  
		FILENAME = 'COORDS'//TRIM(CI)
		OPEN( UNIT = 21, FILE = FILENAME, STATUS = 'UNKNOWN' ) !OPENS FILE FOR RANDOM CONFIGURATION DATA TO BE KEPT WITHIN

		DO J = 1, 3*NMOL !*3 IN ORDER TO GENERATE X, Y, Z COORDINATES
			COORDS(I,J) = RADIUS*(RANF(DUMMY)-0.5D0)*2.0D0 !USE RANF FUNCTION TO GENERATE RANDOM COORDINATE WITHIN RADIUS SELECTED
		END DO

		DO J= 1, NMOL
			CALL RANDQTN(Q) !SUBROUTINE THAT DEFINES RANDOM ORIENTATION OF EACH PARTICLE
			K = 3*NMOL + 4*J
			COORDS(I, K-3:K) = Q(1:4)
		END DO

		DO J = 1, NMOL 
			K = 3*J
			WRITE(21,*) COORDS(I,K-2), COORDS(I,K-1), COORDS(I,K) !RANDOM POSITION OF PARTICLE INTO OPEN OUTPUT FILE
		END DO
		DO J = NMOL+1, NATOMS
			K = 4*J
			WRITE(21,*) COORDS(I,K-3), COORDS(I,K-2), COORDS(I,K-1), COORDS(I,K) !RANDOM QUARTERNION INTO OPEN OUTPUT FILE
		END DO

	END DO

	END PROGRAM GEN_COORDSDM

! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	SUBROUTINE RANDQTN(Q)
	IMPLICIT NONE 
	! GENERATES UNIFORMLY RANDOM QUARTERNION

	DOUBLE PRECISION :: Q(4), RANF, DUMMY, FCTR, XIASQ, XIBSQ, XI1, XI2, XI3, XI4
	
	XIASQ = 1.0
	XIBSQ = 1.0

        ! ITERATICE LOOP
	DO WHILE (XIASQ >= 1.D0)
		XI1  = RANF(DUMMY)*2.D0 - 1.D0
		XI2  = RANF(DUMMY)*2.D0 - 1.D0
		XIASQ = XI1*XI1 + XI2*XI2
	ENDDO
      
	DO WHILE (XIBSQ >= 1.D0)
		XI3  = RANF(DUMMY)*2.D0 - 1.D0
		XI4  = RANF(DUMMY)*2.D0 - 1.D0
		XIBSQ = XI3*XI3 + XI4*XI4
	ENDDO
      
	FCTR = SQRT((1.D0-XIASQ)/XIBSQ)
	Q(1) = XI1
	Q(2) = XI2
	Q(3) = XI3*FCTR
	Q(4) = XI4*FCTR

	END SUBROUTINE RANDQTN


! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	DOUBLE PRECISION FUNCTION RANF(DUMMY)

	!FUNCTION TO GENERATE A UNIFORMLY RANDOM NUMBER BETWEEN 0 AND 1

	INTEGER, PARAMETER :: L = 1029, C = 221591, M = 1048576
	INTEGER ::          SEED
	DOUBLE PRECISION :: DUMMY
	SAVE             SEED
	DATA             SEED / 0 /

	SEED = MOD(SEED * L + C, M)
	RANF = DFLOAT(SEED) / DFLOAT(M)

	END FUNCTION RANF

