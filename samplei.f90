	MODULE SAMPLEI
	USE COMMONS
	CONTAINS
	SUBROUTINE SAMPLE
	IMPLICIT NONE

	DOUBLE PRECISION :: APE
	CHARACTER (LEN=20) :: FILENAME
	APE = 0.0D0
	FILENAME = 'averagePenergy'
	OPEN(UNIT=30, FILE=FILENAME, STATUS='UNKNOWN')

	APE = ENERGY/DBLE(NMOL)

	WRITE(30,*)'NMOL =', NMOL, 'STEP =', ICYCLE, 'AVERAGE PARTICLE ENERGY=', APE

	END SUBROUTINE
	END MODULE
	 
