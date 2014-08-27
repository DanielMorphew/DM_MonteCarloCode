	MODULE COMMONS
	IMPLICIT NONE

	INTEGER :: NMCE !NUMBER OF MC CYCLE DURING EQUILLIBRIUM
	INTEGER :: NMCP !NUMBER OF MC CYCLE DURING PRODUCTION
	INTEGER :: NMCS !NUMBER OF MC CYCLES BETWEEN SAMPLES
	INTEGER :: NMOL !NUMBER OF MOLECULES
	INTEGER :: NDISP !NUMBER OF DISPLACEMENT ATTEMPTS
	INTEGER :: II
	INTEGER :: ICYCLE 
	INTEGER :: NCYCLE
	INTEGER :: IMOVE
	INTEGER :: NMOVE
	INTEGER :: ATTEMPT
	INTEGER :: NACCEPT
	INTEGER :: NATOM
	INTEGER :: FSPACELIMIT
        DOUBLE PRECISION :: TEMP !TEMPERATURE
	DOUBLE PRECISION :: RADIUS !RADIUS OF UNIT CELL
	DOUBLE PRECISION :: MAXSTEP !MAX TRANSLATIONAL STEP
	DOUBLE PRECISION :: MAXROT !MAX ROTATIONAL STEP
	DOUBLE PRECISION :: L ! BOXLENGTH
	DOUBLE PRECISION :: ALPHA
	DOUBLE PRECISION :: U !DIPOLE STRENGHT
	DOUBLE PRECISION :: ENERGY
	DOUBLE PRECISION :: NPENERGY
        DOUBLE PRECISION :: OPENERGY
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: QUARTS
	END MODULE
