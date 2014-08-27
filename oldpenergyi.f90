	MODULE OLDPENERGYI
	USE COMMONS
	USE CALCENERGYI
	USE CALCUNITVECI
	CONTAINS
	
	SUBROUTINE OLDPENERGY(O, JB, IND)
	IMPLICIT NONE

	INTEGER :: O, JB, IND, J, P, Q
	DOUBLE PRECISION :: JENERGY, R2, DX, DY, DZ, SI(3), SJ(3), RIJ(3), RI(3), RJ(3)
	OPENERGY = 0.D0
	P = 3*O
	Q = 4*O
	

	! CALCULAE THE ENERGY OF ALL THE REAL SPACE INEREATCTION
	IF (IND == 1) THEN
		DO J = JB, NMOL
			IF (3*J .NE. P) THEN
				DX = COORDS(P-2) - COORDS(3*J-2)
				DY = COORDS(P-1) - COORDS(3*J-1)
				DZ = COORDS(P) - COORDS(3*J)
				R2 = DX*DX + DY*DY + DZ*DZ
				RI(1) = COORDS(P-2)
				RI(2) = COORDS(P-1)
				RI(3) = COORDS(P)
				RJ(1) = COORDS(3*J-2)
				RJ(2) = COORDS(3*J-1)
				RJ(3) = COORDS(3*J)
				RIJ = RI - RJ
					IF (SQRT(R2)<=(L/2))
						CALL CALCUNITVEC(Q, J, SI(3), SJ(3))
						CALL CALCENERGY(R2, JENERGY, SI(3), SJ(3), RIJ(3))
					ELSE (SQRT(R2)>(L/2))
						JENERGY = 0
					END IF
				OPENERGY = OPENERGY + JENERGY
				!ACCOUNT FOR PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
				IF ((RI(1) - L) > -L) THEN
					DX = (COORDS(P-2) - L) - COORDS(3*J-2)
					DY = COORDS(P-1) - COORDS(3*J-1)
 					DZ = COORDS(P) - COORDS(3*J)
					R2 = DX*DX + DY*DY + DZ*DZ
                                	RI(1) = COORDS(P-2)
                                	RI(2) = COORDS(P-1)
                                	RI(3) = COORDS(P)
                                	RJ(1) = COORDS(3*J-2)
                                	RJ(2) = COORDS(3*J-1)
                                	RJ(3) = COORDS(3*J)
                                	RIJ = RI - RJ
						IF (R2<=(L/2))
                                                	CALL CALCUNITVEC(Q, J, SI(3), SJ(3))
                                                	CALL CALCENERGY(R2, JENERGY, SI(3), SJ(3), RIJ(3))
                                        	ELSE (R2>(L/2))
                                               		JENERGY = 0
                                        	END IF
					OPENERGY = OPENERGY + JENERGY
				END IF
				IF ((RI(1) + L) < L) THEN
                                        DX = (COORDS(P-2) + L) - COORDS(3*J-2)
                                        DY = COORDS(P-1) - COORDS(3*J-1)
                                        DZ = COORDS(P) - COORDS(3*J)
                                        R2 = DX*DX + DY*DY + DZ*DZ
                                        RI(1) = COORDS(P-2)
                                        RI(2) = COORDS(P-1)
                                        RI(3) = COORDS(P)
                                        RJ(1) = COORDS(3*J-2)
                                        RJ(2) = COORDS(3*J-1)
                                        RJ(3) = COORDS(3*J)
                                        RIJ = RI - RJ
                                                IF (R2<=(L/2))                                                                                                  
                                                        CALL CALCUNITVEC(Q, J, SI(3), SJ(3))                                                                    
                                                        CALL CALCENERGY(R2, JENERGY, SI(3), SJ(3), RIJ(3))                                                      
                                                ELSE (R2>(L/2))                                                                                                 
                                                        JENERGY = 0
                                                END IF                                                                                                          
                                        OPENERGY = OPENERGY + JENERGY
				ENDIF
				!ACCOUNT FOR PERODIC BOUNDARY CONDITIONS IN Y DIRECTION
				DX = COORDS(P-2) - COORDS(3*J-2)
                                DY = (COORDS(P-1) - L) - COORDS(3*J-1)
                                DZ = COORDS(P) - COORDS(3*J)
                                R2 = DX*DX + DY*DY + DZ*DZ
                                RI(1) = COORDS(P-2)
                                RI(2) = COORDS(P-1)
                                RI(3) = COORDS(P)
                                RJ(1) = COORDS(3*J-2)
                                RJ(2) = COORDS(3*J-1)
                                RJ(3) = COORDS(3*J)
                                RIJ = RI - RJ
                                        IF (R2<=(L/2))
                                                CALL CALCUNITVEC(Q, J, SI(3), SJ(3))
                                                CALL CALCENERGY(R2, JENERGY, SI(3), SJ(3), RIJ(3))
                                        ELSE (R2>(L/2))
                                                JENERGY = 0
                                        END IF
				OPENERGY = OPENERGY + JENERGY
				!ACCOUNT FOR PERIODIC BOUNDARY CONDITIONS IN Z DIRECTION
				DX = COORDS(P-2) - COORDS(3*J-2)
                                DY = COORDS(P-1) - COORDS(3*J-1)
                                DZ = (COORDS(P) - L) - COORDS(3*J)
                                R2 = DX*DX + DY*DY + DZ*DZ
                                RI(1) = COORDS(P-2)
                                RI(2) = COORDS(P-1)
                                RI(3) = COORDS(P)
                                RJ(1) = COORDS(3*J-2)
                                RJ(2) = COORDS(3*J-1)
                                RJ(3) = COORDS(3*J)
                                RIJ = RI - RJ
                                        IF (R2<=(L/2))
                                                CALL CALCUNITVEC(Q, J, SI(3), SJ(3))
                                                CALL CALCENERGY(R2, JENERGY, SI(3), SJ(3), RIJ(3))
                                        ELSE (R2>(L/2))
                                                JENERGY = 0
                                        END IF
				OPENERGY = OPENERGY + JENERGY
			END IF
		END DO
	
	!CALCULATE THE ENERGY OF ALL THE FOURIER SPACE TERMS
	END IF
	RETURN
	END SUBROUTINE
	END MODULE
