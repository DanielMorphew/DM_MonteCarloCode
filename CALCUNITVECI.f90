	MODULE CALCUNITVECI
	USE COMMONS
	CONTAINS

	SUBROUTINE CALCUNITVEC(Q, J, SI(3), SJ(3), QRIJ(3), DQUART
	IMPLICIT NONE

	INTEGER :: Q, J
	DOUBLE PRECISION :: SI(3), SJ(3), QRIJ(3), QROTMAT(3:3), JROTMAT(3:3), STNDV(3), DQUART

	QROTMAT(1:1) = 1 - 2*(QUARTS(Q-1)*QUARTS(Q-1)) - 2*(QUARTS(Q)*QUARTS(Q))
	QROTMAT(1:2) = 2*(QUARTS(Q-2)*QUARTS(Q-1)) - 2*(QUARTS(Q-3)*QUARTS(Q)
	QROTMAT(1:3) = 2*(QUARTS(Q-2)*QUARTS(Q) + 2*(QUARTS(Q-3)*QUARTS(Q-2))
	QROTMAT(2:1) = 2*(QUARTS(Q-1)*QUARTS(Q-2)) + 2*(QUARTS(Q-3)*QUARTS(Q))
	QROTMAT(2:2) = 1 - 2*(QUARTS(Q)*QUARTS(Q)) - 2*(QUARTS(Q-2)*QUARTS(Q-2))
	QROTMAT(2:3) = 2*(QUARTS(Q-1)*QUARTS(Q)) - 2*(QUARTS(Q-3)*QUARTS(Q-2))
	QROTMAT(3:1) = 2*(QUARTS(Q)*QUARTS(Q-2)) - 2*(QUARTS(Q-3)*QUARTS(Q-1))
	QROTMAT(3:2) = 2*(QUARTS(Q)*QUARTS(Q-1)) + 2*(QUARTS(Q-3)*QUARTS(Q-2))
	QROTMAT(3:3) = 1 - 2*(QUARTS(Q-2)*QUARTS(Q-2)) - 2*(QUARTS(Q-1)*QUARTS(Q-1))

        JROTMAT(1:1) = 1 - 2*(QUARTS(4*J-1)*QUARTS(4*J-1)) - 2*(QUARTS(4*J)*QUARTS(4*J))
        JROTMAT(1:2) = 2*(QUARTS(4*J-2)*QUARTS(4*J-1)) - 2*(QUARTS(4*J-3)*QUARTS(4*J))
        JROTMAT(1:3) = 2*(QUARTS(4*J-2)*QUARTS(4*J)) + 2*(QUARTS(4*J-3)*QUARTS(4*J-2))
        JROTMAT(2:1) = 2*(QUARTS(4*J-1)*QUARTS(4*J-2)) + 2*(QUARTS(4*J-3)*QUARTS(4*J))
        JROTMAT(2:2) = 1 - 2*(QUARTS(4*J)*QUARTS(4*J)) - 2*(QUARTS(4*J-2)*QUARTS(4*J-2))
        JROTMAT(2:3) = 2*(QUARTS(4*J-1)*QUARTS(4*J)) - 2*(QUARTS(4*J-3)*QUARTS(4*J-2))
        JROTMAT(3:1) = 2*(QUARTS(4*J)*QUARTS(4*J-2)) - 2*(QUARTS(4*J-3)*QUARTS(4*J-1))
        JROTMAT(3:2) = 2*(QUARTS(4*J)*QUARTS(4*J-1)) + 2*(QUARTS(4*J-3)*QUARTS(4*J-2))
        JROTMAT(3:3) = 1 - 2*(QUARTS(4*J-2)*QUARTS(4*-2)) - 2*(QUARTS(4*J-1)*QUARTS(4*J-1))

	STNDV(1) = 0
	STNDV(2) = 0
	STNDV(3) = 1

	SI = MATMUL(QROTMAT, STNDV)
	SJ = MATMUL(JROTMAT, STNDV)


	RETURN
	END SUBROUTINE
	END MODULE
