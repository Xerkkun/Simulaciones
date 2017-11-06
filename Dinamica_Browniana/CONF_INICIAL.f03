MODULE INICIAL
CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE REGULAR(N,nd,BOXL1,DENS,x)
!-------------------------------------------------------------------------------
!SUBRUTINA PARA OBTENER POSICIONES REGULARES DE N PARTICULAS
!-------------------------------------------------------------------------------
IMPLICIT REAL (A-H,O-Z)
!DECLARACION DE VARIABLES
CHARACTER(len=20)::ARCHIVO
INTEGER::i,j,k,NPS,l,N,nd
REAL,ALLOCATABLE:: x(:,:)

ALLOCATE(x(N+10,nd))

!NUMERO DE PARTICULAS POR ARISTA
NPS = INT(N**(1/REAL(nd)))

!SEPARACION ENTRE PARTICULAS
DA = BOXL1/(REAL(NPS)-1.0)

!ABRIR ARCHIVO NUEVO DE DATOS
WRITE(ARCHIVO , '(A,F6.3,A)') "init", DENS, ".dat"
OPEN(10,FILE=ARCHIVO,STATUS="REPLACE",ACTION="WRITE")
!===============================================================================
!COLOCAR A LAS PARTICULAS
SELECT CASE(nd)
!ARREGLO TRIDIMENSIONAL
CASE(3)
i=1
DO WHILE (i<=N)
   !COORDENADAS Z
   DO j=0,NPS-1
      !COORDENADAS Y
      DO k=0,NPS-1
         !COORDENADAS X
         DO l=0,NPS-1
            x(i,1)=-(BOXL1/2)+(DA*l)
            x(i,2)=-(BOXL1/2)+(DA*k)
            x(i,3)=-(BOXL1/2)+(DA*j)
            i=i+1
         END DO
      END DO
   END DO
END DO
!-------------------------------------------------------------------------------
!ARREGLO BIDIMENSIONAL
CASE(2)
i=1
DO WHILE (i<=N)
      !COORDENADAS Y
      DO k=0,NPS-1
         !COORDENADAS X
         DO l=0,NPS-1
            x(i,1)=-(BOXL1/2.0)+(DA*l)
            x(i,2)=-(BOXL1/2.0)+(DA*k)
            i=i+1
         END DO
      END DO
   END DO
END SELECT
!===============================================================================
!ESCRIBIR EN ARCHIVO DE SALIDA
SELECT CASE(nd)
CASE(2)
  DO i=1,N
    WRITE(10,*)x(i,1),x(i,2)
  END DO
!-------------------------------------------------------------------------------
CASE(3)
  DO i=1,N
    WRITE(10,*)x(i,1),x(i,2),x(i,3)
  END DO
END SELECT
!===============================================================================
CLOSE(10)
END SUBROUTINE REGULAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ALEATORIA(N,nd,BOXL1,SIGMA,DENS,x)
!-------------------------------------------------------------------------------
!SUBRUTINA PARA OBTENER POSICIONES ARBITRARIAS DE N PARTICULAS
!DENTRO DE UNA CELDA DE DOS O TRES DIMENSIONES A PARTIR DE UN
!GENERADOR DE NUMEROS ALEATORIOS
!-------------------------------------------------------------------------------
IMPLICIT REAL (A-H,O-Z)
!DECLARACION DE VARIABLES
CHARACTER(len=20)::ARCHIVO
INTEGER          :: i,j,k,f,N,nd
REAL,ALLOCATABLE :: x(:,:),rnum(:,:)

ALLOCATE(x(N+10,nd))
ALLOCATE(rnum(N,nd))

!ABRIR ARCHIVO NUEVO DE DATOS
WRITE(ARCHIVO , '(A,F6.3,A)') "init", DENS, ".dat"
OPEN(10,FILE=ARCHIVO,STATUS="REPLACE",ACTION="WRITE")

!LLAMAR SEMILLA Y OBTENER NUMEROS ALEATORIOS
  CALL INIT_RANDOM_SEED()
  CALL RANDOM_NUMBER(rnum)

!COLOCAR LAS PARTICULAS
DO j=2,N
  f = 1

  DO WHILE (f==1)

  CALL RANDOM_NUMBER(rnum(j,:))

  x(j,:) = (rnum(j,:)-0.5)*BOXL1 !LAS POSICIONES DE LAS PARTICULAS
	!VERIFICAR QUE NO SE TRASLAPEN
  f = 0
  DO k=1,j-1
    r =SQRT(SUM( (x(j,:) - x(k,:))**2))
    IF(r <= SIGMA)THEN
      !PRINT*, "TRASLAPE", j, k
      f = 1
    END IF
  END DO
END DO
END DO
!===============================================================================
!ESCRIBIR EN ARCHIVO DE SALIDA
SELECT CASE(nd)
CASE(2)
  DO i=1,N
    WRITE(10,*)x(i,1),x(i,2)
  END DO
!-------------------------------------------------------------------------------
CASE(3)
  DO i=1,N
    WRITE(10,*)x(i,1),x(i,2),x(i,3)
  END DO

END SELECT
!===============================================================================
CLOSE(10)
DEALLOCATE(rnum)
END SUBROUTINE ALEATORIA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!NOTA:CORREGIR PARA UNA SEMILLA DADA POR EL USUARIO
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER,DIMENSION(:),ALLOCATABLE :: seed

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

CALL SYSTEM_CLOCK(COUNT=clock)

seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE INICIAL
