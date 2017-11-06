MODULE ENERGIA
CONTAINS
SUBROUTINE ENER_CONF(N,y,BOXL,V)
  USE INICIAL
IMPLICIT REAL*4 (A-H,O-Z)
!==============================================================================
!SUBRUTINA PARA CALCULAR LA ENERGIA DE LA CONFIGURACION
!MODELO DE POTENCIAL DE INTERACCION: DISCOS DUROS
!==============================================================================
!DECLARACION DE VARIABLES
INTEGER               ::N,i,j
REAL,ALLOCATABLE      ::y(:,:)

V = 0.0
Vij = 0.0
DO i=1,N-1
  rix = y(i,1)
  riy = y(i,2)
  DO j=i+1,N
    rijx = rix - y(j,1)
    rijy = riy - y(j,2)

!CONDICION DE IMAGEN MINIMA
    rijx = rijx - BOXL*(ANINT(rijx/BOXL))
    rijy = rijy - BOXL*(ANINT(rijy/BOXL))
!===============================================================================

    rr = SQRT(rijx**2 + rijy**2)
    IF (rr < (BOXL/2.0)) THEN

!MODELO DE POTENCIAL DE INTERACCION (HD)
      IF (rr < 1.0) THEN
        Vij= 1.0E+13
      ELSE
        Vij=0.0
      END IF

      V = V + Vij
    END IF
  END DO
END DO

END SUBROUTINE ENER_CONF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ENER_PART(N,l,xl,yl,x,BOXL,V) !lyr
  USE INICIAL
IMPLICIT REAL*4 (A-H,O-Z)
!==============================================================================
!SUBRUTINA PARA CALCULAR LA ENERGIA DE LA L-ESIMA PARTICULA
!MODELO DE POTENCIAL DE INTERACCION: DISCOS DUROS
!==============================================================================
!DECLARACION DE VARIABLES

INTEGER               ::N,j,l
REAL,ALLOCATABLE      ::x(:,:)

V = 0.0
rix = xl !lyr
riy = yl !lyr
DO j=1,N
  IF(l.NE.j)THEN
    rijx = rix - x(j,1)
    rijy = riy - x(j,2)

!CONDICION DE IMAGEN MINIMA
    rijx = rijx - BOXL*(ANINT(rijx/BOXL))
    rijy = rijy - BOXL*(ANINT(rijy/BOXL))
!===============================================================================
!MODELO DE POTENCIAL DE INTERACCION (HD)
    rr = SQRT(rijx**2 + rijy**2)

    IF (rr < (BOXL/2.0)) THEN

      !MODELO DE POTENCIAL DE INTERACCION (HD)
            IF (rr < 1.0) THEN
              Vij= 1.0E+13
            ELSE
              Vij=0.0
            END IF

      V = V + Vij
    END IF

  END IF
END DO
END SUBROUTINE ENER_PART
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE ENERGIA
