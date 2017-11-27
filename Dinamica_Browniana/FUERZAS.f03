MODULE FRZAS
CONTAINS
SUBROUTINE FUERZAS(T,RCUT,BOXL,N,ISTEP,ENPOT,x,f,VIR)
IMPLICIT REAL (A-H, O-Z)
REAL,ALLOCATABLE::x(:,:),f(:,:)

ENPOT = 0.0
f=0.0
VIR = 0.0

DO i=1,N-1
  fxi = f(i,1)
  fyi = f(i,2)
  fzi = f(i,3)

  DO j=i+1,N
    rijx = x(i,1)-x(j,1)
    rijy = x(i,2)-x(j,2)
    rijz = x(i,3)-x(j,3)

    !CONDICION DE IMAGEN MINIMA
    rijx = rijx - BOXL*(ANINT(rijx/BOXL))
    rijy = rijy - BOXL*(ANINT(rijy/BOXL))
    rijz = rijz - BOXL*(ANINT(rijz/BOXL))

    rr = SQRT(rijx**2 + rijy**2 + rijz**2)

    !VERIFICANDO TRASLAPES
    IF (rr < 1.0)THEN
      PRINT*, "TRASLAPE",i,j
    END IF

    !INICIA LA ESPECIFICACION DEL MODELO DE POTENCIAL DEL CUAL SE DERIVA LA FUERZA
    !DE INTERACCION ENTRE LAS PARTICULAS

    IF (rr < RCUT)THEN
      !MODELO DE POTENCIAL:YUKAWA
      !U = EXP(-ZK*rr)
      !U2 = AA*U*(ZK*rr+1.0)/(rr**3)
      !ENPOT = (AA*U)/rr+ENPOT

      !MODELO DE POTENCIAL: GAUSSIANO ATRACTIVO
      U = -(1.0/T)*EXP(-rr**2) 
      U2 = 2.0*U
      ENPOT = U + ENPOT
      VIR =  + VIR

      fijx = rijx*U2
      fijy = rijy*U2
      fijz = rijz*U2

      fxi = fxi+fijx
      fyi = fyi+fijy
      fzi = fzi+fijz

      f(j,1) = f(j,1)-fijx
      f(j,2) = f(j,2)-fijy
      f(j,3) = f(j,3)-fijz
    END IF

  END DO
  f(i,1)=fxi
  f(i,2)=fyi
  f(i,3)=fzi
END DO
VIR = VIR(3.0*REAL(N))

END SUBROUTINE FUERZAS
END MODULE FRZAS
