
c
c
c     =====================================================
      subroutine bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      USE solver_eos, ONLY: eos_1d
c      USE properties, ONLY: pressure
      USE Interp_table, ONLY: CO2BLLT_EQUI
      USE derivees    , ONLY: CO2DER
      USE var_const, ONLY: e_const,flag_diagno,
     &                     eguess_out,vguess_in



      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension q_comp(1-mbc:maxmx+mbc, mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension mthbc(4) 
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
c
c      print*, 'bc start for t=',t
c
      go to (100,110,120,130) mthbc(1)+1
c     Rho, e fix and  u,v interpolated 
  100 continue
c
      flag_diagno = 0
      p_bc1   = 9.1d6
c      T_bc   = 310.15
      rho_bc1 = 621.475
      e_bc1   = -211.4917453d3 
c      CALL pressure(T_bc, 1.0/rho_bc, ppp)
c     Print*,'ppp',ppp
c      STOP
      DO j=1-mbc,my+mbc
 
         rho1 = q(1,j,1)
         u1   = q(1,j,2)/rho1
         v1   = q(1,j,3)/rho1
         e1   = q(1,j,4)/rho1 - 0.5*(u1*u1+v1*v1) - e_const
         CALL CO2BLLT_EQUI(p1, T1, c1, x1, a1, res1,
     &                         e1, 1.0/rho1, 0.0d0)
c
         rho2 = q(2,j,1)
         u2   = q(2,j,2)/rho2
         v2   = q(2,j,3)/rho2
         e2   = q(2,j,4)/rho2 - 0.5*(u2*u2+v2*v2) - e_const
         CALL CO2BLLT_EQUI(p2, T2, c2, x2, a2, res2,
     &                         e2, 1.0/rho2, 0.0d0)
c
         rho0 = q(0,j,1)
         u0   = q(0,j,2)/rho0
         v0   = q(0,j,3)/rho0
         e0   = q(0,j,4)/rho0 - 0.5*(u0*u0+v0*v0) - e_const
         CALL CO2BLLT_EQUI(p0, T0, c0, x0, a0, res0,
     &                         e0, 1.0/rho0, 0.0d0)
c
         deltax = aux(1,j,6)*dx
c         dpx0 = (p1-p_bc)/deltax
c         dux0 = (u1-u0)/deltax
c         drx0 = (rho1-rho0)/deltax
c         dvx0 = (v1-v0)/deltax
c
         dpx1 = (p2-p1)/deltax
         dux1 = (u2-u1)/deltax
         drx1 = (rho2-rho1)/deltax
c
c###########
c         print*, p_bc1,vguess_in(j),T_bc
c         CALL eos_1d(5, rhoi_bc, e_bc, resbc, Niter,
c     &               iflag, p_bc1, vguess_in(j), T_bc, out3)
c         vguess_in(j) = rhoi_bc
c         rho_bc = 1.0/rhoi_bc
c         print*, 'rhoBC', rho_bc, 'e_bc',e_bc, resbc
c         STOP
c##############
c
         dist = 10.0
         Z1   = 0.0 !(rho0 - rho_bc1) / dist*u1                   ! u
         Z2   = 0.0 !(v0-0.0) / dist * u1                        ! u
         Z5   = (u1-c1)*(-dux1+1.0/(rho1*c1)*dpx1)          ! u-c 
         Z4   = -Z5 !(p0 - p_bc1) /(rho1*dist)            ! u+c
c
c         D1   = Z1 + rho1/(c1*2.0)*(Z4+Z5)
c         D2   = 0.5*(Z4-Z5)
c         D3   = Z2
c         D5   = rho1*c1*0.5 * (Z4+Z5)
c
c         q(0,j,1) = q(0,j,1) - dt*(D1)
c         q(0,j,2) = q(0,j,2) - dt*(u1*D1+rho1*D2)
c         q(0,j,3) = q(0,j,3) - dt*(v1*D1+rho1*D3)
c         q(0,j,4) = q(0,j,4) - dt*(0.5*(u1*u1+v1*v1)*D1 +   !!!!! dp/drhoe * drhoe/dt = d5 !!!!!!!!! 
c     &                 rho1*u1*D2 + rho1*v1*D3 + D5)
c
         rho_bc = rho_bc1!q(0,j,1) - dt*( (Z4+Z5)*rho1/(2.0*c1)+Z1 )
         v_bc   = v1 !v0 - dt*Z2
c         u_bc   = 60.0
         u_bc   = u1 !q(0,j,2)/q(0,j,1) - dt*(Z4-Z5)*0.5
         p_bc   = p_bc1 !p0 - dt*(rho1*c1/2.0)*(Z4+Z5)
         e_bc   = e_bc1
 
c         CALL eos_1d(3, e_bc, out_2, resbc, Niter,
c     &               iflag, p_bc, e0,1.0/rho_bc,out3)
c
c         IF (resnorm>1e-3) THEN
c           print*, 'e_bc',e_bc,'p,eguess',p_bc,e0
c           print*, 'BC inlet NewRaphson problem',resnorm
c           STOP
c         ENDIF
c
         q(0,j,1) = rho_bc
         q(0,j,2) = u_bc*rho_bc
         q(0,j,3) = v_bc*rho_bc        
         q(0,j,4) = (e_bc + 0.5*(u_bc*u_bc+v_bc*v_bc)+e_const)*rho_bc
      ENDDO
      DO j = 1-mbc, my+mbc
 
         q(-1,j,1) = q(0,j,1)
         q(-1,j,2) = q(0,j,2)
         q(-1,j,3) = q(0,j,3)
         q(-1,j,4) = q(0,j,4)

      ENDDO
c
c
      go to 199
c
  110 continue
c      print*,' # zero-order extrapolation: left'
      do 115 m=1,meqn
         do 115 ibc=1,mbc
            do 115 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(1,j,m)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
            do 125 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
            do 135 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
  135       continue
c     # negate the normal velocity:
      do 136 ibc=1,mbc
         do 136 j = 1-mbc, my+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
      flag_diagno = 1
      p_bc = 3.0d6
c max Mach
         rhomax = q(mx,my/2,1)
         umax   = q(mx,my/2,2) / rhomax
         vmax   = q(mx,my/2,3) / rhomax
         emax   = q(mx,my/2,4) / rhomax -
     &            0.5*(umax*umax+vmax*vmax) - e_const
c
         CALL CO2BLLT_EQUI(pmax, Tmax, cmax, xmax, amax, resmax,
     &                     emax, 1.0d0/rhomax, 0.0d0)
        u1c1    = sqrt(umax*umax+vmax*vmax) / cmax

        IF (u1c1>1.0) THEN
           print*, 'supersonic outlet'
           go to 210
        ENDIF

      DO j=1-mbc,my+mbc

         rho1 = q(mx,j,1)
         u1   = q(mx,j,2)/rho1
         v1   = q(mx,j,3)/rho1
         e1   = q(mx,j,4)/rho1 - 0.5*(u1*u1+v1*v1) - e_const
         CALL CO2BLLT_EQUI(p1, T1, c1, x1, a1, res1,
     &                     e1, 1.0d0/rho1, 0.0d0)
c
         rho2 = q(mx-1,j,1)
         u2   = q(mx-1,j,2)/rho2
         v2   = q(mx-1,j,3)/rho2
         e2   = q(mx-1,j,4)/rho2 - 0.5*(u2*u2+v2*v2) - e_const
         CALL CO2BLLT_EQUI(p2, T2, c2, x2, a2, res2,
     &                     e2, 1.0d0/rho2, 0.0d0)
c
         rho0 = q(mx+1,j,1)
         u0   = q(mx+1,j,2)/rho0
         v0   = q(mx+1,j,3)/rho0
         e0   = q(mx+1,j,4)/rho0 - 0.5*(u0*u0+v0*v0) - e_const
         CALL CO2BLLT_EQUI(p0, T0, c0, x0, a0, res0,
     &                     e0, 1.0d0/rho0, 0.0d0)
c
         deltax = aux(mx,j,6)*dx
c         dpx0 = (p_bc-p1)/deltax
c         dux0 = (u0-u1)/deltax
c         drx0 = (rho0-rho1)/deltax
c
         dpx1 = (p1-p2)/deltax
         dux1 = (u1-u2)/deltax
         drx1 = (rho1-rho2)/deltax
         dvx1 = (v1-v2)/deltax
cc max Mach
c         rhomax = q(mx,my/2,1)
c         umax   = q(mx,my/2,2) / rhomax
c         vmax   = q(mx,my/2,3) / rhomax
c         emax   = q(mx,my/2,4) / rhomax - 
c     &            0.5*(umax*umax+vmax*vmax) - e_const
c
c         CALL CO2BLLT_EQUI(pmax, Tmax, cmax, xmax, amax, resmax,
c     &                     emax, 1.0d0/rhomax, 0.0d0)
c        u1c1    = sqrt(umax*umax+vmax*vmax) / cmax
c
c
         dist = 500.0
         Z1   = u1*(drx1 - 1.0/(c1*c1)*dpx1)
         Z2   = u1*dvx1                                       ! u
         Z4   = (u1+c1)*(dux1 + 1.0/(rho1*c1)*dpx1)                  ! u+c
c         Z5   = 0.25*c0*(1.0-(u0/c0)**2.0)/dist*(p0-p_bc)       ! u-c 
         Z5   = 0.25*c0*(1.0-(u1c1)**2.0)/dist*(p0-p_bc)       ! u-c 
c
c         rho_bc = q(mx+1,j,1) - dt*(Z1 + (Z4+Z5)*rho1*0.5/c1)
c
c         print*, 'before eos1 right',rho_bc
c---------------------------------------------------------------------------
c         CALL eos_1d(3, e_bc, out_2, resbc, Niter,
c     &               iflag, p_bc, eguess_out(j),1.0/rho_bc,out3)
c         eguess_out(j) = e_bc
c         IF (resnorm>1e-3) THEN
c           print*, 'e_bc',e_bc,'p,eguess',p_bc,eguess_out(j)
c           print*, 'BC outlet NewRaphson problem',resnorm
c           STOP
c         ENDIF
c---------------------------------------------------------------------------
         D1     =  Z1 + rho1/(c1*2.0)*(Z4+Z5)
         D2     =  0.5*(Z4-Z5)
         D3     =  Z2
         D5     =  rho1*c1*0.5 * (Z4+Z5)
c
c         v_bc   =  v0 - dt*Z2
c         u_bc   = u1
c         u_bc   = q(mx+1,j,2)/q(mx+1,j,1) - dt*(Z4-Z5)*0.5 
c         q(mx+1,j,1) = rho_bc
c         q(mx+1,j,2) = u_bc*rho_bc
c         q(mx+1,j,3) = v_bc*rho_bc
c         q(mx+1,j,4) = (e_bc + 0.5*(u_bc*u_bc+v_bc*v_bc)+e_const)*rho_bc
c
         CALL CO2DER(dp_dv_u, dp_du_v , e0, 1d0/rho0, T0, p0,
     &               res2,res3,res4)
         dvdr    =  -1d0 / (rho0*rho0)
         dedp_r  =   1d0 / dp_du_v
         dedr_p  =  -dp_dv_u*dvdr / dp_du_v
c         
         q(mx+1,j,1) = q(mx+1,j,1) - dt*(D1)
         q(mx+1,j,2) = q(mx+1,j,2) - dt*(u0*D1+rho0*D2) 
         q(mx+1,j,3) = q(mx+1,j,3) - dt*(v0*D1+rho0*D3) 
         q(mx+1,j,4) = q(mx+1,j,4) - dt*( q(mx+1,j,4)/rho0*D1 +
     &                 rho0*u0*D2 + rho0*v0*D3 +
     &                 rho0*(dedp_r*D5+dedr_p*D1) ) 
      ENDDO
c
      DO j = 1-mbc, my+mbc

         q(mx+2,j,1) = q(mx+1,j,1)
         q(mx+2,j,2) = q(mx+1,j,2)
         q(mx+2,j,3) = q(mx+1,j,3)
         q(mx+2,j,4) = q(mx+1,j,4)

      ENDDO
c
c
      go to 299

  210 continue
c      print*,' # zero-order extrapolation: right'
      do 215 m=1,meqn
         do 215 ibc=1,mbc
            do 215 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
            do 225 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
            do 235 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
  235       continue
c     # negate the normal velocity:
      do 236 ibc=1,mbc
         do 236 j = 1-mbc, my+mbc
            q(mx+ibc,j,2) = -q(mx+1-ibc,j,2)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # noslip
        do jbc=1,mbc
         do  i = 1-mbc, mx+mbc
          q(i,1-jbc,1) = q(i,jbc,1)
          q(i,1-jbc,4) = q(i,jbc,4)
         enddo
        enddo
c projection of u and v components in computational domain
         do jbc=1,mbc
          do i = 1-mbc, mx+1
c    
c            vnx = aux(i,jbc,4)
c            vny = aux(i,jbc,5)
           q(i,1-jbc,2) = -q(i,jbc,2)
           q(i,1-jbc,3) = -q(i,jbc,3)
          enddo
         enddo
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 m=1,meqn
         do 315 jbc=1,mbc
            do 315 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
  315       continue
      go to 399

  320 continue
c     # periodic:  
      do 325 m=1,meqn
         do 325 jbc=1,mbc
            do 325 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
  325       continue
      go to 399

  330 continue

c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      IF (maux == 0) THEN
      do 335 m=1,meqn
         do 335 jbc=1,mbc
            do 335 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
  335       continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
  336    continue
c
c     # Nonuniform grid
c
      ELSEIF (maux==7) THEN
c      print*,'BC nonuniform bottom'
c
        do jbc=1,mbc
         do  i = 1-mbc, mx+mbc
          q(i,1-jbc,1) = q(i,jbc,1)      
          q(i,1-jbc,4) = q(i,jbc,4)
         enddo
        enddo
c projection of u and v components in computational domain
         do jbc=1,mbc
          do i = 1-mbc, mx+1
c    
            vnx = aux(i,jbc,4)
            vny = aux(i,jbc,5)
           q(i,1-jbc,2) = q(i,jbc,2)*(-vny*vny+vnx*vnx) +
     &                    q(i,jbc,3)*2.0*vnx*vny
           q(i,1-jbc,3) = q(i,jbc,2)*2.0*vnx*vny +
     &                     q(i,jbc,3)*(vny*vny-vnx*vnx)
          enddo
         enddo
      ENDIF

      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     noslip
        do jbc=1,mbc
         do  i = 1-mbc, mx+mbc
          q(i,my+jbc,1) = q(i,my-jbc+1,1)
          q(i,my+jbc,4) = q(i,my-jbc+1,4)
         enddo
        enddo
c projection of u and v components in computational domain
c 
         do jbc=1,mbc
          do i = 1-mbc, mx+mbc
c           vnx = aux(i,my-jbc+2,4)
c            vny = aux(i,my-jbc+2,5)
           q(i,my+jbc,2) = -q(i,my-jbc+1,2)
           q(i,my+jbc,3) = -q(i,my-jbc+1,3)
c
          enddo
         enddo
c
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 m=1,meqn
         do 415 jbc=1,mbc
            do 415 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
  415       continue
      go to 499

  420 continue
c     # periodic:  
      do 425 m=1,meqn
         do 425 jbc=1,mbc
            do 425 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      IF (maux == 0) THEN
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
  435       continue
c     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
  436    continue
c
c       # IF nonuniform grid
c
      ELSEIF (maux==7) THEN
c      print*,'BC nonuniform top'
       do jbc=1,mbc
         do  i = 1-mbc, mx+mbc
          q(i,my+jbc,1) = q(i,my-jbc+1,1)
          q(i,my+jbc,4) = q(i,my-jbc+1,4)
         enddo
        enddo
c projection of u and v components in computational domain
c 
         do jbc=1,mbc
          do i = 1-mbc, mx+mbc
            vnx = aux(i,my-jbc+2,4)
            vny = aux(i,my-jbc+2,5)
           q(i,my+jbc,2) = q(i,my-jbc+1,2)*(-vny*vny+vnx*vnx)+
     &                     q(i,my-jbc+1,3)*2.0*vnx*vny
           q(i,my+jbc,3) = q(i,my-jbc+1,2)*2.0*vnx*vny+
     &                      q(i,my-jbc+1,3)*(vny*vny-vnx*vnx)
c
          enddo
         enddo
      ENDIF

      go to 499

  499 continue
c         print*,'bc finish for t=', t
      return
      end SUBROUTINE bc2
