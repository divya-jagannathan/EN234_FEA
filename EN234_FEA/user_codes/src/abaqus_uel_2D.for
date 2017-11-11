!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      
      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  Nbar(9)
      double precision  ::  dNbardxi(9,2)
      double precision  ::  dNbardx(9,2)
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords

      double precision  ::  sol(9), dsol(9)                   ! Sol vector contains [mu, c, dmudx1, dmudx2, dcdx1, dcdx2]
      double precision  ::  q(9)                              ! q vector defined in class
      double precision  ::  D(9,9), Dold(4,4)                 ! D matrix defined in class
      double precision  ::  B(9,24)             ! p = B*U
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      double precision  ::  diffusion_coeft,kappa,theta       ! Material properties
      double precision  ::  E,xnu,Omega,xw                     ! Material properties
      double precision  ::  c,xmu, dc, dxmu,stress(4),s33, enew(4)    ! concentration

    !
    !     Example ABAQUS UEL implementing 2D phase field model

      if (NNODE == 3) n_points = 4
      if (NNODE == 4) n_points = 4
      if (NNODE == 6) n_points = 4
      if (NNODE == 8) n_points = 4
      if (NNODE == 9) n_points = 9

      call abq_UEL_2D_integrationpoints(4, NNODE, xi, w)

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      E = PROPS(1)
      xnu = PROPS(2)
      Omega = PROPS(3)
      xw=PROPS(4)
      kappa=PROPS(5)
      diffusion_coeft=PROPS(6)
      theta=PROPS(7)

      D = 0.d0
      d44 = 0.5D0*E/(1+xnu)
      d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
      d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
      d13=d12
      d23=d12
      d21=d12
      d22=d11
      D(1,1) = d11
      D(2,2) = d11
      D(3,3) = d44
      D(1,2) = d12
      D(2,1) = D(1,2)
      D(1,5)=-(Omega/3.d0)*(d11+d12+d13)
      D(2,5)=-(Omega/3.d0)*(d21+d22+d23)
      D(4,1)=D(1,5)
      D(4,2)=D(2,5)
      D(4,4)=1.d0
    ! D(4,5) defined inside the integration loop 
      D(5,5)= 1/DTIME
      D(6,8)= -kappa
      D(7,9)= -kappa
      D(8,6)= theta*diffusion_coeft
      D(9,7)= theta*diffusion_coeft
  
    !     --  Loop over integration points
      do kint = 1, 4
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),8,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        call abq_inverse_LU(dxdxi,dxidx,2)
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dNdx(1:8,1:2) = matmul(dNdxi(1:8,1:2),dxidx)
        
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,Nbar,dNbardxi)
      !  dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:4,1:2))
       ! call abq_inverse_LU(dxdxi,dxidx,2)
        !determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dNbardx(1:4,1:2) = matmul(dNbardxi(1:4,1:2),dxidx)
        
        B = 0.d0
        B(1,1:16:4) = dNdx(1:4,1)
        B(1,17:24:2)= dNdx(5:8,1)
        B(2,2:16:4) = dNdx(1:4,2)
        B(2,18:24:2)= dNdx(5:8,2)
        B(3,1:16:4) = dNdx(1:4,2)
        B(3,17:24:2)= dNdx(5:8,2)
        B(3,2:16:4) = dNdx(1:4,1)
        B(3,18:24:2)= dNdx(5:8,1)
        B(4,3:16:4) = Nbar(1:4)
        B(5,4:16:4) = Nbar(1:4)
        B(6,3:16:4) = dNbardx(1:4,1)
        B(7,3:16:4) = dNbardx(1:4,2)
        B(8,4:16:4) = dNbardx(1:4,1)
        B(9,4:16:4) = dNbardx(1:4,2)
        
        sol = matmul(B(1:9,1:24),U(1:24))      ! The p vector at the end of the step
        dsol = matmul(B(1:9,1:24),DU(1:24,1))    ! Increment in the p vector

        c = sol(5)
        xmu = sol(4)
       
        
        enew(1)=sol(1)-c*Omega/3
        enew(2)=sol(2)-c*Omega/3
        enew(3)=-c*Omega/3
        enew(4)=sol(3)
        
        Dold=0.d0
        Dold(1,1)=d11
        Dold(2,2)=d11
        Dold(3,3)=d11
        Dold(4,4)=d44
        Dold(1,2)=d12
        Dold(1,3)=d12
        Dold(2,1)=d12
        Dold(2,3)=d12
        Dold(3,1)=d12
        Dold(3,2)=d12
        
        stress(1:4)=matmul(Dold(1:4,1:4),enew(1:4))
        
        q(1) = stress(1)
        q(2) = stress(2)
        q(3) = stress(4)
        q(4) = xmu-2*xw*(c)*(c-1)*(2*(c)-1)-
     1   (Omega/3)*(stress(1)+stress(2)+stress(3))
        q(5) = dsol(5)/DTIME
        q(6) = -kappa*(sol(8))
        q(7) = -kappa*(sol(9))
        q(8) = diffusion_coeft*(sol(6)+(theta-1)*dsol(6))
        q(9) = diffusion_coeft*(sol(7)+(theta-1)*dsol(7))

        D(4,5) = -xw*(12*(c)**2-12*(c)+2) +
     1   ((Omega**2)/9)*(3*d11+6*d12)


        RHS(1:24,1) = RHS(1:24,1) 
     1   - matmul(transpose(B(1:9,1:24)),q)*w(kint)*determinant

        AMATRX(1:24,1:24) = AMATRX(1:2*NNODE,1:2*NNODE) 
     1   + matmul(transpose(B(1:9,1:24)),
     2             matmul(D,B(1:9,1:24)))*w(kint)*determinant
        
        
         if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
             SVARS(4*kint-3:4*kint) = stress(1:4)
        endif
  

      end do
  


      return

      END SUBROUTINE UEL


         

      subroutine abq_UEL_2D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_UEL_2D_integrationpoints



      subroutine abq_UEL_2D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2

            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if

      end subroutine abq_UEL_2D_shapefunctions


      subroutine abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)


      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)



      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in subroutine abq_UEL_1D_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select







      end subroutine ABQ_UEL_1D_integrationpoints



      subroutine abq_facenodes_2D(nelnodes,face,list,nfacenodes)

      implicit none

      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes
    !
    !        Subroutine to return list of nodes on an element face for standard 2D solid elements
    !
      integer :: i3(3)
      integer :: i4(4)

      i3(1:3) = [2,3,1]
      i4(1:4) = [2,3,4,1]

      if (nelnodes == 3) then
        nfacenodes = 2
        list(1) = face
        list(2) = i3(face)
      else if (nelnodes == 4) then
        nfacenodes = 2
        list(1) = face
        list(2) = i4(face)
      else if (nelnodes == 6) then
        nfacenodes = 3
        list(1) = face
        list(2) = i3(face)
        list(3) = face+3
      else if (nelnodes == 8) then
        nfacenodes = 3
        list(1) = face
        list(2) = i4(face)
        list(3) = face+4
      else if (nelnodes == 9) then
        nfacenodes = 3
        if (face==1) list(1:3) = (/1,3,2/)
        if (face==2) list(1:3) = (/3,9,6/)
        if (face==3) list(1:3) = (/9,7,8/)
        if (face==4) list(1:3) = (/7,1,4/)
      endif

      end subroutine abq_facenodes_2d

      subroutine abq_inverse_LU(Ain,A_inverse,n)  ! Compute the inverse of an arbitrary matrix by LU decomposition

        implicit none

        integer, intent(in)  :: n

        double precision, intent(in)    :: Ain(n,n)
        double precision, intent(out)   :: A_inverse(n,n)

        double precision :: A(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
        double precision :: coeff
        integer :: i, j, k

        A(1:n,1:n) = Ain(1:n,1:n)
        L=0.d0
        U=0.d0
        b=0.d0

        do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
            end do
        end do

        forall (i=1:n)  L(i,i) = 1.d0
        forall (j=1:n) U(1:j,j) = A(1:j,j)
        do k=1,n
            b(k)=1.d0
            d(1) = b(1)
            do i=2,n
                d(i)=b(i)
                d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
            end do
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                x(i)=x(i)-dot_product(U(i,i+1:n),x(i+1:n))
                x(i) = x(i)/U(i,i)
            end do
            A_inverse(1:n,k) = x(1:n)
            b(k)=0.d0
        end do

      end subroutine abq_inverse_LU


      subroutine abq_UEL_invert2d(A,A_inverse,determinant)

      double precision, intent(in) :: A(2,2)
      double precision, intent(out) :: A_inverse(2,2)
      double precision, intent(out) :: determinant

      double precision COFACTOR(2,2)

!   Compute inverse and determinant of 2x2 matrix

      determinant =   A(1,1)*A(2,2)
     1   - A(1,2)*A(2,1)
     

      IF (determinant==0.d0) THEN
        write(6,*) ' Error in subroutine abq_UEL_inver3d'
        write(6,*) ' A 3x3 matrix has a zero determinant'
        stop
      endif
      COFACTOR(1,1) = +A(2,2)
      COFACTOR(1,2) = -A(2,1)
      COFACTOR(2,1) = -A(1,2)
      COFACTOR(2,2) = +A(1,1)
 
      A_inverse = transpose(COFACTOR) / determinant


      end subroutine abq_UEL_invert2d
      
      subroutine neohooke(element_properties,n_properties,G,F,J,qvec,
     1 Sigma,H,D)

       implicit none

       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: element_properties(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: J
       double precision, intent(in)  :: G(6,6)
       double precision, intent(out) :: Sigma(6)
       double precision, intent(out) :: D(6,6)
       double precision, intent(out) :: H(6,9)
       double precision, intent(out) :: qvec(9)

       double precision :: C(3,3)
       double precision :: Cvec(6)
       double precision :: Cinv(3,3)
       double precision :: Cinvvec(6)
       double precision :: Cbarvec(6)
       double precision :: Cstarvec(6)
       double precision :: Cbarstarvec(6)
       double precision :: eyevec(6)
       double precision :: mu
       double precision :: K
       double precision :: ss
       double precision :: Qval, q(3,3)
       double precision :: Pvec(6)
       double precision :: Omega(6,6), D1(6,6), D2(6,6), D3(6,6)
       double precision :: P1(6,6),P2(6,6), C1(6,6), Sigmamat(3,3)
       integer :: i

       !  This subroutine calculates the Kirchhoff stress tau = J*cauchy_stress (stored as a vector stress(i) = [tau_11, tau_22, tau_33, etc]
       !  and the tangent matrix D[I,J] = [dtau_11/dB_11, dtau_11/dB_22,... 
       !                                   dtau_22/dB_11, dtau_22/dB_22,
       !                                    etc
       
       mu = element_properties(1)
       K  = element_properties(2)

       C = matmul(transpose(F),F)
       
       Cvec(1) = C(1,1)
       Cvec(2) = C(2,2)
       Cvec(3) = C(3,3)
       Cvec(4) = C(1,2)
       Cvec(5) = C(1,3)
       Cvec(6) = C(2,3)
       
       Cstarvec(1:6)=Cvec(1:6)
       Cstarvec(4:6)=2.d0*Cvec(4:6)
       
       call abq_UEL_invert3d(C,Cinv,ss)      ! ss is just a dummy variable here

       Cinvvec(1) = Cinv(1,1)
       Cinvvec(2) = Cinv(2,2)
       Cinvvec(3) = Cinv(3,3)
       Cinvvec(4) = Cinv(1,2)
       Cinvvec(5) = Cinv(1,3)
       Cinvvec(6) = Cinv(2,3)
       
       ss = J**(-2.d0/3.d0)
       
       do i = 1,6
         Cbarvec(i) = Cvec(i)*ss
         Cbarstarvec(i)=Cstarvec(i)*ss
       end do
       
       eyevec(1:3) = 1.d0
       eyevec(4:6) = 0.d0
      
       Qval=(1/4)*dot_product(Cbarvec-
     1  eyevec,matmul(G,Cbarvec(1:6)-eyevec(1:6)))

       Pvec=0.5*ss*(matmul(G,(Cbarstarvec-
     1    eyevec))-(1/3)*dot_product(Cstarvec,
     2    matmul(G,(Cbarstarvec-eyevec)))*Cinvvec)
       
       Sigma=mu*exp(Qval)*Pvec+K*J*(J-1)*Cinvvec
       
        Sigmamat(1:3,1:3)=0.d0
        Sigmamat(1,1)=Sigma(1)
        Sigmamat(2,2)=Sigma(2)
        Sigmamat(3,3)=Sigma(3)
        Sigmamat(1,2)=Sigma(4)
        Sigmamat(2,1)=Sigma(4)
        Sigmamat(1,3)=Sigma(5)
        Sigmamat(3,1)=Sigma(5)
        Sigmamat(2,3)=Sigma(6)
        Sigmamat(3,2)=Sigma(6)
       
       q=0.d0
      
       q=matmul(Sigmamat,transpose(F))
     
       qvec(1)=q(1,1)
       qvec(2)=q(2,2)
       qvec(3)=q(3,3)
       qvec(4)=q(2,1)
       qvec(5)=q(1,2)
       qvec(6)=q(3,1)
       qvec(7)=q(1,3)
       qvec(8)=q(3,2)
       qvec(9)=q(2,3)
       
       
       H(1:6,1:9)=0.d0
       
       H(1,1)=F(1,1)
       H(1,5)=F(2,1)
       H(1,7)=F(3,1)
       H(2,2)=F(2,2)
       H(2,4)=F(1,2)
       H(2,9)=F(3,2)
       H(3,3)=F(3,3)
       H(3,6)=F(1,3)
       H(3,8)=F(2,3)
       H(4,1)=F(1,2)
       H(4,2)=F(2,1)
       H(4,4)=F(1,1)
       H(4,5)=F(2,2)
       H(4,7)=F(3,2)
       H(4,9)=F(3,1)
       H(5,1)=F(1,3)
       H(5,3)=F(3,1)
       H(5,5)=F(2,3)
       H(5,6)=F(1,1)
       H(5,7)=F(3,3)
       H(5,8)=F(2,1)
       H(6,2)=F(2,3)
       H(6,3)=F(3,2)
       H(6,4)=F(1,3)
       H(6,6)=F(1,2)
       H(6,8)=F(2,2)
       H(6,9)=F(3,3)
       
      Omega(1:6,1:6)=0.d0
      
      do i=1,3
          Omega(i,i:3)=Cinv(i,i:3)*Cinv(i,i:3)  
      end do
      
      Omega(1:3,4)=Cinv(1:3,1)*Cinv(1:3,2)
      Omega(1:3,5)=Cinv(1:3,1)*Cinv(1:3,3)
      Omega(1:3,6)=Cinv(1:3,2)*Cinv(1:3,3)
     
      Omega(4,4)=(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(1,2))/2
      Omega(4,5)=(Cinv(1,1)*Cinv(2,3)+Cinv(1,3)*Cinv(1,2))/2
      Omega(4,6)=(Cinv(1,2)*Cinv(2,3)+Cinv(1,3)*Cinv(2,2))/2
      Omega(5,5)=(Cinv(1,1)*Cinv(3,3)+Cinv(1,3)*Cinv(1,3))/2
      Omega(5,6)=(Cinv(1,2)*Cinv(3,3)+Cinv(1,3)*Cinv(2,3))/2
      Omega(6,6)=(Cinv(2,2)*Cinv(3,3)+Cinv(2,3)*Cinv(2,3))/2
      
      Omega(1:6,1:6)=-Omega(1:6,1:6)-transpose(Omega(1:6,1:6))
      
      do i=1,6
          Omega(i,i)=Omega(i,i)/2
      end do
      
      
      
      D1(1:6,1:6)=0.d0
      D2(1:6,1:6)=0.d0
      D3(1:6,1:6)=0.d0
      P1(1:6,1:6)=0.d0
      P2(1:6,1:6)=0.d0
      C1(1:6,1:6)=0.d0
      
      D1=matmul(G,spread(Cstarvec,dim=2,
     1 ncopies=6)*spread(Cinvvec,dim=2,ncopies=6))+
     2 spread(Cinvvec,dim=1, ncopies=6)*spread(matmul(
     3 G,Cstarvec),dim=1,ncopies=6) 
      
      D2=dot_product(Cstarvec(1:6),matmul(G(1:6,1:6),
     1 (Cbarstarvec(1:6)-eyevec(1:6))))*Omega(1:6,1:6)
      
      D3=dot_product(Cstarvec,matmul(G,Cstarvec))*spread(Cinvvec,
     1 dim=2,ncopies=6)*spread(Cinvvec,dim=1,ncopies=6) 
      
      P1=spread(Pvec,dim=2,ncopies=6)*spread((Pvec-(1/3)*
     1 Cinvvec),dim=1,ncopies=6)
      
      P2=spread(Cinvvec,dim=2,ncopies=6)*spread(matmul(G,(Cbarstarvec-
     1  eyevec)),dim=1,ncopies=6) 
      
      C1=spread(Cinvvec,dim=2,ncopies=6)*
     1  spread(Cinvvec,dim=1,ncopies=6)
      
      
      D(1:6,1:6)=0.d0
      
      D(1:6,1:6)=mu*exp(Qval)*(ss**2)*(G(1:6,1:6)-(1/3)*D1(1:6,1:6)-
     1 ((ss**(-1))/3)*D2(1:6,1:6)+(1/9)*D3(1:6,1:6)) +
     2 mu*exp(Qval)*(2*P1(1:6,1:6)-(ss/3)*P2(1:6,1:6)) + 
     3  K*J*((2*J-1)*C1(1:6,1:6)+2*(J-1)*Omega(1:,1:6)) 
      
       return

      end subroutine neohooke


