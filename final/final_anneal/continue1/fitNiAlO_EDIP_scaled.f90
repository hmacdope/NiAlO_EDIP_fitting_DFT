module fit
   implicit none
   integer, parameter::np = 51
   integer, parameter::MaxAt = 600
   integer, parameter::MaxConf = 110

   logical, parameter::NoGas = .false., RefN2 = .false., TwoBodyOnly = .false., NoDissoc = .false.
   real(kind=8), dimension(3, MaxAt, MaxConf)::x
   real(kind=8), dimension(3, 3, MaxConf)::Q
   real(kind=8), dimension(MaxConf)::En, model, twobody, threebody, target, EAM, VN, ScalFac
   character(len=20), dimension(MaxConf)::confnames
   logical, parameter::debug = .false.
   integer, dimension(MaxConf)::NAt
   character(len=2), dimension(MaxAt, MaxConf)::lab
   integer::NConf

contains
   ! Read in the configurations
   subroutine LoadConfigs
      integer::i, k, l, badness
      character(len=8)::str
      character::sk
      NConf = 0
      NAt = 0
      open (unit=10, file='master_configs.xyz')
      do NConf = NConf + 1, MaxConf
         read (10, *, iostat=badness) l
         if (badness /= 0) exit
         read (10, *) str, En(NConf)
         read (10, *)
         read (10, *) str, Q(:, 1, NConf)
         read (10, *) str, Q(:, 2, NConf)
         read (10, *) str, Q(:, 3, NConf)
         read (10, *)
         read (10, *)
         read (10, *)
         read (10, *)
         NAt(NConf) = l - 8
         do i = 1, NAt(NConf)
            read (10, *) lab(i, NConf), x(:, i, NConf)
         end do
         print *,  NConf
      end do
      NConf = NConf - 1
      close (10)
      
      ! How many confs?
      print *, "NCONF FOUND"
      print *, NConf
      print *, ""

      ! Read in traget energies
      open (unit=10, file='V_energies.csv')
      read (10, *)
      do i = 1, NConf 
         read (10, *) confnames(i), target(i)
      end do
      close (10)

      print *, " NAME, DFT ENERGY, TARGET ENERGY"
      do i = 1, NConf
         print *, confnames(i),  En(i),  target(i) 
      end do
      print *, ""

      print *, "STARTING FIT ROUTINE"
      print *, ""
      print *, ""

      
   end subroutine LoadConfigs

   subroutine Energies(P2b, PZ, P3b)
      real(kind=8), dimension(5, 3), intent(in)::P2b  ! A:1, B:2, rho:3, beta:4, sigma:5
      real(kind=8), dimension(3, 3), intent(in)::PZ   ! alpha:1, cutoffA:2, cutoffCfrac:3
      real(kind=8), dimension(9, 3), intent(in)::P3b  ! gamma:1, lambda:2, eta:3, Q0:4, mu:5, u1-u4:6-9
      real(kind=8), dimension(3)::rv1, rv2
      real(kind=8)::r, r1, r2, cth, tau, Qz, h
      real(kind=8)::Z, f
      integer::i, l, m, n, i1, i2, i3
      integer::l1, l2, m1, m2, n1, n2
      integer::ty, ty2, lim
      character(len=10)::dts
      model = 0
      twobody = 0
      threebody=0
      if (debug) print *, "Computing Energies", ""
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(ty,ty2,lim,r,r1,r2,rv1,rv2,cth,dts,tau,Qz,Z,f,h)
      do i = 1, NConf
      if (debug) print *, "CONF", i
         ty = 0
         lim = 1
         do i1 = 1, NAt(i)
            if (lab(i1, i) /= 'O') cycle
            Z = 0
            do i2 = 1, NAt(i)
               if (lab(i2, i) == 'Ni') ty = 1
               if (lab(i2, i) == 'Al') ty = 2
               if (lab(i2, i) == 'O') ty = 3 
               do l = -lim, lim
                  do m = -lim, lim
                     do n = -lim, lim
                        if ((i2 == i1) .and. (all((/l, m, n/) == 0))) cycle
                        r = sum((x(:, i1, i) - x(:, i2, i) + matmul(Q(:, :, i), real((/l, m, n/), 8)))**2)
                        if (r > PZ(2, ty)**2) cycle
                        r = sqrt(r)
                        if (r <= PZ(3, ty)*PZ(2, ty)) then
                           Z = Z + 1
                        else
                           f = (r - PZ(3, ty)*PZ(2, ty))/(PZ(2, ty) - PZ(3, ty)*PZ(2, ty))
                           Z = Z + exp(PZ(1, ty)/(1 - (1/f)**3))
                        end if
                     end do
                  end do
               end do
            end do
            ! two-body
            do i2 = 1, NAt(i)
               if (lab(i2, i) == 'Ni') ty = 1
               if (lab(i2, i) == 'Al') ty = 2
               if (lab(i2, i) == 'O') ty = 3 
               do l = -lim, lim
                  do m = -lim, lim
                     do n = -lim, lim
                        if ((i2 == i1) .and. (all((/l, m, n/) == 0))) cycle
                        r = sum((x(:, i1, i) - x(:, i2, i) + matmul(Q(:, :, i), real((/l, m, n/), 8)))**2)
                        if (r > PZ(2, ty)**2) cycle
                        r = sqrt(r)
                        model(i) = model(i) + P2b(1, ty)*((P2b(2, ty)/r)**P2b(3, ty) - exp(-P2b(4, ty)*Z**2)) &
                                   *exp(P2b(5, ty)/(r - PZ(2, ty)))
                        ! twobody(i) = twobody(i) + P2b(1, ty)*((P2b(2, ty)/r)**P2b(3, ty) - exp(-P2b(4, ty)*Z**2)) &
                        !            *exp(P2b(5, ty)/(r - PZ(2, ty)))
                     end do
                  end do
               end do
            end do
            if (TwoBodyOnly) cycle
            ! three body
            do i2 = 1, NAt(i)
               if (i1 == i2) cycle
               if (lab(i2, i) == 'Ni') ty = 1
               if (lab(i2, i) == 'Al') ty = 2
               if (lab(i2, i) == 'O') ty = 3
               do l1 = -lim, lim
                  do m1 = -lim, lim
                     do n1 = -lim, lim
                        r = sum((x(:, i1, i) - x(:, i2, i) - matmul(Q(:, :, i), real((/l1, m1, n1/), 8)))**2)
                        if (r > PZ(2, ty)**2) cycle ! FIXED
                        do i3 = i2 + 1, NAt(i)
                           if (i1 == i3) cycle
                           if (lab(i2, i) /= lab(i3, i)) cycle
                           do l2 = -lim, lim
                              do m2 = -lim, lim
                                 do n2 = -lim, lim
                                    r = sum((x(:, i1, i) - x(:, i3, i) - matmul(Q(:, :, i), real((/l2, m2, n2/), 8)))**2)
                                    if (r > PZ(2, ty)**2) cycle ! FIXED
                                    rv1 = x(:, i2, i) + matmul(Q(:, :, i), real((/l1, m1, n1/), 8)) - x(:, i1, i)
                                    rv2 = x(:, i3, i) + matmul(Q(:, :, i), real((/l2, m2, n2/), 8)) - x(:, i1, i)
                                    ! we have X-N-X vectors
                                    r1 = sqrt(sum(rv1**2))
                                    r2 = sqrt(sum(rv2**2))
                                    cth = sum(rv1*rv2)/r1/r2
                                    Qz = P3b(4, ty)*exp(-P3b(5, ty)*Z)
                                    tau = P3b(6, ty) + P3b(7, ty)*(P3b(8, ty)*exp(-P3b(9, ty)*Z) - exp(-2*P3b(9, ty)*Z))
                                    h = P3b(2, ty)*((1 - exp(-Qz*(cth + tau)**2)) + P3b(3, ty)*Qz*(cth + tau)**2)
                                    model(i) = model(i) + exp(P3b(1, ty)/(r1 - PZ(2, ty)))*exp(P3b(1, ty)/(r2 - PZ(2, ty)))*h
                                    ! threebody(i) = threebody(i)  +&
                                    !         exp(P3b(1, ty)/(r1 - PZ(2, ty)))*exp(P3b(1, ty)/(r2 - PZ(2, ty)))*h
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         if (debug) then
           print *, confnames(i), model(i), twobody(i), threebody(i) 
         end if
      end do
!$OMP END PARALLEL DO
   end subroutine Energies

   function Cost(pa) result(R)
      real(kind=8), dimension(:), intent(in)::pa
      real(kind=8)::R
      real(kind=8), dimension(5, 3)::P2b  ! A:1, B:2, rho:3, beta:4, sigma:5
      real(kind=8), dimension(3, 3)::PZ   ! alpha:1, cutoffA:2, cutoffCfrac:3
      real(kind=8), dimension(9, 3)::P3b  ! gamma:1, lambda:2, eta:3, Q0:4, mu:5, u1-u4:6-9
      real(kind=8), dimension(MaxConf)::dE
      P2b = reshape(pa(1:15), (/5, 3/))
      PZ = reshape(pa(16:24), (/3, 3/))
      if (.not. TwoBodyOnly) then
         P3b = reshape(pa(25:51), (/9, 3/))
      else
         P3b = 0
      end if
      call Energies(P2b, PZ, P3b)

      ! generate a scaling factor for DeltaE, scale by 2/num_ox for the large systems 
      ! to get approximately the same scale as the scans, otherwise one is penalised in favour of other
      
      ! set all to one 
      ScalFac(1:NConf) = 1.0
      
      ! specific scaling 

      ! Al-Slab-Olayer1
      !ScalFac(1) = 2.0/64.0 
      ! Al-Slab-Olayer4
      !ScalFac(2) = 2.0/128.0
      ! Al2O3-O
      !ScalFac(3) = 2.0/163.0
      ! Al2O3-O2
      !ScalFac(4) = 2.0/164.0
      ! Al2O3-Olayer1
      !Scalfac(5) = 2.0/226.0
      ! Al2O3-Olayer2
      !ScalFac(6) = 2.0/290.0 
      !Al2O3-trigonal
      !ScalFac(7) = 2.0/162.0
      !Ni-Slab-Olayer1
      !ScalFac(10) = 2.0/64.0
      ! Ni-Slab-Olayer2
      !ScalFac(11) = 2.0/128.0
      ! NiO-Slab
      !ScalFac(16) = 2.0/108.0
      ! NiO-supercell
      !ScalFac(17) = 2.0/108.0
      !PureAl2O3
      !ScalFac(22) = 2.0/18.0


      VN(1:NConf) = model(1:NConf)
      dE(1:Nconf) = target(1:Nconf) - VN(1:Nconf)
      dE(1:NConf) = dE(1:NConf) * ScalFac(1:NConf) 
      R = sum(dE(1:NConf)**2)
      if (debug) print *, "R2", R

   end function Cost

end module fit

SUBROUTINE FCN(ln, lx, F)
   use fit
   implicit none
   INTEGER ln
   DOUBLE PRECISION lx(ln), F
   F = Cost(lx)
END

program fitter
   use fit
   implicit none
   integer, parameter::neps = 10
   real(kind=8), dimension(np)::par, lb, ub, cc, VM, po, xp
   real(kind=8), dimension(neps)::fstar
   integer, dimension(np)::iw
   real(kind=8)::RT, etol, T, mf
   integer::ns, nt, mev, s1, s2, na2, nfe, nte, badness

   par = 0
   par(17:15 + 3*3:3) = 4



   par = (/ &
  215.00,  1.3732,       3.6041,      0.18335E-01,   7.2763,       50.895,       1.3430,     5.7048,      0.31194E-01,   4.1795, &
198.10,      0.60332, &
 0.95085E-01,  0.79011E-02,   1.2834,       67.868,       3.8360,      0.48651,       6.0779,       3.8287,      0.46536, &
26.117,       2.3828,      0.65982, &
  10.588,       46.761,       25.695,       7219.2,       1.5919,      -1.7088,       72.223,      0.45687,      0.85828,&
9.4833,       90.176,       28.008, &
  6048.6,       3.5554,       2.6398,       280.20,      0.66633,      0.12086,       1.9585,       84.518,       28.874,&
3899.7,       6.6687,       1.6829, &
  12.402,      0.27533,      0.63484/)

   lb = 0
!real(kind=8),dimension(5,3),intent(in)::P2b  ! A:1, B:2, rho:3, beta:4, sigma:5
   ub(1:15:5) = 350
   ub(2:15:5) = 9
   ub(3:15:5) = 10
   ub(4:15:5) = 5
   ub(5:15:5) = 8
!real(kind=8),dimension(3,3),intent(in)::PZ   ! alpha:1, cutoffA:2, cutoffCfrac:3
   ub(16:15 + 3*3:3) = 150
   ub(17:15 + 3*3:3) = 9
   ub(18:15 + 3*3:3) = 0.95_8
   if (.not. TwoBodyOnly) then
      ! real(kind=8),dimension(9,3)::P3b  ! gamma:1, lambda:2, eta:3, Q0:4, mu:5, u1-u4:6-9
      ub(25:51:9) = 11
      ub(26:51:9) = 100
      ub(27:51:9) = 30
      ub(28:51:9) = 10000
      ub(29:51:9) = 15
      lb(30:51:9) = -3
      ub(30:51:9) = 3
      ub(31:51:9) = 350
      ub(32:51:9) = 2
      ub(33:51:9) = 2
   end if

   par = max(par, lb + 1.0e-8)
   par = min(par, ub - 1.0e-8)


   call LoadConfigs

   mf = Cost(par)

   print *, "TARGET", "VN", "DIFF", "SCALED DIFF"
   do nt = 1, NConf
      write (*,*) confnames(nt), target(nt), VN(nt), VN(nt) - target(nt), ((VN(nt) -target(nt))*ScalFac(nt))
   enddo

   T = 5.0
   RT = 0.85
   etol = 2.0e-3
   ns = 10 
   nt = 40
   mev = huge(1)
   cc = 2
   s1 = 1011
   s2 = 1782
   VM = 1.0
   call SA(np, par, .false., RT, etol, ns, nt, neps, mev, lb, ub, cc, 1, s1, s2, T, VM, &
           po, mf, na2, nfe, nte, badness, fstar, xp, iw)
   write (*, *)
   write (*, *) badness
   write (*, *) real(po)
   par = po

   mf = Cost(par)
   print *, "TARGET", "VN", "DIFF"
   do nt = 1, NConf
      write (*,*) confnames(nt), target(nt), VN(nt), VN(nt) - target(nt)
   enddo

end program fitter