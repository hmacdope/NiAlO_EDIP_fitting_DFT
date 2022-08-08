module fit
   implicit none
   integer, parameter::np = 51
   integer, parameter::MaxAt = 300
   integer, parameter::MaxConf = 110

   logical, parameter::NoGas = .false., RefN2 = .false., TwoBodyOnly = .false., NoDissoc = .false.
   real(kind=8), dimension(3, MaxAt, MaxConf)::x
   real(kind=8), dimension(3, 3, MaxConf)::Q
   real(kind=8), dimension(MaxConf)::En, model, target, EAM, VN
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
      if (debug) print *, "Computing Energies", ""
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(ty,ty2,lim,r,r1,r2,rv1,rv2,cth,dts,tau,Qz,Z,f,h)
      do i = 1, NConf
         ty = 0
         lim = 1
!if (.not.any(i==(/31,32,33,34,35/))) cycle
         if (any(i == (/13, 14, 15, 16, 25, 30/))) cycle
!  call date_and_time(time=dts)
!  write(*,*)i,'start ',dts
!  if (NoGas.and.((i==23).or.(i==25).or.(i>=31))) cycle
         if (NoGas .and. ((i == 23) .or. (i == 25) .or. ((i >= 31) .and. (i < 83)) .or. (i >= 79))) cycle
         if (NoDissoc .and. (i >= 36) .and. (i <= 82)) cycle
         if ((i == 83) .or. (i == 84)) lim = 6
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
                        if (r > PZ(2, ty)) cycle
                        do i3 = i2 + 1, NAt(i)
!              if (lab(i2,i)=='Ni') ty2=1
!              if (lab(i2,i)=='Al') ty2=2
!              if (lab(i2,i)=='N') ty2=3
                           if (i1 == i3) cycle
                           if (lab(i2, i) /= lab(i3, i)) cycle
                           do l2 = -lim, lim
                              do m2 = -lim, lim
                                 do n2 = -lim, lim
                                    r = sum((x(:, i1, i) - x(:, i3, i) - matmul(Q(:, :, i), real((/l2, m2, n2/), 8)))**2)
                                    if (r > PZ(2, ty)) cycle
                                    rv1 = x(:, i2, i) + matmul(Q(:, :, i), real((/l1, m1, n1/), 8)) - x(:, i1, i)
                                    rv2 = x(:, i3, i) + matmul(Q(:, :, i), real((/l2, m2, n2/), 8)) - x(:, i1, i)
                                    ! we have X-N-X vectors
                                    r1 = sqrt(sum(rv1**2))
                                    r2 = sqrt(sum(rv2**2))
                                    cth = sum(rv1*rv2)/r1/r2
!real(kind=8),dimension(9,3),intent(in)::P3b  ! gamma:1, lambda:2, eta:3, Q0:4, mu:5, u1-u4:6-9
                                    Qz = P3b(4, ty)*exp(-P3b(5, ty)*Z)
                                    tau = P3b(6, ty) + P3b(7, ty)*(P3b(8, ty)*exp(-P3b(9, ty)*Z) - exp(-2*P3b(9, ty)*Z))
                                    h = P3b(2, ty)*((1 - exp(-Qz*(cth + tau)**2)) + P3b(3, ty)*Qz*(cth + tau)**2)
                                    model(i) = model(i) + exp(P3b(1, ty)/(r1 - PZ(2, ty)))*exp(P3b(1, ty)/(r2 - PZ(2, ty)))*h
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
!  call date_and_time(time=dts)
!  write(*,*)i,' done ',dts
         if (debug)  print *, confnames(i), model(i) 
      end do
!$OMP END PARALLEL DO
!write(*,*)model(31:35)
!stop
   end subroutine Energies

   function Cost(pa) result(R)
      real(kind=8), dimension(:), intent(in)::pa
      real(kind=8)::R
      real(kind=8), dimension(5, 3)::P2b  ! A:1, B:2, rho:3, beta:4, sigma:5
      real(kind=8), dimension(3, 3)::PZ   ! alpha:1, cutoffA:2, cutoffCfrac:3
      real(kind=8), dimension(9, 3)::P3b  ! gamma:1, lambda:2, eta:3, Q0:4, mu:5, u1-u4:6-9
!real(kind=8),dimension(7)::B,ctb,C,g,r0
      real(kind=8), dimension(MaxConf)::dE
      P2b = reshape(pa(1:15), (/5, 3/))
      PZ = reshape(pa(16:24), (/3, 3/))
      if (.not. TwoBodyOnly) then
         P3b = reshape(pa(25:51), (/9, 3/))
      else
         P3b = 0
      end if
      call Energies(P2b, PZ, P3b)

      VN(1:NConf) = model(1:NConf)
      dE(1:Nconf) = target(1:Nconf) - VN(1:Nconf)

      R = sum(dE(1:NConf)**2)
      if (debug) print *, "R", R

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



   par = (/6.270725, 0.9256532, 1.529872, 1.3633298E-02, 0.2077398, &
           17.27714, 1.177238, 10.00000, 3.2836724E-02, 2.173588, &
           77.14455, 0.8241709, 9.054739, 0.1729968, 4.370404, &
           58.50896, 2.870127, 0.8597835, 4.436456, 3.453213, &
           0.4651635, 4.199774, 3.297550, 0.1998602, 6.712620, &
           97.43101, 24.25369, 41.67469, 0.2340235, -0.6057560, &
           61.09523, 1.224794, 1.293125, 3.812594, 0.7535967, &
           8.970005, 8499.998, 0.8680952, -9.8761562E-03, 349.9999, &
           5.6327466E-02, 1.336010, 6.535113, 8.229236, 1.579276, &
           5596.842, 2.876048, -2.923217, 156.8400, 9.0139635E-02, &
           0.9942751/)


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

   print *, "TARGET", "VN", "DIFF"
   do nt = 1, NConf
      write (*,*) confnames(nt), target(nt), VN(nt), VN(nt) - target(nt)
   enddo

   T = 10.0
   RT = 0.85
   etol = 2.0e-3
   ns = 25
   nt = 150
!  mev=2000000
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
