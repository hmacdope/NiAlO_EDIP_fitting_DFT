module fit
   implicit none
   integer, parameter::np = 51
   integer, parameter::MaxAt = 300
   integer, parameter::MaxConf = 110

   logical, parameter::NoGas = .false., RefN2 = .false., TwoBodyOnly = .false., NoDissoc = .false.
   real(kind=8), dimension(3, MaxAt, MaxConf)::x
   real(kind=8), dimension(3, 3, MaxConf)::Q
   real(kind=8), dimension(MaxConf)::En, model, target, EAM, VN
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

      ! DFT energies
      print *, "DFT ENERGIES"
      do i = 1, NConf
         print *, En(i)
      end do
      print *, ""

      ! Read in EAM energies
      open (unit=10, file='EAM_energies.csv')
      read (10, *)
      do i = 1, NConf 
         read (10, *) str, EAM(i)
      end do
      close (10)

      print *, " EAM ENERGIES"
      do i = 1, NConf
         print *, EAM(i)
      end do
      print *, ""

      print *, "STARTING FIT ROUTINE"
      print *, ""
      print *, ""

      target = 0


      if (RefN2) then
         ! Ni+N
         target(1:4) = En(9:12) - (En(14) + En(31)*0.5) - (EAM(14) - EAM(9:12))
         ! Al+N
         target(5:8) = En(1:4) - (En(13) + En(31)*0.5) - (EAM(13) - EAM(1:4))
         ! NiAl+N
         target(9:12) = En(5:8) - (En(15) + En(31)*0.5) - (EAM(15) - EAM(5:8))
         ! N2
         target(13:16) = En(31:34) - (/2, 4, 4, 4/)*En(31)*0.5
         ! slabs: Al+N
         target(17:20) = En(17:20) - (En(16) + En(31)*0.5) - (EAM(16) - EAM(17:20))
         target(21:24) = En(21:24) - (En(16) + 2*En(31)*0.5) - (EAM(16) - EAM(21:24))
         ! slabs: Ni+N
         target(25:28) = En(26:29) - (En(25) + En(31)*0.5) - (EAM(25) - EAM(26:29))
         ! N2-short
         target(29) = En(35) - 2*En(31)*0.5
         ! slab-N2 scans:
         target(30:53) = En(36:59) - (En(16) + 2*En(31)*0.5)
         target(54:76) = En(60:82) - (En(25) + 2*En(31)*0.5)
         ! AlN:
         target(77:78) = En(83:84) - (2*En(13)/256.0_8 + 2*En(31)*0.5) - (2*EAM(13)/256.0_8 - EAM(30:31))
         ! N2 clusters
!  target(79:98)=En(85:104)- ...
      else
         ! Ni+N
         target(1:4) = En(9:12) - (En(14) + En(30)) - (EAM(14) - EAM(9:12))
         ! Al+N
         target(5:8) = En(1:4) - (En(13) + En(30)) - (EAM(13) - EAM(1:4))
         ! NiAl+N
         target(9:12) = En(5:8) - (En(15) + En(30)) - (EAM(15) - EAM(5:8))
         ! N2
         target(13:16) = En(31:34) - (/2, 4, 4, 4/)*En(30)
         ! slabs: Al+N
         target(17:20) = En(17:20) - (En(16) + En(30)) - (EAM(16) - EAM(17:20))
         target(21:24) = En(21:24) - (En(16) + 2*En(30)) - (EAM(16) - EAM(21:24))
         ! slabs: Ni+N
         target(25:28) = En(26:29) - (En(25) + En(30)) - (EAM(25) - EAM(26:29))
         ! N2-short
         target(29) = En(35) - 2*En(30)
!write(*,*)target(13:16),target(29)
!stop
         ! slab-N2 scans:
         target(30:53) = En(36:59) - (En(16) + 2*En(30))
         target(54:76) = En(60:82) - (En(25) + 2*En(30))
         ! AlN:
         target(77:78) = En(83:84) - (2*En(13)/256.0_8 + 2*En(30)) - (2*EAM(13)/256.0_8 - EAM(30:31))
         ! N2 clusters
         target(79:98) = En(85:104) - 10*En(30)
      end if
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
            if (lab(i1, i) /= 'N') cycle
            Z = 0
            do i2 = 1, NAt(i)
               if (lab(i2, i) == 'Ni') ty = 1
               if (lab(i2, i) == 'Al') ty = 2
               if (lab(i2, i) == 'O') ty = 3 ! NOTE EDITED
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
               if (lab(i2, i) == 'O') ty = 3 ! NOTE EDITED
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
      VN(1:4) = model(9:12)
      VN(5:8) = model(1:4)
      VN(9:12) = model(5:8)
      VN(13:16) = model(31:34)
      VN(17:24) = model(17:24)
      VN(25:28) = model(26:29)
      VN(29) = model(35)
      VN(30:76) = model(36:82)
      VN(77:78) = model(83:84)
      VN(79:98) = model(85:104)
      dE(1:98) = target(1:98) - VN(1:98)
!dE(3)=dE(3)/2.0
!dE(7)=dE(7)/20.0
      dE(7) = dE(7)/5.0
!dE(9)=dE(9)/5.0
!dE(11)=dE(11)/5.0
      dE(80) = dE(11)/30
      dE(1:78) = dE(1:78)/2.0
      if (NoGas) then
!  R=sum(dE(1:12)**2)+sum(dE(17:22)**2)+sum(dE(24:28)**2)
         R = sum(dE(1:12)**2) + sum(dE(17:22)**2) + sum(dE(24:28)**2) + dE(77)**2 + dE(78)**2
!  R=sum(dE(1:13)**2)+sum(dE(17:22)**2)+sum(dE(24:28)**2)
      else
!  R=sum(dE(13:16)**2+dE(29)**2)  ! wrong
         if (NoDissoc) then
            R = sum(dE(1:29)**2) + sum(dE(77:98)**2)
         else
            R = sum(dE(1:98)**2)
         end if
!  R=sum(dE(29:76)**2)+sum(dE(13:16)**2)
      end if
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

!par(1:24)=(/3.63208961,15.7711592,3.39307761,4.89977002E-02,8.99143410,6.43206263,1.83402944,&
!2.43297362,14.9088411,6.89796070E-05,3.18342257,2.52107453,0.796471894,2.06707790E-03,6.54625082,&
!8.81447697,2.52026415,0.905645967,29.9908390,1.87918830,0.218855426,19.1351109,7.70084858,2.47109160E-02/)

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

!par(11:15)=(/13.853,0.76143,4.7753,0.13113,1.1231/)
!par(22:24)=(/1.4956,7.2177,0.25882/)
!par(31:33)=(/28.877,3.3538,0.17964/)
!par(52:60)=(/2.3767,8.7378,0.47089,10.766,0.59636,-0.31442,9.9665,0.68064,0.23313/)

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

!do nt=1,3
!  do i=0,100
!    rr=1+(7-1)/100.0_8*i
!    write(20+nt,*)rr,par(nt)/rr**(7+floor(par(3+nt))) &
!                                +par(6+nt)*par(9)/rr*exp(-rr/par(10)) &
!                                -par(10+nt)/rr**4*exp(-rr/par(14)) &
!                                -par(14+nt)/rr**6
!  enddo
!enddo

   call LoadConfigs

   mf = Cost(par)
   do nt = 1, 29
      if (NoGas .and. any(nt == (/13, 14, 15, 16, 23, 29/))) then
         write (*, *)
      else
         write (*, *) target(nt), VN(nt), VN(nt) - target(nt)
      end if
   end do
   if (.not. NoGas) then
      write (*, *)
      do nt = 30, 76
         write (*, *) target(nt), VN(nt), VN(nt) - target(nt)
      end do
   end if
   do nt = 77, 78
      write (*, *) target(nt), VN(nt), VN(nt) - target(nt)
   end do
   write (*, *)
   do nt = 79, 98
      write (*, '(i4,3g16.8)') nt, target(nt), VN(nt), VN(nt) - target(nt)
   end do
!stop

   T = 10.0
   RT = 0.92
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
   do nt = 1, 29
      if (NoGas .and. any(nt == (/13, 14, 15, 16, 23, 29/))) then
         write (*, *)
      else
         write (*, '(i4,3g16.8)') nt, target(nt), VN(nt), VN(nt) - target(nt)
      end if
   end do
   if (.not. NoGas) then
      write (*, *)
      do nt = 30, 76
         write (*, '(i4,3g16.8)') nt, target(nt), VN(nt), VN(nt) - target(nt)
      end do
   end if
   do nt = 77, 78
      write (*, '(i4,3g16.8)') nt, target(nt), VN(nt), VN(nt) - target(nt)
   end do
   write (*, *)
   do nt = 79, 98
      write (*, '(i4,3g16.8)') nt, target(nt), VN(nt), VN(nt) - target(nt)
   end do

!do nt=1,3
!  do i=0,100
!    rr=1+(5-1)/100.0_8*i
!    write(23+nt,*)rr,par(nt)/rr**(7+floor(par(3+nt))) &
!                                +14.3996434*par(6+nt)*par(9)/rr*exp(-rr/par(10)) &
!                                -par(10+nt)/rr**4*exp(-rr/par(14)) &
!                                -par(14+nt)/rr**6
!  enddo
!enddo

end program fitter
