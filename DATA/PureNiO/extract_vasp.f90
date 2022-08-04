program ff
implicit none
integer,parameter::MaxSpec=10,maxat=1000
character(len=100)::inline
character(len=2),dimension(MaxSpec)::Labels
character(len=10)::a,b
integer,dimension(MaxSpec)::NTemp
real,dimension(3,maxat)::x0
real,dimension(3)::x,a1,a2,a3
real(kind=8)::en,edisp
integer::i,j,m
integer::ns,n,badness
logical::vdw,anyX,goodX
vdw=.false.
open(unit=10,file='OUTCAR')
 ! this was for <6.2
!do i=1,15
!  read(10,'(a)')inline
!  if (inline(2:6)=='INCAR') exit
!enddo
!if (i==16) then
!  write(*,*)'Missed INCAR'
!  stop
!endif
!do i=1,MaxSpec
!  read(10,*)a,b,Labels(i)
!  if (Labels(i)(2:2)=='_') Labels(i)(2:2)=' '
!  if ((i>1).and.(Labels(i)==Labels(1)).or.(a/='POTCAR:')) exit
!enddo
do
  read(10,'(a100)')inline
  if (inline(2:7)=='POTCAR') exit
enddo
read(inline,*)a,b,Labels(1)
if (Labels(1)(2:2)=='_') Labels(1)(2:2)=' '
do i=2,MaxSpec
  read(10,*)a,b,Labels(i)
  if (Labels(i)(2:2)=='_') Labels(i)(2:2)=' '
  if ((Labels(i)==Labels(1)).or.(a/='POTCAR:')) exit
enddo

ns=i-1
do
  read(10,'(a)')inline
  if (inline(4:16)=='ions per type') exit
enddo
read(inline(25:),*)NTemp(1:ns)
n=sum(NTemp(1:ns))
if (n>maxat) stop 'maxat too small'
do i=1,ns
  write(*,*)Labels(i),NTemp(i)
enddo
open(unit=11,file='out.dat')
j=0
edisp=0
anyX=.false.
goodX=.false.
en=0
do
  read(10,'(a)',iostat=badness)inline
  if (badness/=0) exit
  if (inline(18:22)=='TOTEN') then
    if (inline(39:40)=='**') then
      en=0
    else
      read(inline(27:43),*)en
    endif
  endif
  if (inline(3:6)=='IVDW') vdw=.true.
  if (inline(2:7)=='Edisp') read(inline(14:),*)edisp
  if (inline(7:28)=='direct lattice vectors') then
    read(10,'(3X,3f13.9)')a1
    read(10,'(3X,3f13.9)')a2
    read(10,'(3X,3f13.9)')a3
  endif
  if (inline(2:42)=='position of ions in cartesian coordinates') then
    do i=1,n
      read(10,'(X,3f12.8)')x0(:,i)
    enddo
    anyX=.true.
  endif
  if (inline(2:9)=='POSITION') then
    j=j+1
    write(11,*)n+8
    write(11,*)'Energy:',en+edisp
    read(10,'(a)')inline
    write(11,*)'X',0.0,0.0,0.0
    write(11,*)'X',a1
    write(11,*)'X',a2
    write(11,*)'X',a3
    write(11,*)'X',a1+a2
    write(11,*)'X',a1+a3
    write(11,*)'X',a2+a3
    write(11,*)'X',a1+a2+a3
    m=1
    do i=1,n
      read(10,*)x
      if (i>sum(NTemp(1:m))) m=m+1
      write(11,*)Labels(m),x
    enddo
    goodX=.true.
  endif
enddo
if (anyX.and..not.goodX) then
  write(11,*)n+8
  write(11,*)'Energy:',en+edisp
  write(11,*)'X',0.0,0.0,0.0
  write(11,*)'X',a1
  write(11,*)'X',a2
  write(11,*)'X',a3
  write(11,*)'X',a1+a2
  write(11,*)'X',a1+a3
  write(11,*)'X',a2+a3
  write(11,*)'X',a1+a2+a3
  m=1
  do i=1,n
    if (i>sum(NTemp(1:m))) m=m+1
    write(11,*)Labels(m),x0(:,i)
  enddo
  j=-1
endif
close(10)
close(11)
if (vdw) write(*,*)'Warning: extra vdW contribution truncated'
write(*,'(i4,'' x '',i3,'' ('',i3,'')'')')j,n,n+8
end
