
program pdb2angles
implicit none
integer, parameter :: mat=1000  ! atom number limit
character(80) atmp,inf,outf
integer io,i,iat,s2i,it(mat),resid(mat),nat,mres,is,strands(100),k,send,sstart,ctype(mat)
real(8) xyz(3,mat),s2r,dih,dihgrad,FK
character(4) atype
integer phos(mat),o3p(mat),o5p(mat),c3p(mat),c4p(mat),c5p(mat),n1(mat),n9(mat),c2(mat),c4(mat),o4p(mat),c1p(mat)
real(8), allocatable :: alpha(:),beta(:),gama(:),eps(:),zeta(:),delta(:),chi(:)
character(3) resname(100)

call getarg(1,inf)
inf=trim(inf)
!call getarg(2,outf)
!outf=trim(outf)
!print*,trim(inf),trim(outf)

call printhead

is=1
strands(is)=1

io=33
open(unit=io,file=inf)
do
 read(io,'(a)',end=666) atmp
 if(index(atmp,'ATOM').ne.0) then
! atom number
   iat=s2i(atmp(7:11)) 
! coordinates
   xyz(1,iat)=s2r(atmp(31:38))
   xyz(2,iat)=s2r(atmp(39:46))
   xyz(3,iat)=s2r(atmp(47:54))
! residue sequence number
   resid(iat)=s2i(atmp(23:26))
! residue name
   resname(resid(iat))=atmp(18:20)
  ! pyrimidines
  if(index(atmp(18:20),"C").ne.0) ctype(resid(iat))=0
  if(index(atmp(18:20),"T").ne.0) ctype(resid(iat))=0
  if(index(atmp(18:20),"U").ne.0) ctype(resid(iat))=0
  ! purines
  if(index(atmp(18:20),"A").ne.0) ctype(resid(iat))=1
  if(index(atmp(18:20),"G").ne.0) ctype(resid(iat))=1
  if(index(atmp(18:20),"X").ne.0) ctype(resid(iat))=1
! atom type
   atype=atmp(13:16)
   if(atype==" O3'".or.atype==" O3*") then
     it(iat)=1
     o3p(resid(iat))=iat
   endif
   if(atype==" P  ") then
      it(iat)=2
      phos(resid(iat))=iat
   endif
   if(atype==" O5'".or.atype==" O5*") then
      it(iat)=3
      o5p(resid(iat))=iat
   endif
   if(atype==" C5'".or.atype==" C5*") then
      it(iat)=4
      c5p(resid(iat))=iat
   endif
   if(atype==" C4'".or.atype==" C4*") then
      it(iat)=5
      c4p(resid(iat))=iat
   endif
   if(atype==" C3'".or.atype==" C3*") then
       it(iat)=6
       c3p(resid(iat))=iat
   endif
   if(atype==" N1 ") then
       it(iat)=7
       n1(resid(iat))=iat
   endif
   if(atype==" N9 ") then
       it(iat)=8
       n9(resid(iat))=iat
   endif
   if(atype==" C2 ") then
       it(iat)=9
       c2(resid(iat))=iat
   endif
   if(atype==" C4 ") then
       it(iat)=10
       c4(resid(iat))=iat
   endif
   if(atype==" O4'".or.atype==" O4*") then
      it(iat)=11
      o4p(resid(iat))=iat
   endif
   if(atype==" C1'".or.atype==" C1*") then
      it(iat)=12
      c1p(resid(iat))=iat
   endif
!   print*,trim(atmp)
!   print*, iat,it(iat),resid(iat)
!   print*,'--'
 endif
 if(index(atmp,'TER').ne.0) then
  strands(is)=resid(iat) !end of last strand
 print*,strands(is)
 is=is+1
 endif
enddo
666 close(io)

if(is==1) then 
 strands(is)=resid(iat) !end of last strand
 is=is+1
endif

nat=iat
mres=maxval(resid)
write(*,*)''
write(*,*)''
write(*,*) nat,' number of atoms / ',mres,' residues / ',is-1,' strands'
write(*,*)''
write(*,*)''
allocate(alpha(mres),beta(mres),eps(mres),zeta(mres),delta(mres),gama(mres),chi(mres))

! From X3DNA:
! Note: alpha:   O3'(i-1)-P-O5'-C5'
!      beta:    P-O5'-C5'-C4'
!      gamma:   O5'-C5'-C4'-C3'
!      delta:   C5'-C4'-C3'-O3'
!      epsilon: C4'-C3'-O3'-P(i+1)
!      zeta:    C3'-O3'-P(i+1)-O5'(i+1)

!      chi for pyrimidines(Y): O4'-C1'-N1-C2
!          chi for purines(R): O4'-C1'-N9-C4
! Y= C,T,U
! R=A,G,X

write(*,'(7x,a,1x,a,1x,7(a,2x) )') 'resid','resname',' alpha-1',' beta-1',' gamma ',' delta ','epsilon','  zeta ','| chi'
sstart=1
do k=1,is ! loop over strands
send=strands(k)
do i=sstart,send ! loop over res in strands

if(i.ne.sstart)    call dihed(xyz,o3p(i-1),phos(i),o5p(i),c5p(i),dih,alpha(i)) !alpha
if(i.ne.sstart)    call dihed(xyz,phos(i),o5p(i),c5p(i),c4p(i),dih,beta(i)) !beta
!    call dihed(xyz,o3p(i-1),phos(i),o5p(i),c5p(i),dih,alpha(i)) !alpha
!    call dihed(xyz,phos(i),o5p(i),c5p(i),c4p(i),dih,beta(i)) !beta
              call dihed(xyz,o5p(i),c5p(i),c4p(i),c3p(i),dih,gama(i)) !gamma
              call dihed(xyz,c5p(i),c4p(i),c3p(i),o3p(i),dih,delta(i)) !delta
if(i.ne.send) call dihed(xyz,c4p(i),c3p(i),o3p(i),phos(i+1),dih,eps(i)) !epsilon
if(i.ne.send) call dihed(xyz,c3p(i),o3p(i),phos(i+1),o5p(i+1),dih,zeta(i)) !zeta
! call dihed(xyz,c4p(i),c3p(i),o3p(i),phos(i+1),dih,eps(i)) !epsilon
! call dihed(xyz,c3p(i),o3p(i),phos(i+1),o5p(i+1),dih,zeta(i)) !zeta


if(ctype(i)==0) call dihed(xyz,o4p(i),c1p(i),n1(i),c2(i),dih,chi(i)) !chi Y
if(ctype(i)==1) call dihed(xyz,o4p(i),c1p(i),n9(i),c4(i),dih,chi(i)) !chi R

!call dihed(xyz,o4p(i),c1p(i),n1(i),c2(i),dih,chi(i)) !chi Y
!print*,chi(i)
!call dihed(xyz,o4p(i),c1p(i),n9(i),c4(i),dih,chi(i)) !chi R
!print*,chi(i)


write(*,'(7x,2x,I3,3x,a3,3x,7(F7.2,2x) )') i,resname(i),alpha(i),beta(i),gama(i),delta(i),eps(i),zeta(i),chi(i)
enddo
sstart=strands(k)+1
write(*,*)''
enddo


!-----------------------------------------------------------------------
! print output
io=44
open(unit=io,file='atomindex.dat')
FK=0.01d0
sstart=1
do k=1,is ! loop over strands
send=strands(k)
do i=sstart,send ! loop over res in strands
write(io,'(a,1x,I3,2x,a)')'# resid/name',i,resname(i)
if(i.ne.sstart) then
 write(io,*) '# alpha'
 write(io,'(a,4(I4),F5.2)') 'dihed ', o3p(i-1),phos(i),o5p(i),c5p(i),FK
endif
if(i.ne.sstart) then
 write(io,*) '# beta'
 write(io,'(a,4(I4),F5.2)') 'dihed ', phos(i),o5p(i),c5p(i),c4p(i),FK
endif
 write(io,*) '# gamma'
write(io,'(a,4(I4),F5.2)') 'dihed ', o5p(i),c5p(i),c4p(i),c3p(i),FK
 write(io,*) '# delta'
write(io,'(a,4(I4),F5.2)') 'dihed ', c5p(i),c4p(i),c3p(i),o3p(i),FK
if(i.ne.send) then
 write(io,*) '# eps'
 write(io,'(a,4(I4),F5.2)') 'dihed ', c4p(i),c3p(i),o3p(i),phos(i+1),FK
endif
if(i.ne.send) then
 write(io,*) '# zeta'
 write(io,'(a,4(I4),F5.2)') 'dihed ', c3p(i),o3p(i),phos(i+1),o5p(i+1),FK
endif
write(io,*) '# chi'
if(ctype(i)==0)write(io,'(a,4(I4),F5.2)') 'dihed ',o4p(i),c1p(i),n1(i),c2(i),FK !chi Y
if(ctype(i)==1)write(io,'(a,4(I4),F5.2)') 'dihed ',o4p(i),c1p(i),n9(i),c4(i),FK !chi R
write(io,*)''
enddo
sstart=strands(k)+1
write(io,'(a)')''
enddo

close(io)

end

!-----------------------------------------------------------------------
! String to integer
integer pure function s2i(a)
implicit none
character(*), intent(in):: a
read(a,'(i)') s2i
return
end function


! String to real(8)
real(8) pure function s2r(a)
implicit none
character(*), intent(in):: a
read(a,'(f)') s2r
return
end function



subroutine dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
implicit none
real(8), parameter:: pi = 3.141592653589793d0
integer aa,bb,cc,dd
real(8) b1(3),b2(3),b3(3),n1(3),n2(3)
real(8) un1(3),un2(3),ub2(2),m(3),um(3),dix,diy
real(8) dih,dihgrad,xyz(3,*)

 b1=xyz(1:3,aa)-xyz(1:3,bb)
 b2=xyz(1:3,bb)-xyz(1:3,cc)
 b3=xyz(1:3,cc)-xyz(1:3,dd)

 ! normal of the planes
 call cross_prod(n1,b1,b2)
 call cross_prod(n2,b2,b3)

 call unitvec(n1,un1)
 call unitvec(n2,un2)
 call unitvec(b2,ub2)

 call cross_prod(m,un1,ub2)
 call unitvec(m,um)

 dix=DOT_PRODUCT(un1,un2)
 diy=DOT_PRODUCT(um,un2)

 dih=atan2(diy,dix)

!  Quadrant    Angle              sin    cos    tan
!----------------------------------------------------
!  I           0    < α < π/2     > 0    > 0    > 0       
!  II          π/2  < α < π       > 0    < 0    < 0
!  III         π    < α < 3π/2    < 0    < 0    > 0
!  IV          3π/2 < α < 2π      < 0    > 0    < 0
! atan2(0,1) =   0
! atan2(1,0) =   pi/2
! atan2(-1,0) = -pi/2
! atan2(0,-1) =  pi


! give results in 0 to 360 degree
 if(dih<0.0d0) dih=dih+pi*2
 dihgrad=dih*180.0d0/pi
end subroutine

real*8 function di360(x)
implicit none
real(8), parameter:: pi = 3.141592653589793d0
real(8) x
if(x<0.0d0) x=x+pi*2
di360=x*180.0d0/pi
end function


subroutine unitvec(x,e)
implicit none
real(8) x(3),e(3),t(3)
t=DOT_PRODUCT(x,x)
!t=x(1)**2+x(2)**2+x(3)**2
e=x/sqrt(t)
end

subroutine veclen(x,v)
implicit none
real(8) x(3),v
v=dsqrt(dot_product(x,x))
end subroutine


subroutine veclen2(a,b,v)
implicit none
real(8) a(3),b(3),v,x(3)
x=a-b
v=dsqrt(dot_product(x,x))
end subroutine


subroutine cross_prod(y,x2,x3)
implicit none
real(8) y(3),x2(3),x3(3)
  y(1) =  x2(2)*x3(3) - x3(2)*x2(3)
  y(2) = -x2(1)*x3(3) + x3(1)*x2(3)
  y(3) =  x2(1)*x3(2) - x3(1)*x2(2)
end subroutine


subroutine printhead
write(*,'(2x,a)') ''
write(*,'(2x,a)') ''
write(*,'(2x,a)') '****************************'
write(*,'(2x,a)') '* Backbone angle calculator *'
write(*,'(2x,a)') '****************************'
write(*,'(2x,a)') ''
!write(*,'(2x,a)') "NOTE: convert atomtypes from * to '"


end subroutine
