program histogram
 implicit none
 integer            :: n_datos = 0
 integer            :: n_boxs  = 60
 integer            :: ierr,i,j,k
 real               :: max_,min_,suma
 character (len=80) :: line
 real, allocatable  :: values(:),delta(:)
 logical            :: normalise = .false.
!
 open(111,file="input",status='old',action='read',iostat=ierr)
 fileopen: if( ierr == 0) then
  read_: do
    read (111,'(a)',iostat=ierr) line
    if( ierr /= 0 ) exit
    n_datos = n_datos + 1
  end do read_
  rewind( 111 )
 endif fileopen
 allocate(values(1:n_datos) ,stat=ierr)
 allocate(delta(0:n_boxs+1) ,stat=ierr)
 if(ierr/=0) stop '[error] variables sin alicatar en memoria.'
 datas: do i=1,n_datos
   read(111,*) values(i)
 end do datas
 max_  = maxval(values)
 min_  = minval(values)
 delta(0) = min_
 do j=1,n_boxs+1
  delta(j)=delta(j-1)+(max_-min_)/real(n_boxs)
 enddo
 call make_histogram(values,delta,n_datos,n_boxs+1)
 close(111)
 deallocate(values)
 deallocate(delta)
 contains
 subroutine make_histogram(data,bound,j,k)
   implicit none
   integer :: j,k
   real    :: ave,adev,sdev,var,skew,curt
   real,   intent(in) :: data(1:j)
   real,   intent(in) :: bound(0:k)
   real    :: histo(0:k), suma 
   integer :: i
   suma=0.0
   do i = 1,k
    histo(i) = count( data <= bound(i) .and. data >= bound(i-1))
   enddo
   ave = sum(histo)
   do i=1,k
    histo(i)=histo(i)/ave
    suma = bound(i)*histo(i) + suma
    write(6,*) bound(i),histo(i), suma
   end do
   call moment(data,j,ave,adev,sdev,var,skew,curt)
   write(6,*)'# ave,sdev,skew,curt:'
   write(6,*)'# moments:',ave,sdev,skew,curt
   write(6,*)'# adev,var:'
   write(6,*)'# deviation:',adev,var
   ! <E>= sum(E*p(E))/sum(p(E))
   write(6,*)'# Boltzmann:', suma
 end subroutine make_histogram
!
 subroutine moment(data,n,ave,adev,sdev,var,skew,curt)
  implicit none
! numerical recipes (fortran 90), pp 607-608.
! given an array of data(1:n), its returns its mean ave, average deviation adev,
! standar deviation sdev, variance var, skewness skew, and kurtosis curt.
  integer :: n,j
  real    :: adev,ave,curt,sdev,skew,var,data(n)
  real    :: p,ep
  real    :: s = 0.0
  if (n<=1) print*,'n must be at least 2 in moment'
  do j=1,n
   s=s+data(j)
  end do
  ave  = s/real(n)
  adev = 0.0
  var  = 0.0
  skew = 0.0
  curt = 0.0
  ep   = 0.0
  storage: do j=1,n
     s    = data(j) - ave
     ep   = ep + s
     adev = adev + abs(s)
     p    = s*s
     var  = var + p
     p    = p*s
     skew = skew + p
     p    = p*s
     curt = curt + p
  end do storage
  adev = adev/real(n)
  var  = (var-ep*ep/real(n))/(n-1)
  sdev = sqrt(var)
  if(var/=0.0)then
    skew = skew/(n*adev**3)
    curt = curt/(n*var*var)-3.0
  else
   skew=0.0
   curt=0.0
  end if
  return
 end subroutine
end program histogram
