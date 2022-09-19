program xyzensemble2defs
 use iso_fortran_env
 integer             :: i,j,k,l,h,hh,m,n,ierr,u,uu
 integer             :: ii,jj,kk,ll
 real                :: r = 1.0e12
 integer             :: num_args
 integer             :: n_atoms = 0
 integer             :: n_files = 0
 real,parameter      :: r_min_criteria_connectivity=0.15, lambda = 0.855
 integer,parameter   :: max_n_componets=1
 character(len=100)  :: CIFFilename=" ", DEFFilename="xxx"
 character(len=100)  :: filename=" "
 character(len=10)   :: string
 character(len=100)  :: line
 type  :: particle
  integer           :: element
  integer           :: type_
  character(len=2)  :: label_element
  character(len=4)  :: label
  character(len=4)  :: new_label="Xxxx"
  character(len=6)  :: label_from_CIFFile="Xxxxxx"
  character(len=50) :: topol_label="Xxxx"
  character(len=10) :: hybridization="Unknown"
  integer           :: degree
  real              :: charge
  real              :: radius
  real              :: mass
  integer           :: n_components
  character(5)      :: clabel(max_n_componets)
  real              :: xyzs(1:3,max_n_componets)
  real              :: xyzc(1:3,max_n_componets)
 end type
 type(particle),allocatable,dimension(:,:) :: atom
!
! arguments in line
 character(len=100),dimension(:), allocatable :: args
!
 real,allocatable    :: DistanceMatrix(:,:), Energy(:), MolFraction(:)
 logical,allocatable :: ConnectedAtoms(:,:)
 integer :: bond_types_max=100, bend_types_max=100, tors_types_max=100, impr_types_max=100
!
 logical              :: modify_topology_flag
! Element table and Matrix Topology:
 integer,allocatable  :: TopologyMatrix(:,:)
 integer,parameter    :: max_number_of_elements = 12
 integer              :: element_table(1:max_number_of_elements)
 ! C  H  O  N  P  S Zn Cd He Ar Xe ...
 ! 1  2  3  4  5  6  7  8  9 10 11 ...
 ! 6  1  8  7 15 16 30 48  2 18 54 ...
!C                      H                      O                      N
 element_table(1) = 6 ; element_table(2) = 1 ; element_table(3) = 8 ; element_table(4) = 7
!P                      S                      Zn                     Cd
 element_table(5) =15 ; element_table(6) =16 ; element_table(7) =30 ; element_table(8) =48
!He                     K                     Xe                     Cl
 element_table(9) = 2 ; element_table(10)=19 ; element_table(11)=54 ; element_table(12)=17
 !
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
    stop
   case ('-i','--xyz','-in')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
   case ('-S','--search-and-modify-topology')
    modify_topology_flag = .true.
  end select
 end do
 open(100,file=CIFfilename,status='old',iostat=ierr)
 if(ierr/=0) stop 'Error opening CIF file'
! 
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  n_files = n_files + 1
  read(line,*) n_atoms
  read(100,*,iostat=ierr) line
  do i=1,n_atoms
   read(100,'(a)',iostat=ierr) line
  end do
 end do read_cif
 allocate( atom(1:n_files,1:n_atoms) )
 allocate(ConnectedAtoms(n_atoms,n_atoms))
 allocate(DistanceMatrix(n_atoms,n_atoms))
 allocate( TopologyMatrix(n_atoms,11) )
 allocate(Energy(1:n_files))
 allocate(MolFraction(1:n_files))
 write(6,*) n_files, 'configurations'
 rewind(100)
 do i=1,n_files
  read(100,'(a)') line
  read(100,*) Energy(i), MolFraction(i) 
  select case(i)
   case(1)
    do j=1,n_atoms
     read(100,'(a)') line
     read(line,*) atom(i,j)%label, (atom(i,j)%xyzc(k,1), k=1,3), atom(i,j)%charge
     call CheckAtom(atom(i,j)%label,atom(i,j)%mass,atom(i,j)%radius,atom(i,j)%element,atom(i,j)%label_element)
    end do
   case default
    do j=1,n_atoms
     read(100,'(a)') line
     read(line,*) atom(i,j)%label, (atom(i,j)%xyzc(k,1), k=1,3)
    end do
  end select
 end do
 close(100)
! topology:
 DistanceMatrix=0.0
 ConnectedAtoms=.false.
 do i=1,n_atoms
  k=0
  do j=1,n_atoms
   r = 0.0
   if ( i/=j ) call Distance(atom(1,i)%xyzc(1:3,1),atom(1,j)%xyzc(1:3,1),r)
   DistanceMatrix(i,j)=r
   if( r>0.05 .and. r <= (atom(1,i)%radius + atom(1,j)%radius)/lambda ) then
    k=k+1
    ConnectedAtoms(i,j)=.true.
   end if
  end do
  atom(1,i)%degree=k
 end do
 kk=0 ! number of bonds
 do i=1,n_atoms
  do j=i+1,n_atoms
   if(ConnectedAtoms(i,j)) kk=kk+1
  end do
  write(string,'(a,i2)')atom(1,i)%label,i ; call StripSpaces(string)
  atom(1:n_files,i)%new_label=string(1:4)
 end do
 write(6,'(a)') 'Connectivity of each atom:'
 write(6,'(100(i2,1x))') ( atom(1,i)%degree, i=1,n_atoms )
 write(6,'(a)')' '
 write(6,'(a)')'Simulation input File'
 open(newunit=uu,file='simulation.input')
 write(uu,'(a)')'SimulationType                   MC'
 write(uu,'(a)')'NumberOfCycles                   2500000'
 write(uu,'(a)')'NumberOfInitializationCycles     0'
 write(uu,'(a)')'NumberOfEquilibrationCycles      0'
 write(uu,'(a)')'PrintPropertiesEvery             10'
 write(uu,'(a)')'PrintEvery                       10'
 write(uu,'(a)')' '
 write(uu,'(a)')'ChargeMethod                     Ewald'
 write(uu,'(a)')'Forcefield                       Local'
 write(uu,'(a)')'Framework                        0'
 write(uu,'(a)')'FrameworkName                    LOL_scaled'
 write(uu,'(a)')'UseChargesFromCIFFile            yes'
 write(uu,'(a)')'RemoveAtomNumberCodeFromLabel    yes'
 write(uu,'(a)')'UnitCells                        1 1 1'
 write(uu,'(a)')'ReplicaUnitCells                 2 2 2'
 write(uu,'(a)')'ExternalTemperature              300.0'
 write(uu,'(a)')' '
 !
 write(6,'(a)') 'Pseudoatoms:'
 open(newunit=u,file='pseudo_atoms.def')
 write(u,'(a)') '# number of pseudo atoms'
 write(u,'(i3)') n_atoms
 write(u,'(a)') '#type'
 do i=1,n_atoms
  write(u,'(a,1x,a,1x,a,1x,a,1x,i1,1x,5(f14.7,1x),i3,1x,i3,1x,a,1x,i1)') &
   atom(1,i)%new_label,'yes',atom(1,i)%label_element,atom(1,i)%label,0,atom(1,i)%mass,&
   atom(1,i)%charge,0.0,1.0,atom(1,i)%radius,atom(1,i)%degree,0,'absolute',0
 end do
 close(u)
 !
 write(6,'(a)') 'Molecule Definiton (RASPA):'
 do ii=1,n_files
  write(DEFFilename,'(a,a,i3,a)')filename(1:clen_trim(filename)),'_',ii-1,'.def' ; call StripSpaces(DEFFilename)
  write(6,'(a,1x,a)') 'writing:', DEFFilename
  open(newunit=u,file=DEFFilename)
  write(u,'(a,/,a,/,a,/,a,/,a,/,i3,/,a,/,i1,/,a,/,a,/,a,/,i3,/,a)')&
  '# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]',&
  '0.0',&
  '0.0',&
  '1.0',&
  '# Number Of Atoms',&
  n_atoms,&
  '# Number of groups',&
  1,&
  '# Flexibility (no)',&
  'rigid',&
  '# number of atoms',&
  n_atoms,&
  '# atomic positions'
  do i=1,n_atoms
   write(u,'(i3,1x,a,1x,3(f14.7))') i-1,atom(ii,i)%new_label,(atom(ii,i)%xyzc(k,1), k=1,3)
  end do
  !
  write(u,'(a)')'#'
  write(u,'(14i4)') 0, kk, (0, i=1,12)
  write(u,'(a)')'Bond stretch:'
  do i=1,n_atoms
   do j=i+1,n_atoms
    if(ConnectedAtoms(i,j)) write(u,'(i4,1x,i4,1x,a,1x,a,1x,a)') &
      i-1,j-1,'RIGID_BOND #',atom(1,i)%new_label,atom(1,j)%new_label
   end do
  end do
  write(u,'(a)')'Number of config moves'
  write(u,'(a)')'0'
  close(u)
  !
  write(uu,'(a,1x,i4,1x,a,a)')'Component',ii-1,'MoleculeName ',DEFFilename(1:Clen_trim(DEFFilename)-4)
  write(uu,'(a)')'  MoleculeDefinition             Local'
  write(uu,'(a)')'  WidomProbability               1.0'
  write(uu,'(a)')'  FugacityCoefficient            1.0'
  write(uu,*)    '# MolFraction', MolFraction(ii), '#', Energy(ii)  ! Only for Widom Test calculations
  write(uu,'(a)')'  ExtraFrameworkMolecule         no'
  write(uu,'(a)')'  CreateNumberOfMolecules        0'
  write(uu,'(a)')' '
 end do
 deallocate(Energy)
 deallocate(MolFraction)
 stop
 contains
!
 subroutine Distance(r_1,r_2,r)
  implicit none
  real,intent(in)  :: r_1(1:3), r_2(1:3)
  real,intent(out) :: r
  real             :: d(1:3) 
  integer          :: jj
  forall ( jj=1:3 )
   d(jj) = r_2(jj) - r_1(jj)
  end forall
  r = sqrt( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )
 end subroutine Distance
 !
 subroutine CheckAtom(Label,m,s,Z,Zlabel)
  implicit none
  character(len=4),intent(in)  :: Label
  real,intent(out)             :: m,s
  real,parameter               :: dummy_mass = 1.00794 ! (hydrogen atom mass)
  integer,intent(out)          :: Z
  character(len=2),intent(out) :: ZLabel
  ! Covalent Radious from:
  ! B. Cordero et al., Dalton Trans., 2008, 2832–2838
  ! doi: 10.1039/b801115j
  ! r_ij <= ( s_i + s_j ) / lambda
  ! lambda = 0.9  for all organic interactions
  ! lambda = 0.85 for metal-organic interactions
  ! we have selected: lambda=0.855 and
  ! the mass like the "conventional weight"
  select case(Label)
   case('H   ','H0  ':'H999')
    Z=1  ; m=dummy_mass ; s=0.31   ; Zlabel= ' H'
   case('K   ','K0  ':'K999')
    Z=19 ; m=39.0983 ; s=0.0       ; Zlabel= ' K'
   case('S   ','S0  ':'S999')
    Z=16 ; m=32.06 ; s=1.05        ; Zlabel= ' S'
   case('O   ','O0  ':'O999')
    Z=8  ; m=15.999 ; s=0.66       ; Zlabel= ' O'
   case('C   ','C0  ':'C999')
    Z=6  ; m=12.0107 ; s=0.75      ; ZLabel= ' C'
   case('N   ','N0  ':'N999')
    Z = 7 ; m = 14.007  ; s = 0.71 ; Zlabel= ' N'
   case('Zn  ','Zn0 ':'Zn99')
    Z = 30 ; m = 65.382 ; s = 1.22 ; Zlabel= 'Zn'
   case('P   ','P0  ':'P999')
    Z = 15 ; m = 39 ; s=1.07 ; Zlabel = ' P'
   case('Cl  ',' Cl ','Cl0 ':'Cl99')
    Z = 17 ; m = 35.453 ; s = 1.02 ; Zlabel= 'Cl'
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    STOP 'Atom unknowed'
  end select
 end subroutine checkatom
!
 subroutine StripSpaces(string)
  character(len=*) :: string
  integer :: stringLen
  integer :: last, actual
  stringLen = len (string)
  last = 1
  actual = 1
  do while (actual < stringLen)
   if (string(last:last) == ' ') then
    actual = actual + 1
    string(last:last) = string(actual:actual)
    string(actual:actual) = ' '
   else
    last = last + 1
    if (actual < last) actual = last
   endif
  end do
 end subroutine
!
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
!
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 subroutine print_help()
    print '(a)', '  -h, --help   print usage information and exit'
    print '(a)', '  -i, --xyz    input XYZ file of the conformer ensemble'
 end subroutine print_help
end program xyzensemble2defs
