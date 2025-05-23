program pep_loci
  implicit none

  character(len=200) :: pdbfile,outfile,rwpdbfile,name,num
  integer :: narg,cptArg,nframes,natom,nres,pepnr,pepna,lipna
  character(len=1) :: connn
  character(len=4) :: lipc,waterc,pepc
  logical :: lookForpdbfile=.false.
  logical :: lookForoutfile=.false.
  logical :: lookForpepc=.false.
  logical :: lookForlipc=.false.
  logical :: lookForwaterc=.false.
  logical :: lookFornframes=.false.
  logical :: fileExist
	 
!Check if arguments are found
  narg=command_argument_count()
	 
  if(narg>0)then
!loop across options
    do cptArg=1,narg
      call get_command_argument(cptArg,name)
      print*, name
      select case(adjustl(name))
!First known args
        case("-i")
          lookForpdbfile = .true. !change logical value
        case("-o")
          lookForoutfile = .true.
        case("-nf")
          lookFornframes = .true.
        case("-pepc")
          lookForpepc = .true.
        case("-lipc")
          lookForlipc = .true.
        case("-waterc")
          lookForwaterc = .true.
        case default
!Treat the second arg of a serie
          if(LookForpdbfile)then
            pdbfile=adjustl(name) !assign a value to pedfile
            inquire(file=pdbfile,exist=fileExist)!check if it exist
            if(.not.fileExist)then
              write(*,*)'file ',trim(pdbfile),' not found'
              stop
            endif
            LookForpdbfile=.false. !put the logical variable to its initial value
          elseif(LookForoutfile)then
            outfile=adjustl(name)
            inquire(file=outfile,exist=fileExist)
            if(fileExist)then
              write(*,*)'file ',trim(outfile),' already exists'
              print*, "continue?(y/n)"
              read*, connn
              if (connn=="n") stop
            endif
            LookForoutfile=.false.
          elseif(LookFornframes)then
            num=adjustl(name)
            read(num,*) nframes
            LookFornframes=.false.
          elseif(LookForpepc)then
            pepc=adjustl(name)
            LookForpepc=.false.
          elseif(LookForwaterc)then
            waterc=adjustl(name)
            LookForwaterc=.false.
          elseif(LookForlipc)then
            lipc=adjustl(name)
            LookForlipc=.false.
          else
            write(*,*)"Option ",adjustl(name),"unknown"
          endif
      end select
    end do
  endif
  
  print*,"pepc",pepc
  print*,"lipc",lipc
  print*,"waterc",waterc

  call pdbinform(pdbfile,natom,rwpdbfile,pepc,pepnr,pepna,lipc,lipna)

!  call interactions(rwpdbfile,natom,pepc,waterc,lipc,pepnr,pepna,nframes,outfile)

  call cmdistances(rwpdbfile,natom,pepc,lipc,pepna,lipna,pepnr,nframes,outfile)

end program

subroutine pdbinform(pdb,j,rwpdb,pepc,pepnr,pepna,lipc,lipna)
  implicit none

! leitura do pdb referência
  character(len=200), intent(in) :: pdb
  integer :: n_atom, nres,pos
  real :: xr, yr, zr, O, B
  character(len=4) :: nameatom,namechain,pepc,lipc
  character(len=6) :: ATOM
  character(len=3) :: res
  character(len=1) :: chain, charge, atomo,conf
  character(len=200) :: remark,cryst,rwpdb
  character(len=6) :: atomtype
! cálculo interno
  integer, intent(out) :: j,pepnr,pepna,lipna

  pepnr=0
! Abre o pdb de referência
  open (unit=9, file=pdb)
  rwpdb=pdb
  pos=index(rwpdb,".pdb")
  rwpdb=rwpdb(1:pos)//"rw.pdb"
  atomtype = ""

! Abre o novo pdb-só coordenadas
  open (unit=7, file=rwpdb)

  j=0
  do while (atomtype.ne."END   ")
    read(unit=9, fmt='(A6)'), atomtype
    if (atomtype.eq."CRYST1") then
      backspace(9)
      read(unit=9, fmt='(A)'), cryst
      write(unit=7, fmt='(A)'), cryst
    end if
    if (atomtype.eq."ATOM  " .or. atomtype.eq."HETATM") then 
      backspace(9)
      read(unit=9,fmt='(A6, I5, 1X, A4, A1, A3, 1X, A1, I4, 4X, F8.3, F8.3, F8.3, &
        &F6.2, F6.2, 6x,A4, A2, A2)'), ATOM, n_atom, nameatom,conf,&
        res, chain, nres,xr, yr, zr, O,&
        B, namechain,atomo, charge
      if (namechain.eq.pepc) pepnr=nres
      if (namechain.eq.pepc) pepna=n_atom
      if (namechain.eq.lipc) lipna=n_atom
      if (conf.eq."A") then
         j=j+1
        conf=" "
        n_atom=j
        write(unit=7,fmt='(A6, I5, 1X, A4, A1, A3, 1X, A1, I4, 4X, F8.3, F8.3, F8.3, &
          &F6.2, F6.2, 6x,A4, A2, A2)'), atomtype, n_atom, nameatom,conf,&
          res, chain, nres,xr, yr, zr, O,&
          B, namechain,atomo, charge
      elseif (conf.ne."B") then
        j=j+1
        n_atom=j
        write(unit=7,fmt='(A6, I5, 1X, A4, A1, A3, 1X, A1, I4, 4X, F8.3, F8.3, F8.3, &
          &F6.2, F6.2, 6x,A4, A2, A2)'), atomtype, n_atom, nameatom,conf,&
          res, chain, nres,xr, yr, zr, O,&
          B, namechain,atomo, charge
      end if
    end if


  end do
  write(unit=7, fmt='(A)'), "END   "
  
  close (9)
  close (7)

end subroutine

subroutine cmdistances(pdbfile,natom,pepc,lipc,pepna,lipna,pepnr,nframes,outfile)
  implicit none
  
! para a leitura do pdb
  character(len=200) :: remark,pdbfile,dcdfile,outfile
  character(len=4),dimension(natom) :: ATOM,nameatom,namechain
  character(len=3),dimension(natom) :: res
  character(len=1),dimension(natom) :: chain,nchain,atomtype
  integer,dimension(natom) :: n_atom,nres
  real,dimension(natom) :: O,B,x,y,z
  integer :: i
! leitura do dcd
  logical :: unitcell=.true.
  double precision :: sides(3)
  integer :: error,ioerror,stride=1,firstatom,lastatom,firstframe,lastframe
  character(len=4) :: dummyc
  integer :: dummyi,index,j,ntotframes,ntotat
  real :: dummyr
  double precision :: dummyd
! calculo intero
  integer :: a,pepna,lipna,pepnr,nlip,natom,nframes,k,l,nrot
  real :: soma,f,time
  character(len=4) :: dummy,pepc,lipc
  real,dimension(natom) :: m,distlipcm
  real :: masspep,masslip,cmlipx,cmlipy,cmlipz,cmpepx,cmpepy,cmpepz,micrad,cmdist,micd,var,dp,time
  logical,dimension(pepna) :: mainchain
  real,dimension(pepnr) :: cadist
  real,dimension(3) :: d,cem
  real,dimension(natom) :: xc,yc,zc
 
  natom=lipna

  i=0
  print*, "lendo pdb de referência..."
  open (unit=9, file=pdbfile, status="old")
!  read(unit=9, fmt='(A)'),  remark    ! lê a primeira linha de observação do pdb
  do i=1, natom, 1
    read(unit=9,fmt='(A4,1X,I6,1X,A4,1X,A3,1X,A1,1X,I3,4X,F8.3,F8.3,F8.3,2X,&
      &F4.2, 2X, F4.2, 6X, A4, A2,A2 )') ATOM(i),n_atom(i),nameatom(i),&
      res(i),chain(i),nres(i),x(i),y(i),z(i),O(i),B(i),namechain(i),nchain(i),atomtype(i)
  end do
  close (9)

  do a=1, natom
    dummy = nameatom(a)
    dummy = ADJUSTL(dummy)
    if (dummy(1:1) == "O") then
      m(a) = 15.99940
    else if (dummy(1:1) == "N") then
      m(a) = 14.00670
    else if (dummy(1:1) == "C") then
      m(a) = 12.01100
    else if (dummy(1:1) == "S") then
      m(a) = 32.06000
    else if (dummy(1:1) == "H") then
      m(a) = 1.00800
    end if
  end do

  masspep=0.
  masslip=0.
  micrad=0.
  nlip=0
  do a=1, pepna
    if (nameatom(a).eq." N  " .or. &
       &nameatom(a).eq." CA " .or. &
       &nameatom(a).eq." C  " .or. &
       &nameatom(a).eq." O  ") then
      mainchain(a)=.true.
      masspep=masspep+m(a)
    end if
  end do
  do i=pepna+1,lipna
    dummy = nameatom(i)
    dummy = ADJUSTL(dummy)
    if (dummy(1:1) .eq. "S" .or. dummy(1:1).eq."P") nlip=nlip+1
    if (dummy(1:1) .ne. "H") masslip=masslip+m(i)
  end do
  
  
  soma=0
  open (unit=8,file=outfile)
  write(unit=8,fmt='(a)'), "Distances to the micelle centre of mass."
  write(unit=8,fmt='(a)'), "micelle ray - peptide centre of mass - each residue"
  
!!loop para concatenar coordenadas na mesma linha na trajetória final
  print*, "lendo coordenadas..."
  open(unit=7,file="dcdfiles.dat")
  read(unit=7,fmt='(A)'),dcdfile
!!loop para ler vários dcds
  do while((dcdfile.ne."exit") .and. (soma<nframes))
!! checando número de frames e número de átomos da trajetória
    open(11,file=dcdfile,action='read',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      print*, "erro ao abrir o arquivo dcd"
      return
    end if
    read(11) dummyc, ntotframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(11) dummyi, dummyr
    read(11) ntotat
    print*, ntotframes,ntotat,dcdfile,pepnr,pepna
!! checado

    firstatom = 1 
    lastatom = natom
!    index = 0
    do f=1, ntotframes
      print*,f+soma
      cmlipx=0.
      cmlipy=0.
      cmlipz=0.
      micrad=0.
      cmpepx=0.
      cmpepy=0.
      cmpepz=0.
      
    
!      index = index + 1
!      j = 3*(index-1)
      if(unitcell) read(11) sides(1), dummyd, sides(2), dummyd, dummyd, sides(3)
!      iatom = natom*(index - 1)                                           
      read(11) (dummyr, j = 1, firstatom - 1 ), (x(j), j = 1, natom)
      read(11) (dummyr, j = 1, firstatom - 1 ), (y(j), j = 1, natom)            
      read(11) (dummyr, j = 1, firstatom - 1 ), (z(j), j = 1, natom) 
      
      
      do i=1,pepna
        if (mainchain(i)) then
          cmpepx=cmpepx+m(i)*x(i)
          cmpepy=cmpepy+m(i)*y(i)
          cmpepz=cmpepz+m(i)*z(i)
        end if
      end do
      cmpepx=cmpepx/masspep
      cmpepy=cmpepy/masspep
      cmpepz=cmpepz/masspep
            
      do i=pepna+1,lipna
        dummy = nameatom(i)
        dummy = ADJUSTL(dummy)
        if (dummy(1:1) .ne. "H") then
          cmlipx=cmlipx+m(i)*x(i)
          cmlipy=cmlipy+m(i)*y(i)
          cmlipz=cmlipz+m(i)*z(i)
        end if
      end do
      cmlipx=cmlipx/masslip
      cmlipy=cmlipy/masslip
      cmlipz=cmlipz/masslip

      cmdist=sqrt((cmpepx-cmlipx)**2 + (cmpepy-cmlipy)**2 + (cmpepz-cmlipz)**2)

      do i=1,pepna
        if (nameatom(i) .ne." CA ") then
          cadist(nres(i))=sqrt((x(i)-cmlipx)**2 + (y(i)-cmlipy)**2 + (z(i)-cmlipz)**2)
        end if
      end do

      var=0.
      do i=pepna+1,lipna
        if (nameatom(i) .eq." P  " .or.nameatom(i) .eq." S  ") then
          distlipcm(i)=(x(i)-cmlipx)**2 + (y(i)-cmlipy)**2 + (z(i)-cmlipz)**2
          micrad=micrad+distlipcm(i)
        end if
      end do
      micrad=sqrt(micrad/nlip)
      do i=pepna+1,lipna
        if (nameatom(i) .eq." P  " .or.nameatom(i) .eq." S  ") var = var + (sqrt(distlipcm(i))-micrad)**2
      end do
      var=var/nlip
      dp=sqrt(var)
      time=(soma+f)/100

      write(unit=8,fmt='(f6.2,2x,f6.2,2x,f6.2,2x,f6.2,*(2x,f6.2))'), time,micrad,dp,cmdist,(cadist(j),j=1,pepnr)
         
    end do
    soma=soma+ntotframes
    read(unit=7,fmt='(A)'),dcdfile
    close(11)
    
  end do
  soma=0
  close(7)
  close(8)
      
end subroutine

subroutine interactions(pdbfile,natom,pepc,waterc,lipc,pepnr,pepna,nframes,outfile)
  implicit none

! para a leitura do pdb
  character(len=200) :: remark,pdbfile,dcdfile,outfile
  character(len=4),dimension(natom) :: ATOM,nameatom,namechain
  character(len=3),dimension(natom) :: res
  character(len=1),dimension(natom) :: chain,nchain,atomtype
  integer,dimension(natom) :: n_atom,nres
  real,dimension(natom) :: O,B,x,y,z
  integer :: i
! leitura do dcd
  logical :: unitcell=.true.
  double precision :: sides(3)
  integer :: error,ioerror,stride=1,firstatom,lastatom,firstframe,lastframe
  character(len=4) :: dummyc
  integer :: dummyi,index,j,ntotframes,ntotat
  real :: dummyr
  double precision :: dummyd
! cálculo interno
  integer :: f,nframes,pepna,h,natom,pepnr,a,soma
  integer,dimension(pepnr,nframes) :: pwin,plin
  real :: dist
  character(len=4) :: pepc,waterc,lipc,dummy
  real,dimension(pepnr) :: mintwpres,mintlpres
  real,dimension(nframes) :: mintwpframe,mintlpframe

  i=0
  print*, "lendo pdb de referência..."
  open (unit=9, file=pdbfile, status="old")
!  read(unit=9, fmt='(A)'),  remark    ! lê a primeira linha de observação do pdb
  do i=1, natom, 1
    read(unit=9,fmt='(A4,1X,I6,1X,A4,1X,A3,1X,A1,1X,I3,4X,F8.3,F8.3,F8.3,2X,&
      &F4.2, 2X, F4.2, 6X, A4, A2,A2 )') ATOM(i),n_atom(i),nameatom(i),&
      res(i),chain(i),nres(i),x(i),y(i),z(i),O(i),B(i),namechain(i),nchain(i),atomtype(i)
  end do
  close (9)
  
  do f=1,nframes
    do i=1,pepnr
      pwin(i,f)=0
      plin(i,f)=0
    end do
  end do
  soma=0
  
!!loop para concatenar coordenadas na mesma linha na trajetória final
  print*, "lendo coordenadas..."
  open(unit=7,file="dcdfiles.dat")
  read(unit=7,fmt='(A)'),dcdfile
!!loop para ler vários dcds
  do while((dcdfile.ne."exit") .and. (soma<nframes))
!! checando número de frames e número de átomos da trajetória
    open(11,file=dcdfile,action='read',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      print*, "erro ao abrir o arquivo dcd"
      return
    end if
    read(11) dummyc, ntotframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(11) dummyi, dummyr
    read(11) ntotat
    print*, ntotframes,ntotat,dcdfile,pepnr,pepna
!! checado

    firstatom = 1 
    lastatom = natom
!    index = 0
    do f=1, ntotframes
      print*,f+soma
    
!      index = index + 1
!      j = 3*(index-1)
      if(unitcell) read(11) sides(1), dummyd, sides(2), dummyd, dummyd, sides(3)
!      iatom = natom*(index - 1)                                           
      read(11) (dummyr, j = 1, firstatom - 1 ), (x(j), j = 1, natom)
      read(11) (dummyr, j = 1, firstatom - 1 ), (y(j), j = 1, natom)            
      read(11) (dummyr, j = 1, firstatom - 1 ), (z(j), j = 1, natom)           
      
      do i=1,pepna
        dummy = nameatom(i)
        dummy = ADJUSTL(dummy)
        if (dummy(1:1) .ne. "H") then
          do h=pepna+1,natom
            dummy = nameatom(h)
            dummy = ADJUSTL(dummy)
            if (dummy(1:1) .ne. "H") then
              if ((x(i)-x(h))<4) then
                if ((y(i)-y(h))<4) then
                 if ((z(i)-z(h))<4) then
                    dist = sqrt( (x(i)-x(h))**2 + (y(i)-y(h))**2 + (z(i)-z(h))**2 )
                    if (dist <= 4.) then
                      if (namechain(h).eq.waterc) pwin(nres(i),soma+f)=pwin(nres(i),soma+f)+1
                      if (namechain(h).eq.lipc) plin(nres(i),soma+f)=plin(nres(i),soma+f)+1
                    end if
                  end if
                end if
              end if
            end if
          end do 
        end if
      end do
      print*, pwin(1,f),plin(1,f)
         
    end do
    soma=soma+ntotframes
    read(unit=7,fmt='(A)'),dcdfile
    close(11)
    
  end do
  soma=0
  close(7)
  
  do i=1,pepnr
    mintwpres(i)=0
    mintlpres(i)=0
    do f=1,nframes
      mintwpres(i)=mintwpres(i)+pwin(i,f)
      mintlpres(i)=mintlpres(i)+plin(i,f)
    end do
    mintwpres(i)=mintwpres(i)/nframes
    mintlpres(i)=mintlpres(i)/nframes
  end do
  
  do f=1,nframes
    mintwpframe(f)=0
    mintlpframe(f)=0
    do i=1,pepnr
      mintwpframe(f)=mintwpframe(f)+pwin(i,f)
      mintlpframe(f)=mintlpframe(f)+plin(i,f)
    end do
!    mintwpframe(f)=mintwpframe(f)/pepnr
!    mintlpframe(f)=mintlpframe(f)/pepnr
  end do
  
  open (unit=8,file=outfile)
  write(unit=8,fmt='(a)'), "Interactions per residue with water or lipid"
  do i=1,pepnr
    write(unit=8,fmt='(i4,2x,f6.2,2x,f6.2)'), i,mintwpres(i),mintlpres(i)
  end do
  write(unit=8,fmt='(a)'), "Interactions per frame with water or lipid"
  do f=1,nframes
    write(unit=8,fmt='(i5,2x,f6.2,2x,f6.2)'), f,mintwpframe(f),mintlpframe(f)
  end do
  close(8)
    
end subroutine

