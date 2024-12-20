module globals
  implicit none

  integer kminC,kmaxC,kminI,kmaxI
  integer N_node,N_nodeC,N_nodeI
  integer Itot,CStot
  integer tmax

  real(8) lambdaC,lambdaI
  real(8) kmedC,kmedI  
  real(8) r  !random number
  real(8) ff
  real(8) eps
  real(8) beta

  integer,allocatable,dimension(:)::edge,edgep
  integer,allocatable,dimension(:)::node,nodep
  integer,allocatable,dimension(:)::kk,kkp
  integer,allocatable,dimension(:)::State
  integer,allocatable,dimension(:)::Svec,Ivec
  
  real(8),allocatable,dimension(:)::PkC,PkI
    
end module globals

module random
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)
   real(8)::nmax=2*(2**30-1)+1
   real(8) sig
end module random


Program MainProgram
  use globals
  use random

  implicit none

  integer i,j,rea,nrea
  integer intento
  integer nptosf
  real(8) product
  
  character(8) Var01
  character(3) Var02
  character(8) Var03
  character(8) Var04


  call initialize_random
  
  !Parameters

  N_nodeI=10000000

  nrea   =100

  eps    =0.000001d0
  

  kminC  =2!0
  kmaxC  =20!10
  lambdaC=7d0
  
  kminI  =1!3
  kmaxI  =20!3 
  lambdaI=3d0
  
  beta   =1d0

  ff     =0.3d0
  
  nptosf =1

  
  write(Var01,'(I8)') N_nodeI
  call zeros(Var01)
  write(Var03,'(F8.6)') ff
  call zeros(Var03)
  write(Var04,'(F8.6)') beta
  call zeros(Var04)

  allocate(PkC(kminC:kmaxC),PkI(kminI:kmaxI))
  allocate(Svec(0:N_nodeI),Ivec(0:N_nodeI))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Degree distribution Cliques
  PkC           =0d0
  if(kminC==0)PkC(0)        =exp(-lambdaC)

  product      =1d0
  do i=1,kminC-1
     product=product*i
  enddo  
  do i=kminC,kmaxC
     if(i==0)cycle
     product=product*i         
     PkC(i)=1.d0*exp(-lambdaC)*lambdaC**(i)/dble(product)
     !PkC(i)=(1d0*i)**(-lambdaC)
  enddo  
  !PkC  =1d0
  
  PkC  =PkC/sum(PkC)
  
  kmedC=0d0
  do i=kminC,kmaxC
     kmedC=kmedC+i*PkC(i)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !Degree distribution Individuals
  PkI           =0d0
  if(kminI==0)PkI(0)        =exp(-lambdaI)

  product      =1d0
  do i=1,kminI-1
     product=product*i
  enddo  
  do i=kminI,kmaxI
     if(i==0)cycle

     product=product*i         
     PkI(i)=1.d0*exp(-lambdaI)*lambdaI**(i)/dble(product)
     !PkI(i)=(1d0*i)**(-lambdaI)
  enddo  
  PkI  =PkI/sum(PkI)
     
  kmedI=0d0
  do i=kminI,kmaxI
     kmedI=kmedI+i*PkI(i)
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  N_nodeC=kmedI*(N_nodeI/kmedC)
  N_node=N_nodeC+N_nodeI
  print*,PkC
  print*,PkI  
  print*,"Number of cliques", N_nodeC,"Number of individuals",N_nodeI, "N_nodeC+N_nodeI",N_node


  allocate(kk(N_node),node(N_node),State(N_nodeI),kkp(N_nodeI),nodep(N_nodeI))  
  

  rea=0
  intento=0
  do while(rea<100)     
     intento=intento+1
     call ConfModel     
           
     
     call SIQ
     if(1d0-Svec(tmax)/dble(N_nodeI)>0.05d0)then
        if(mod(rea,1)==0)print*,rea
        rea=rea+1
        write(Var02,'(I3)') rea
        call zeros(Var02)
        open(2,file='evolInf_NI_'//Var01//'_f_'//Var03//'_bet_'//Var04//"_"//Var02//'.dat')
        write(2,*)'#',rea,intento
        write(2,*)'#','kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI'
        write(2,*)'#',kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI
        do i=1,tmax
           write(2,*) i,Ivec(i)/dble(N_nodeI)
        enddo
        close(2)
        
        open(2,file='evolSus_NI_'//Var01//'_f_'//Var03//'_bet_'//Var04//"_"//Var02//'.dat')
        write(2,*)'#',rea,intento
        write(2,*)'#','kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI'
        write(2,*)'#',kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI
        do i=1,tmax
           write(2,*) i,Svec(i)/dble(N_nodeI)
        enddo
        close(2)   
             
     endif
     deallocate(edge,edgep)          
  enddo
  
  

end Program MainProgram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   SIQ
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SIQ
  use globals
  use random
  implicit none
  integer i,j
  integer sel,nod,vec
  integer cont
  integer Nindex,Nleaders,Nremanent,Aux
  integer Ninf,NinfAux,NisoInew
  integer tt,Ntot,Nsus
  integer, dimension(N_nodeI):: ListIndInf,ListIndIsolInf
  integer, dimension(N_nodeI):: ListIndInfNew,ListIndIsolInfNew
  integer, dimension(N_nodeI):: Leader
  integer, dimension(N_nodeI):: Received
  integer, dimension(N_nodeI):: ListSel

  
  !State=0  :susceptible
  !State=1  :infected
  !State=2  :quarantined
  
  Leader        =0
  
  Ivec          =0
  Svec          =N_nodeI
  
  State         =0
  Ninf          =0   !number of infected individuals

  ListIndInf    =0   !List of infected individuals
  ListIndIsolInf=0   !List of isolated individuals


  Itot          =0   !number of individuals ever infected
  tt            =0

  ListSel       =0
  do i=1,N_nodeI
     ListSel(i)=i
  enddo

  Nleaders=ff*N_nodeI
  Nremanent=N_nodeI

  do i=1,Nleaders
     call rand
     sel =r*Nremanent+1
     Aux =ListSel(sel)
     
     Leader(Aux) =1
     ListSel(sel)=ListSel(Nremanent)
     Nremanent=Nremanent-1
  enddo

  !INDEX-CASES
  
  Nindex        =Nremanent*eps
  cont          =0
  do i=1,Nindex
     call rand
     sel =r*Nremanent+1
     Aux =ListSel(sel)
     State(Aux)=1
     Ninf =Ninf+1     
     ListIndInf(Ninf) =Aux     

     ListSel(sel)=ListSel(Nremanent)
     Nremanent=Nremanent-1
  enddo
  
  Itot      =Ninf

  Received  =0
  do while(Ninf>0)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !New time step
     tt    =tt+1  
     !STEP 1: TRANSMISSION
     NinfAux       =0
     ListIndInfNew =0

     do i=1,Ninf
        nod =ListIndInf(i)                             
        do j=1,kkp(nod)
           vec =edgep(nodep(nod)+j-1)
           call rand
           if(r>beta)cycle
           if(State(vec)==0.and.Received(vec)==0)then
              Received(vec)  =1
              NinfAux     =NinfAux+1              
              ListIndInfNew(NinfAux)=vec              
           endif
        enddo      
     enddo
     Itot            =Itot+NinfAux
       
     !STEP 2a: QUARANTINE (people with access to testint)
     NisoInew         =0
     ListIndIsolInfNew=0
     do i=1,NinfAux
        nod  =ListIndInfNew(i)
        if(Leader(nod)==1)then
           NisoInew =NisoInew+1
           ListIndIsolInfNew(NisoInew)=nod
           State(nod)=2
        endif
     enddo
     
     !STEP 2b: QUARANTINE (NEIGHBORS)

     do i=1,NisoInew
        nod =ListIndIsolInfNew(i)
        do j=1,kkp(nod)
           vec =edgep(nodep(nod)+j-1)
           State(vec) =2           
        enddo
     enddo

     
     !STEP 3: UPDATE
     !Those individuals who have received the infection and have not been isolated, will now 
     !move to the infected compartment
     
     Ninf       =0
     ListIndInf =0


     do i=1,NinfAux
        nod =ListIndInfNew(i)
        if(State(nod)==0.and.Received(nod)==1)then
           State(nod)      =1
           Ninf            =Ninf+1
           ListIndInf(Ninf)=nod
        endif        
     enddo
     
     Nsus=0
     do i=1,N_nodeI
        if(State(i)==0)Nsus=Nsus+1
     enddo
     
     
     Ivec(tt)=NinfAux
     Svec(tt)=Nsus
     print*,NinfAux/dble(N_nodeI),Nsus/dble(N_nodeI),tt
     !read(*,*)
     tmax =tt

  enddo 
  

end subroutine SIQ




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Configuration Model
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine creates three arrays: kk, edge, and node, which contain all the network topology information.
!kk(i) is the degree of node "i".
!The subset edge(node(i):node(i)+kk(i)-1) contains the list of neighbors of node "i".

subroutine ConfModel
  use globals
  use random
  implicit none
  integer i,j,k,ll
  integer nod1,lug,nod2,vec
  integer Aux,flag
  integer m,n,stubmaxpos,stubminpos
  integer counting, nstubs,nstubsC,nstubsI
  integer nstubAux,stubposC,stubposI

  real(8) ws

  integer, allocatable, dimension(:)::listStubC,listStubI,kkAux

  real(8),dimension(kminC-1:kmaxC)::CumulativePkC
  real(8),dimension(kminI-1:kmaxI)::CumulativePkI
  
  
  allocate(kkAux(n_node))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution group 1 (cliques)
  
  CumulativePkC   =0d0
  ws             =0d0
  do i=kminC,kmaxC
     ws                = ws+PkC(i)
     CumulativePkC(i)   = ws
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution group 2 (Individuals)
  
  CumulativePkI   =0d0
  ws             =0d0
  do i=kminI,kmaxI
     ws                = ws+PkI(i)
     CumulativePkI(i)   = ws
  enddo  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Assigning Connectivity to all nodes
2 kk=0 
  do i=1,N_nodeC
     call rand
     do j=kminC,kmaxC
        if(CumulativePkC(j-1)<r.and.r<=CumulativePkC(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo
  do i=N_nodeC+1,N_node
     call rand
     do j=kminI,kmaxI
        if(CumulativePkI(j-1)<r.and.r<=CumulativePkI(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo  

  nstubsC =sum(kk(1:N_nodeC))
  nstubsI =sum(kk(N_nodeC+1:N_node))
  nstubs  =nstubsC+nstubsI
  
  
  if(abs(nstubsC-nstubsI)>0.01d0*(kmedI*N_nodeI)) goto 2
  
  allocate(listStubC(nstubsC),listStubI(nstubsI))  
  allocate(edge(nstubs))

  listStubC     =0
  listStubI     =0
  nstubAux      =0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!List of stubs

  do i=1,N_nodeC
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStubC(nstubAux)=i
     enddo
  enddo
  
  nstubAux      =0
  do i=N_nodeC+1,N_node
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStubI(nstubAux)=i
     enddo
  enddo  

  node(1)       =1
  do i=1,n_node-1
     node(i+1)  = node(i)+kk(i)
  enddo

  edge        =0
  kkAux       =0
  counting    =0
  do while(nstubsC>0.and.nstubsI>0)
3    if(counting==100) then
	deallocate(listStubC,listStubI,edge)
	go to 2
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Randomly choose two stubs
     
     call rand
     stubposC    =r*nstubsC+1
     call rand
     stubposI    =r*nstubsI+1

     m        =listStubC(stubposC)   !the first stub corresponds to node "m"
     n        =listStubI(stubposI)   !the second stub corresponds to node "n"

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Checking whether these two stubs can be joined or not

     if(m.ne.n) then
	do k=1,kkAux(m)
           if(edge(node(m)+k-1)==n)then
              counting     =counting+1
              go to 3   !if nodes "m" and "n" are already connected, we have to choose another pair of stubs
           endif
	enddo
	edge(node(m)+kkAux(m))   =n
	edge(node(n)+kkAux(n))   =m

	kkAux(m)                 =kkAux(m)+1
	kkAux(n)                 =kkAux(n)+1
	listStubC(stubposC)      =listStubC(nstubsC)
	listStubI(stubposI)      =listStubI(nstubsI)
	
	nstubsC             =nstubsC-1
	nstubsI             =nstubsI-1
	counting            =0
     else
        counting           =counting+1
        go to 3 !self-loops are forbidden
     endif
  enddo
  
  kk   =kkAux
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!One-mode Projection
  
  do i=N_nodeC+1,N_node
     nod1 =i-N_nodeC
     Aux  =0
     do j=1,kk(i)
        lug      =edge(node(i)+j-1)
        Aux      =Aux+kk(lug)
     enddo
     kkp(nod1)  =Aux-kk(i)
  enddo    
  
  nodep    =0
  nodep(1) =1

  do i=2,N_nodeI
     nodep(i)=nodep(i-1)+kkp(i-1)
  enddo
  
    
  allocate(edgep(sum(kkp)))
  
  kkp    =0
  edgep  =0

  do i=N_nodeC+1,N_node
     nod1 =i-N_nodeC
     
     do j=1,kk(i)
        lug      =edge(node(i)+j-1)
        do k=1,kk(lug)
           nod2  =edge(node(lug)+k-1)-N_nodeC
           if(nod2==nod1)cycle
           flag  =0
           do ll=1,kkp(nod1)
              vec   =edgep(nodep(nod1)+ll-1)
              if(vec==nod2)then
                 flag =1
                 exit
              endif
           enddo
           if(flag==0)then
              edgep(nodep(nod1)+kkp(nod1))   =nod2
              edgep(nodep(nod2)+kkp(nod2))   =nod1
              kkp(nod1)  =kkp(nod1)+1
              kkp(nod2)  =kkp(nod2)+1
           endif
        enddo
     enddo  
  enddo    
  
  
  deallocate(kkAux,listStubC,listStubI)
end subroutine ConfModel



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Seed for the Pseudorandom number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)
  !hour=1300826913
  print*,'hour',hour
  see=hour
  CALL RANDOM_SEED(put=see)		
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)   =r*nmax
     ir(i)    =i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Pseudo-random number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rand
  use globals
  use random
  implicit none

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Zero Character
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zeros(b)
 implicit none
 
 character(*):: b
 integer i
 
 do i=1,len(b)
   if (b(i:i)==' ') then  
      b(i:i)='0'
   endif
 end do
end subroutine zeros


