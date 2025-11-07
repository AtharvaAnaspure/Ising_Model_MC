Program ising_model
implicit none
integer :: i,j,k,L,p,a,b,c,d,niter,time,mm,nn,N
real :: r,h,q,E,M,mag,Ei,Ef,dE,o,u
real :: T=2.0, J_ising=1.0 ! Assigning values to relevant parameters
integer, dimension(:,:),allocatable :: spin
integer :: seed 
character(len=30):: charac_a,charac_b
charac_b = 'store_config'

print*,'Enter number of lattice points in 1D'
read*,L
print*,'Enter number of iterations'
read*,niter

allocate(spin(L,L))
E = 0.0 ! Instantaneous Energy of Lattice
M = 0.0 ! Instantaneous Magnetization of Lattice
N = L*L ! Total number of spins in lattice

call random_seed

! Initialize the lattice
open(10,file='initial_ising.dat')
p=0
do i=1,L
    do j=1,L  
        call random_number(r)
        spin(j,i) = 1
        ! Writing initial configuration
        write(10,'(I5,1X,I5,1X,I5,1X,I5)') i, j, p, spin(j,i)
    end do
end do
close(10)

! Calculate initial energy and magnetization
do i=1,L
    do j=1,L 
        a = i+1
        b = i-1
        c = j+1
        d = j-1
        if(i==L)a=1
        if(i==1)b=L
        if(j==1)d=L
        if(j==L)c=1
        M = M + spin(i,j)
        E = E - J_ising * real(spin(i,j) * (spin(a,j) + spin(b,j) + spin(i,c) + spin(i,d)))
    end do
end do

mag = M / real(N)
E = E * 0.5
print*,'Initial Energy E, E per spin =',E,E/real(N)
print*,'Initial Magnetization M, M per spin =',M,mag

! Initialization Complete
! Evolve it to reach equilibrium
open(11,file='ising_t2_N40_init_random.dat')
do time=1,niter
    do mm = 1,L
        do nn = 1,L
            call random_number(r); i = int(r * real(L)) + 1 ! Choosing a lattice site 
            call random_number(q); j = int(q * real(L)) + 1 
            a = i+1
            b = i-1
            c = j+1
            d = j-1 ! Neighbours of the spin(i,j)
            if(i==L)a=1; if(i==1)b=L; if(j==1)d=L; if(j==L)c=1 ! PBC

            Ei = -J_ising * real(spin(i,j) * (spin(a,j) + spin(b,j) + spin(i,c) + spin(i,d))) ! Before trial flip
            spin(i,j) = -spin(i,j) ! Trial flip
            Ef = -J_ising * real(spin(i,j) * (spin(a,j) + spin(b,j) + spin(i,c) + spin(i,d))) ! After trial flip
            dE = Ef - Ei ! Difference in energies
            if(dE <= 0.0) then
                E = E + dE ! Updating Energy and Magnetisation of lattice
                M = M + (2.0 * real(spin(i,j)))
            else
                u = exp(-dE / T)
                call random_number(h)
                if(h < u) then
                    E = E + dE
                    M = M + (2.0 * real(spin(i,j))) ! Instantaneous Magnetization of the whole lattice
                else
                    spin(i,j) = -spin(i,j) ! Reverting the trial flip
                end if
            end if
        end do
    end do
    write(11,*) time, M / real(N), E / real(N)   
end do
close(11)

End Program ising_model