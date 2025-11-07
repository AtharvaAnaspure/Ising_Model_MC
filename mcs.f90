program mcs_fixed
    implicit none
    ! Declarations of integers and reals (double precision)
    integer :: i, j, k, L, a, b, c, d, ee, f
    integer :: niter, time, mm, nn, oo, N
    integer :: T_temp, n_equil, n_stat, Tf, Ti, dT, n_meas
    real*8 :: r, h, q, E, M, mag, Ei, Ef, dE, u,abs_mag
    real*8 :: av_m, av_e, av_m_N, av_e_N,abs_av_m_N, av_m2, av_e2, cv, chi, av_abs_m
    real*8 :: T, J_ising, av_m4, binder
    integer, dimension(:,:,:), allocatable :: spin

    ! Set the Ising coupling constant
    J_ising = 1.0d0

    ! Ask for simulation conditions
    print*, "Enter the number of lattice points per side (L):"
    read*, L
    print*, "Enter the Temperature range (Tf, Ti, dT) as integers (e.g., 400 500 1 for 4.00 to 5.00):"
    read*, Tf, Ti, dT
    print*, "Enter the number of Monte Carlo sweeps (niter):"
    read*, niter

    ! Set equilibration and measurement parameters
    n_equil = 10000
    n_stat  = 100

    ! Allocate the 3D spin lattice and initialize energy/magnetization counters
    allocate(spin(L,L,L))
    E = 0.0d0
    M = 0.0d0
    N = L * L * L

    call random_seed()

    ! Initialize the lattice with random spins and save configuration to a file
    open(10, file='MCS_lattice_3D.dat', status='replace')
    do i = 1, L
        do j = 1, L
            do k = 1, L
                call random_number(r)
                if (r < 0.5d0) then
                    spin(i,j,k) = -1
                else
                    spin(i,j,k) = 1
                end if
                write(10, '(I5,1X,I5,1X,I5,1X,I5)') i, j, k, spin(i,j,k)
            end do
        end do
    end do
    close(10)

    ! Compute initial total magnetization
    do i = 1, L
        do j = 1, L
            do k = 1, L
                M = M + spin(i,j,k)
            end do
        end do
    end do

    ! Compute initial energy using periodic boundary conditions
    do i = 1, L
        do j = 1, L
            do k = 1, L
                a = i + 1; if (i == L) a = 1
                b = i - 1; if (i == 1) b = L
                c = j + 1; if (j == L) c = 1
                d = j - 1; if (j == 1) d = L
                ee = k + 1; if (k == L) ee = 1
                f = k - 1; if (k == 1) f = L

                E = E - J_ising * real(spin(i,j,k) * ( &
                        spin(a,j,k) + spin(b,j,k) + spin(i,c,k) + &
                        spin(i,d,k) + spin(i,j,ee) + spin(i,j,f) ))
            end do
        end do
    end do
    E = E * 0.5d0  ! Correct for double counting

    mag = M / real(N)
    print*, 'Initial Energy E, E per spin =', E, E/real(N)
    print*, 'Initial Magnetization M, M per spin =', M, mag

    ! Open output file for temperature loop results
    open(11, file='T_loop_L32.dat', status='replace')

    ! Loop over the temperature range (T_temp in integer; convert to real T)
    do T_temp = Tf, Ti, dT
        T = dfloat(T_temp) / 100.0d0  ! E.g., 400 becomes 4.00

        ! Reset accumulators for the current temperature
        av_m   = 0.0d0
        av_e   = 0.0d0
        av_m_N = 0.0d0
        av_e_N = 0.0d0
        av_m2  = 0.0d0
        av_m4  = 0.0d0
        av_e2  = 0.0d0
        av_abs_m = 0.0d0
        ! Perform niter Monte Carlo sweeps.
        ! Each "sweep" consists of L^3 single-spin update attempts.
        do time = 1, niter
            do mm = 1, L
                do nn = 1, L
                    do oo = 1, L
                        ! Choose a random lattice site
                        call random_number(r)
                        i = int(r * real(L)) + 1
                        call random_number(q)
                        j = int(q * real(L)) + 1
                        call random_number(h)
                        k = int(h * real(L)) + 1

                        ! Determine neighbors with periodic boundary conditions
                        a = i + 1; if (i == L) a = 1
                        b = i - 1; if (i == 1) b = L
                        c = j + 1; if (j == L) c = 1
                        d = j - 1; if (j == 1) d = L
                        ee = k + 1; if (k == L) ee = 1
                        f = k - 1; if (k == 1) f = L

                        ! Calculate energy before the trial flip
                        Ei = -J_ising * real( spin(i,j,k) * ( &
                                spin(a,j,k) + spin(b,j,k) + spin(i,c,k) + &
                                spin(i,d,k) + spin(i,j,ee) + spin(i,j,f) ) )
                        ! Perform the trial flip
                        spin(i,j,k) = -spin(i,j,k)
                        ! Calculate energy after the flip
                        Ef = -J_ising * real( spin(i,j,k) * ( &
                                spin(a,j,k) + spin(b,j,k) + spin(i,c,k) + &
                                spin(i,d,k) + spin(i,j,ee) + spin(i,j,f) ) )
                        dE = Ef - Ei

                        ! Metropolis acceptance/rejection
                        if (dE <= 0.0d0) then
                            E = E + dE
                            M = M + (2.0d0 * real(spin(i,j,k)))
                        else
                            u = exp(-dE / T)
                            call random_number(h)  ! Reuse h for a new random number
                            if (h < u) then
                                E = E + dE
                                M = M + (2.0d0 * real(spin(i,j,k)))
                            else
                                ! Reject the trial flip: revert the spin
                                spin(i,j,k) = -spin(i,j,k)
                            end if
                        end if
                    end do
                end do
            end do

            ! After equilibration, accumulate data every n_stat sweeps
            if (time > n_equil) then
                if (mod(time, n_stat) == 0) then
                    mag = M / real(N)    ! Magnetization per spin
                    abs_mag = abs(M) / real(N)    ! Absolute value Magnetization per spin
                    av_m   = av_m   + mag
                    av_abs_m = av_abs_m + abs_mag
                    av_e   = av_e   + E / real(N)
                    av_m_N = av_m_N + M
                    !abs_av_m_N = av_m_N + abs(M)
                    av_e_N = av_e_N + E
                    av_m2  = av_m2  + (M * M)
                    av_m4  = av_m4  + (M * M * M * M)
                    av_e2  = av_e2  + (E * E)
                end if
            end if
        end do  ! End of Monte Carlo sweeps

        ! Calculate the number of measurement points
        n_meas = (niter - n_equil) / n_stat
        if (n_meas <= 0) then
            print*, "Error: No measurement points collected. Increase niter or decrease n_equil."
            stop
        end if

        ! Compute the averages (per measurement point)
        av_m = av_m / dfloat(n_meas)
        av_abs_m = av_abs_m / dfloat(n_meas)
        av_e = av_e / dfloat(n_meas)

        ! Compute the variance-based estimates (specific heat and susceptibility)
        cv  = ((av_e2 / dfloat(n_meas)) - ( (av_e_N / dfloat(n_meas))**2 )) / (T * T)
        chi = ((av_m2 / dfloat(n_meas)) - ( (av_m_N / dfloat(n_meas))**2 )) / T
        !binder's cumulant
        binder = 1.0d0 - (av_m4 / dfloat(n_meas)) / (3.0d0 * ( (av_m2 / dfloat(n_meas))**2 ))
        ! Write the results for the current temperature:
        ! T, average magnetization per spin, average energy per spin, specific heat, susceptibility
        write(11,*) T, av_m, av_e, cv, chi,av_abs_m,binder
    end do  ! End temperature loop

    close(11)
    print*, "Final values (for last temperature):"
    print*, "Magnetization per spin =", av_m
    print*, "Energy per spin        =", av_e
    print*, "Specific heat          =", cv
    print*, "Susceptibility         =", chi
    print*, "Absolute value Magnetization per spin =", av_abs_m
    print*, "Binder's cumulant      =", binder
    deallocate(spin)
end program mcs_fixed
