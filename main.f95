! A fortran95 program for G95
program main
    implicit none

    integer, parameter :: N_SIZE = 5, &
                          N_DIMENSION = N_SIZE
    real, parameter    :: vector_beta(4) = (/ -0.9, -1.1, -1.19, -1.199 /)

    integer            :: row, &
                          column, &
                          i, &
                          j, &
                          k, &
                          l, &
                          IPVT(N_SIZE)

    real               :: left_sum, &
                          right_sum, &
                          matrix(N_SIZE,N_SIZE)

    real               :: vector_a(N_SIZE), &
                          vector_b(N_SIZE), &
                          vector_g(N_SIZE), &
                          vector_z(N_SIZE), &
                          vector_n(N_SIZE), &
                          vector_x(N_SIZE) ! Solution vector

    real               :: WORK(N_SIZE), COND, CONDP1, max_g

do l = 1, size(vector_beta) ! Spin 4 time

    write (*, '(A, 1I1, A, F6.3, A)') "=====[", l, "]=====Рассчет при значении Beta: ", vector_beta(l), "======="
    left_sum = 0
    right_sum = 0
    COND = 0
    CONDP1 = 0
    max_g = 0

    do row = 1, N_SIZE
        k = row
        vector_a(row) = k ** 2
        vector_b(row) = k + vector_beta(l)
        vector_g(row) = abs(k - 3)
    end do


    ! Calculate the free vector
    do i = 1, N_SIZE
        if ( i == 1 ) then
            left_sum = 0
        else
            do k = 1, i - 1
                left_sum = left_sum + (vector_a(k) * vector_g(k))
            end do
            left_sum = vector_b(i) * left_sum
        end if

        vector_z(i) = left_sum

        left_sum = 0 ! reset
    end do


    do i = 1, N_SIZE
        do k = i, N_SIZE
            right_sum = right_sum + (vector_b(k) * vector_g(k))
        end do

        right_sum = vector_a(i) * right_sum
        vector_z(i) = vector_z(i) + right_sum

        right_sum = 0 ! reset
    end do


    ! Fill the matrix
    do column = 1, N_SIZE
        do row = 1, N_SIZE
            if ( column == 1 ) then
                matrix(column, row) = vector_a(column)
            else if (column == 2 .and. row >= 2) then
                matrix(column, row) = vector_a(column)
            else if (column == 3 .and. row >= 3) then
                matrix(column, row) = vector_a(column)
            else if (column == 4 .and. row >= 4) then
                matrix(column, row) = vector_a(column)
            else if (column == 5) then
                matrix(column, row) = vector_a(row)
            else
                matrix(column, row) = vector_a(row)
            end if
        end do
    end do


    ! Print the filled matrix (DEBUG)
    !do column = 1, N_SIZE
    !    do row = 1, N_SIZE
    !        print '(f8.4 $)', matrix(column, row)
    !    end do
    !    print "(/)"
    !end do


    ! Rewrite number An * Bn
    do column = 1, N_SIZE
        do row = 1, N_SIZE
            if ( column == 1 ) then
                matrix(column, row) = matrix(column, row) * vector_b(row)
            else if (column == 2 .and. row >= 2) then
                matrix(column, row) = matrix(column, row) * vector_b(row)
            else if (column == 3 .and. row >= 3) then
                matrix(column, row) = matrix(column, row) * vector_b(row)
            else if (column == 4 .and. row >= 4) then
                matrix(column, row) = matrix(column, row) * vector_b(row)
            else if (column == 5) then
                matrix(column, row) = matrix(column, row) * vector_b(column)
            else
                matrix(column, row) = matrix(column, row) * vector_b(column)
            end if
            !print '(f9.4 $)', matrix(column, row)
        end do
        !print "(/)"
    end do


    PRINT 101, ((matrix(i,j), j = 1, N_SIZE), vector_z(i), i = 1, N_SIZE)
    Call decomp(N_SIZE, N_DIMENSION, matrix, cond, ipvt, work)

    !print "(/)"
    PRINT 102, COND
    CONDP1 = COND + 1.0
    IF(CONDP1.EQ.COND) PRINT 103
    IF(CONDP1.EQ.COND) STOP

    Call SOLVE(N_SIZE, N_DIMENSION, matrix, vector_z, IPVT)

    vector_x = vector_z ! Fill solution vector

    ! MATRIX AFTER DECOMP (DEBUG)
    !do column = 1, N_SIZE
    !    do row = 1, N_SIZE
    !        print '(f8.4 $)', matrix(column, row)
    !    end do
    !    print "(/)"
    !end do

    PRINT 104, (vector_x(i), i = 1, N_SIZE)
    !STOP
    101 FORMAT(19X, 'Q', 28X, 'Z', 5(//1X, 5F7.2, 5X, F12.7))
    102 FORMAT(/'Число обусловленности матрицы Q (COND): ', E12.5)
    103 FORMAT(5X,'MATPИЦA KЛACCИФИЦИPУETCЯ KAK BЫPOЖДEHHAЯ')
    104 FORMAT(/5X,'BEKTOP PEШEHИЯ         X', 5(/18X, F15.7))


    max_g = abs( maxval(vector_g) )
    vector_n = abs( vector_x - vector_g )


    write (*, '(/A, F4.2)') "Норма MAX вектора G: ", max_g
    write (*, '(A, F9.7)') "Норма MAX вектора X - G: ", abs( maxval(vector_n) )
    write (*, '(A, F9.7/)') "Погрешность: ", abs( abs( maxval(vector_n) ) / max_g )
end do

end
