module physics
    implicit none
    save
    real(kind = 8), allocatable, dimension(:, :) :: s, t, v, p, g, oldp, h, temp, & 
        f, configuration, a, d, x, xt, u_mh, e, fprime, cprime, c
    real(kind = 8), allocatable, dimension(:, :, :, :) :: v2
    real(kind = 8), allocatable, dimension(:) :: zetas, coef, expo, zs
    character(len = 2), allocatable :: atom_label(:)
    real(kind = 8) :: en = 0.0d0, ent = 0.0d0
end module physics

module constants
    implicit none
    save
    ! integer constants
    integer, parameter :: num_basis = 2, &  ! total number of basis (if num_basis .neq. num_atom, 
        ! the user must change the init()function to generate other basis functions)
        num_occupied_orb = 1, & ! number of occupied orbitals / 2 because of closed shell
        num_atom = 2, & ! number of atoms
        max_iter = 1e8 ! maximum of iterations
    ! real constants
    real(kind = 8), parameter :: pi = 3.1415926535898d0, crit = 1.0d-4 ! convergence criterion
end module constants

program main
    implicit none
    call hfcalc()
end program main

subroutine hfcalc()
    implicit none
    call parse_input_and_alloc_men("geometry_rep.xyz") ! change input file name here!
    call init()
    call init_mat()
    call scf()
end subroutine hfcalc

subroutine scf()
    use physics
    use constants
    implicit none
    integer :: i = 0, j = 0, k = 0, iter = 0
    real(kind = 8) :: delta = 0.0d0, distance

    do iter = 1, max_iter
    write(*, *) "========iteration number: ", iter, "======="
    
    write(*, *) "========================================================="
    write(*, *) "we output p here:"
    write(*, *) "========================================================="
    call matcat(p, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of P"
    write(*, *) "========================================================="

    call formg()
    write(*, *) "========================================================="
    write(*, *) "we output G here:"
    write(*, *) "========================================================="
    call matcat(G, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of G"
    write(*, *) "========================================================="

    do i = 1, num_basis
    do j = 1, num_basis
        F(i,j) = H(i,j) + G(i,j)
    end do
    end do


    EN = 0.0d0
    do i = 1, num_basis
    do j = 1, num_basis
        EN = EN+0.5D0*P(I,J)*(H(I,J)+F(I,J))
    end do
    end do

    ENT = 0.0d0
    do i = 1, num_basis
    do j = 1, num_basis
    if (i .gt. j) then
        ENT = EN + zs(i)*zs(j)/distance(i, j)
    end if
    end do
    end do

    write(*, *) "========================================================="
    write(*, *) "we output F here:"
    write(*, *) "========================================================="
    call matcat(F, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of F"
    write(*, *) "========================================================="
    write(*, *) "EN = ", EN

    CALL MULT(F,X,G,num_basis,num_basis)
    CALL MULT(XT,G,FPRIME,num_basis,num_basis)
    ! subroutine diag(fp, cp, e, n)
    CALL DIAG(FPRIME, CPRIME, E, num_basis)
    CALL MULT(X,CPRIME,C,num_basis,num_basis)

    write(*, *) "========================================================="
    write(*, *) "we output F' here:"
    write(*, *) "========================================================="
    call matcat(FPRIME, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of F'"
    write(*, *) "========================================================="
    write(*, *) "========================================================="
    write(*, *) "we output C' here:"
    write(*, *) "========================================================="
    call matcat(CPRIME, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of C'"
    write(*, *) "========================================================="
    write(*, *) "========================================================="
    write(*, *) "we output E here:"
    write(*, *) "========================================================="
    call matcat(E, num_basis)
    write(*, *) "========================================================="
    write(*, *) "End of E"
    write(*, *) "========================================================="

    ! form new P
    do i = 1, num_basis
    do j = 1, num_basis
    OLDP(I,J)=P(I,J)
    P(I,J)=0.0D0
    DO k = 1, num_occupied_orb
    P(I,J) = P(I,J) + 2.0d0*C(I,K)*C(J,K)
    end do
    end do
    end do

    DELTA=0.0D0
    do i = 1, num_basis
    do j = 1, num_basis
    DELTA=DELTA+(P(I,J)-OLDP(I,J))**2
    end do
    end do
    DELTA=DSQRT(DELTA/4.0D0)

    write(*, *) "ENT"
    write(*, *) ENT
    write(*, *) "delta"
    write(*, *) delta

    if (delta .lt. crit) then
        write(*, *) " Converged! "
        stop
    end if

    end do
end subroutine scf

subroutine mult(A,B,C,IM,M)
    implicit none
    real(kind = 8), intent(inout) :: A(IM,IM),B(IM,IM),C(IM,IM)
    integer :: i, j, k, IM, M
    c = 0.0d0
    do i = 1, M
    do j = 1, M
    do k = 1, M
        C(I,J)=C(I,J)+A(I,K)*B(K,J)
    end do
    end do
    end do
end subroutine mult

subroutine formg()
    use physics
    use constants
    implicit none
    integer :: i, j, k, l
    G = 0.0d0
    do i = 1, num_basis
    do j = 1, num_basis
    do k = 1, num_basis
    do l = 1, num_basis
        G(I,J)=G(I,J)+P(K,L)*(v2(I,J,K,L)-0.5D0*v2(I,L,K,J))
    end do
    end do
    end do
    end do
end subroutine formG

subroutine init_mat()
    use physics
    use constants
    implicit none
    integer :: i = 1, j = 1, k = 1, l = 1, m = 1, n = 1, &
        ii = 1, mm = 1, nn = 1
    real(kind = 8) S_f, distance, T_f, &
        c_coord(3), p_coord(3), q_coord(3), RCP, RCP2, v_f, TWOE, RPQ, RPQ2, F0
    write(*, *) "========================================================="
    write(*, *) "=====Initialize S, G, H .etc============================="
    write(*, *) "========================================================="
    do j = 1, num_basis
    do i = 1, num_basis
    do k = 1, 3
    do l = 1, 3
    s(i, j) = s(i, j) + S_f(a(i, k), a(j, l), distance(i, j)**2)*d(i, k)*d(j, l)
    !real(kind = 8) function T_f(A, B, RAB2)
    t(i, j) = t(i, j) + T_f(a(i, k), a(j, l), distance(i, j)**2)*d(i, k)*d(j, l)
    do m = 1, num_atom ! Calculate V
    !real(kind = 8) function V_f(A, B, RAB2, RCP2, ZC)
    c_coord = configuration(m, :)
    p_coord = (a(i, k)*configuration(i, :) + a(j, l)*configuration(j, :))/(a(i, k) + a(j, l))
    RCP2 = 0.0d0
    do n = 1, 3
    RCP2 = RCP2 + (c_coord(n) - p_coord(n))**2
    end do !n
    RCP = dsqrt(RCP2)
    v(i, j) = v(i, j) + v_f(a(i, k), a(j, l), distance(i, j)**2, RCP2, zs(m))*d(i, k)*d(j, l)
    end do !m
    end do
    end do
    end do
    end do

    do j = 1, num_basis
    do i = 1, num_basis
    do k = 1, 3
    do l = 1, 3
    do m = 1, num_atom ! Calculate V
    !real(kind = 8) function V_f(A, B, RAB2, RCP2, ZC)
    c_coord = configuration(m, :)
    p_coord = (a(i, k)*configuration(i, :) + a(j, l)*configuration(j, :))/(a(i, k) + a(j, l))
    RCP2 = 0.0d0
    do n = 1, 3
    RCP2 = RCP2 + (c_coord(n) - p_coord(n))**2
    end do !n
    RCP = dsqrt(RCP2)
    v(i, j) = v(i, j) + v_f(a(i, k), a(j, l), distance(i, j)**2, RCP2, zs(m))*d(i, k)*d(j, l)
    end do !m
    end do
    end do
    end do
    end do
    v = v / 2

    ! H matrix
    h = t + v

    ! two_components
    do i = 1, num_basis
    do j = 1, num_basis
    do k = 1, num_basis
    do l = 1, num_basis
    do m = 1, 3
    do n = 1, 3
    do mm = 1, 3
    do nn = 1, 3
    p_coord = (a(i, m)*configuration(i, :) + a(j, n)*configuration(j, :))/(a(i, m) + a(j, n))
    q_coord = (a(k, mm)*configuration(k, :) + a(l, nn)*configuration(l, :))/(a(k, mm) + a(l, nn))
    RPQ2 = 0.0d0
    do ii = 1, 3
    RPQ2 = RPQ2 + (p_coord(ii) - q_coord(ii))**2
    end do !ii
    RPQ = dsqrt(RPQ2)
    v2(i, j, k, l) = v2(i, j, k, l) + TWOE(a(i, m),a(j, n),a(k, mm),a(l, nn), distance(i, j)**2, distance(k, l)**2, RPQ2) &
        *d(i, m)*d(j, n)*d(k, mm)*d(l, nn)
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

    ! X matrix
    ! subroutine diag(fp, cp, e, n)
    call diag(s, x, e, num_basis)
    write(*, *) "======  X  matrix ======================================="
    call matcat(x, num_basis)
    write(*, *) "========End of X========================================="
    do i = 1, num_basis
    e(i, i) = 1 / dsqrt(e(i, i))
    end do
    !subroutine mult(A,B,C,IM,M)
    call mult(x, e, xt, num_basis, num_basis)
    x = xt
    xt = transpose(x)

    write(*, *) "======  S  matrix ======================================="
    call matcat(s, num_basis)
    write(*, *) "========End of S========================================="
    write(*, *) "======  T  matrix ======================================="
    call matcat(t, num_basis)
    write(*, *) "========End of T========================================="
    write(*, *) "======  V  matrix ======================================="
    call matcat(v, num_basis)
    write(*, *) "========End of V========================================="
    write(*, *) "======  H  matrix ======================================="
    call matcat(h, num_basis)
    write(*, *) "========End of H========================================="
    write(*, *) "======  X  matrix ======================================="
    call matcat(x, num_basis)
    write(*, *) "========End of X========================================="
    write(*, *) "======  two components =================================="
    do l = 1, num_basis
    do k = l, num_basis
    do j = 1, num_basis
    do i = j, num_basis
    write(*, 1040) "V2(", i, ", ", j, ", ", k, ", ", l, ") = ", v2(i, j, k, l)
    write(*, 1040) "V2(", j, ", ", i, ", ", k, ", ", l, ") = ", v2(j, i, k, l)
    write(*, 1040) "V2(", i, ", ", j, ", ", l, ", ", k, ") = ", v2(i, j, l, k)
    write(*, 1040) "V2(", j, ", ", i, ", ", l, ", ", k, ") = ", v2(j, i, l, k)
1040 format(A3, 3(I2, A2), I2, A4, F15.8)
    end do
    end do
    end do
    end do
    write(*, *) "======  End of two components ==========================="
    write(*, *) "========================================================="

end subroutine init_mat

subroutine diag(fp, cp, e, n)
    ! subroutine diag(fp, cp, e, n)
    implicit none
    ! 引入LAPACK函数
    external dsyev

    ! 定义参数
    integer :: i, n
    real(kind = 8), intent(inout) :: fp(n,n), e(1:n, 1:n), cp(n, n)
    real(kind = 8) :: w(n)
    integer :: info, lda, lwork
    real(kind = 8), allocatable :: work(:)
    character :: jobz, uplo

    ! 初始化矩阵

    lda = n
    lwork = 3*n-1
    allocate(work(lwork))

    ! 设置参数
    jobz = 'V'     ! 计算特征值和特征向量
    uplo = 'U'     ! 使用上三角矩阵

    ! 调用dsyev函数
    e = fp
    call dsyev(jobz, uplo, n, fp, lda, w, work, lwork, info)
    cp = fp
    fp = e

    e = 0.0d0
    do i = 1, n
    e(i, i) = w(i)
    end do

    deallocate(work)

end subroutine diag

subroutine matcat(mat, n)
    implicit none
    integer :: i, n
    real(kind = 8) mat(1:n, 1:n)
    do i = 1, n
    write(*, "(*(f8.3))") mat(i, 1:n)
    end do
end subroutine matcat

function distance(i, j)
    use physics
    implicit none
    integer i, j, k
    real(kind = 8) distance
    distance = 0.0d0
    do k = 1, 3
    distance = distance + (configuration(i, k) - configuration(j, k))**2
    end do
    distance = dsqrt(distance)
end function distance

real(kind = 8) function F0(ARG)
	!! DECLARACTIONS
    use constants
    implicit none
    real(kind = 8) :: ARG
    real :: DERFOTHER
	!! EXECUTABLES
    if (ARG .lt. 1.0d-6) then
        F0 = 1.0d0 - ARG/3.0d0
    else
        F0 = DSQRT(PI/ARG)*DERFOTHER(DSQRT(ARG))/2.0d0
    end if
    return
end function F0

real(kind = 8) function DERFOTHER(ARG)
    use constants
    implicit none
    integer :: ii
    real(kind = 8) :: ARG, A(5), PP, POLY, TN, T
    PP = 0.3275911D0
    A(1) = 0.254829592D0
    A(2) = -0.284496736D0
    A(3) = 1.421413741D0
    A(4) = -1.453152027D0
    A(5) = 1.061405429D0
    T=1.0D0/(1.0D0+PP*ARG)
    TN=T
    POLY=A(1)*TN
    do ii=2,5
    TN=TN*T
    POLY=POLY+A(ii)*TN
    end do
    DERFOTHER=1.0D0-POLY*dexp(-ARG*ARG)
    return
end function DERFOTHER

real(kind = 8) function S_f(A, B, RAB2)
    use constants
    implicit none
    real(kind = 8) :: A, B, RAB2
    S_f = (PI/(A+B))**1.5D0*DEXP(-A*B*RAB2/(A+B))
    return
end function

real(kind = 8) function T_f(A, B, RAB2)
    use constants
    implicit none
    real(kind = 8) :: A, B, RAB2
    T_f=A*B/(A+B)*(3.0D0-2.0D0*A*B*RAB2/(A+B))*(PI/(A+B))**1.5D0 & 
        *DEXP(-A*B*RAB2/(A+B))
    return
end function

real(kind = 8) function V_f(A, B, RAB2, RCP2, ZC)
    use constants
    implicit none
    real(kind = 8) :: A, B, RAB2, RCP2, ZC, F0
      V_f = 2.0D0*PI/(A+B)*F0((A+B)*RCP2)*DEXP(-A*B*RAB2/(A+B))
      V_f = -V_f*ZC
    return
end function

real(kind = 8) function TWOE(A,B,C,D,RAB2,RCD2,RPQ2)
    use constants
    implicit none
    real(kind = 8) :: A, B, C, D, RAB2, RCD2, RPQ2, F0
    TWOE = 2.0D0*(PI**2.5D0)/((A+B)*(C+D)*DSQRT(A+B+C+D)) &
    *F0((A+B)*(C+D)*RPQ2/(A+B+C+D)) & 
    *DEXP(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D))
    return
end function

subroutine init()
    use physics
    use constants
    implicit none
    integer :: i = 1, j = 1
    write(*, *) "========================================================="
    write(*, *) "=====Initialize zetas===================================="
    write(*, *) "========================================================="
    do i = 1, num_basis
    if (atom_label(i) .eq. "H") then
        zetas(i) = 1.24d0
        zs(i) = 1.0d0
    else
        zetas(i) = 2.0925d0
        zs(i) = 2.0d0
    end if
    write(*, 1010) " zeta(", i, ") = ", zetas(i)
1010 format((A6, I2, A4, F12.8))
    end do
    write(*, *) "========================================================="

    coef(1) = 0.444635D0
    coef(2) = 0.535328D0
    coef(3) = 0.154329D0
    expo(1) = 0.109818D0
    expo(2) = 0.405771D0
    expo(3) = 2.22766D0

    write(*, *) "========================================================="
    write(*, *) "=====Initialize a and d=================================="
    write(*, *) "========================================================="
    do j = 1, 3
    do i = 1, num_basis
    a(i, j) = expo(j)*(zetas(i)**2)
    d(i, j) = coef(j)*((2*a(i, j)/pi))**0.75
    write(*, 1020) "  a(", i, ', ', j, ') = ', a(i, j), "; d(", i, ', ', j, ') = ', d(i, j)
1020 format(2(A4, I2, A2, I2, A4, F12.8))
    end do
    end do
    write(*, *) "========================================================="
end subroutine init

subroutine parse_input_and_alloc_men(file_name)
    use physics
    use constants
    character(len = *) :: file_name
    ! *********************************************************************
    ! read configuration and store them in configuartion(1:n, 1:3)
    ! *********************************************************************
    allocate(configuration(1:num_basis, 1:3))
    allocate(atom_label(1:num_atom))
    write(*, *) "========================================================="
    write(*, *) "=======Reading Geomtry Configuration====================="
    write(*, *) "========================================================="
    do i = 1, num_basis
    open(unit = 10, file = file_name)
    read(10, '(A2,3F15.8)') atom_label(i), configuration(i, :)
    end do
    configuration = configuration * 1e-10 / 5.2918e-11
    do i = 1, num_basis
    write(*, 1030) "label = ", atom_label(i), "x = ", configuration(i, 1), "y = ", configuration(i, 2), & 
        "z = ", configuration(i, 3)
1030 format(A8, A2, 3(' ', A4,F12.8))
    end do
    write(*, *) "========================================================="
    close(unit = 10)
    ! *********************************************************************
    ! END: Initialization of configuation
    ! *********************************************************************

    !real(kind = 8), allocatable, dimension(:, :) :: s, t, v, p, h, temp& 
    !    f, configuration, a, d
    allocate(s(1:num_basis, 1:num_basis))
    allocate(t(1:num_basis, 1:num_basis))
    allocate(v(1:num_basis, 1:num_basis))
    allocate(p(1:num_basis, 1:num_basis))
    allocate(g(1:num_basis, 1:num_basis))
    allocate(e(1:num_basis, 1:num_basis))
    allocate(oldp(1:num_basis, 1:num_basis))
    allocate(x(1:num_basis, 1:num_basis))
    allocate(xt(1:num_basis, 1:num_basis))
    allocate(u_mh(1:num_basis, 1:num_basis))
    allocate(h(1:num_basis, 1:num_basis))
    allocate(temp(1:num_basis, 1:num_basis))
    allocate(f(1:num_basis, 1:num_basis))
    allocate(fprime(1:num_basis, 1:num_basis))
    allocate(cprime(1:num_basis, 1:num_basis))
    allocate(c(1:num_basis, 1:num_basis))
    allocate(a(1:num_basis, 1:3))
    allocate(d(1:num_basis, 1:3))
    allocate(zetas(1:num_basis))
    allocate(zs(1:num_basis))
    allocate(coef(1:3))
    allocate(expo(1:3))
    allocate(v2(1:num_basis,1:num_basis,1:num_basis,1:num_basis))
    s = 0.0d0
    t = 0.0d0
    v = 0.0d0
    p = 0.0d0
    g = 0.0d0
    h = 0.0d0
    temp = 0.0d0
    f = 0.0d0
    zs = 0.0d0
    v2 = 0.0d0
    oldp = 0.0d0
    x = 0.0d0
    xt = 0.0d0
    u_mh = 0.0d0
    e = 0.0d0
    fprime = 0.0d0
    cprime = 0.0d0
    c = 0.0d0

end subroutine parse_input_and_alloc_men
