subroutine gr_engine(n_part, n_bins, r_bin, z1, z2, Lx, Ly, r0, g_r)
implicit none
integer, intent(in) :: n_part, n_bins
real(kind=8), intent(in) :: r0(3,n_part), r_bin, z1, z2, Lx, Ly
real(kind=8), intent(out) :: g_r(n_bins)
integer i_part, j_part, i_dim, bin_indx, i_bin
real(kind=8) :: dist, delta_r(2), boundary(2), r_bin2, pi
!real(kind==8), parameter :: 
pi = 4 * atan (1.0)

r_bin2 = r_bin**2
boundary(1) = Lx
boundary(2) = Ly

do i_part = 1, n_part - 1
    if ( ( r0(3,i_part) > z2 ) .or. ( r0(3,i_part) < z1 ) ) cycle
    do j_part = i_part + 1, n_part
        if ( ( r0(3,j_part) > z2 ) .or. ( r0(3,j_part) < z1 ) ) cycle
        dist = 0.
        do i_dim = 1, 2
            delta_r(i_dim) = r0(i_dim,j_part) - r0(i_dim,i_part)
            delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim) * int( delta_r(i_dim) * 2 / boundary(i_dim) )
            dist = dist + delta_r(i_dim)**2
        end do

        dist = sqrt(dist)
        bin_indx = dist / r_bin
        g_r(bin_indx) = g_r(bin_indx) + 1.
    end do
end do

do i_bin = 1, n_bins
    g_r(i_bin) = g_r(i_bin) / ( pi * r_bin2 * ( 1. + 2. * i_bin ) )
end do

end subroutine



subroutine gr_test(n_part, n_bins, r0, g_r)
implicit none
integer, intent(in) :: n_part, n_bins
real(kind=8), intent(in) :: r0(3,n_part)
real(kind=8), intent(out) :: g_r(n_bins)
integer i_part, j_part, i_dim, bin_indx, i_bin
real(kind=8) :: dist, delta_r(2), boundary(2), r_bin2, pi

do i_bin = 1, n_bins
    g_r(i_bin) = float(i_bin)
end do


end subroutine


