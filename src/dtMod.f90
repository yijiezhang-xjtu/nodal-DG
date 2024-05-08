module dtMod
    ! calculate scaling for dt in tri
    use constantsMod

    implicit none

    contains

    subroutine scaleDT(x,z,elem,nelem,sDT)
        !Input
        real(kind=CUSTOM_REAL), dimension(:) :: x,z
        integer, dimension(:,:) :: elem
        integer :: nelem
        !Output
        real(kind=CUSTOM_REAL), dimension(nelem) :: sDT
        !Local
        real(kind=CUSTOM_REAL) :: len1,len2,len3,sper,A
        integer :: ie

        do ie=1,nelem
           len1=sqrt( (x(elem(1,ie)) - x(elem(2,ie)))**2 + (z(elem(1,ie)) - z(elem(2,ie)))**2)
           len2=sqrt( (x(elem(2,ie)) - x(elem(3,ie)))**2 + (z(elem(2,ie)) - z(elem(3,ie)))**2)
           len3=sqrt( (x(elem(3,ie)) - x(elem(1,ie)))**2 + (z(elem(3,ie)) - z(elem(1,ie)))**2)
           sper = (len1+len2+len3)/2.0
           a= sqrt(sper*(sper-len1)*(sper-len2)*(sper-len3))
           sDT(ie) = a/sper
        end do
  end subroutine scaleDT
end module dtMod

