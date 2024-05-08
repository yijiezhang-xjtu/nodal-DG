module fileunitMod
    use constantsMod

    implicit none

    contains

    function getFileUnit(fu_max)
        ! get_file_unit returns a unit number that is not in use
        integer :: getfileunit
        integer fu, fu_max, m, iostat
        logical opened

        m = fu_max ; if (m < 1) m = 97
        do fu = m,1,-1
            inquire (unit=fu, opened=opened, iostat=iostat)
            if (iostat.ne.0) cycle
            if (.not.opened) exit
        end do
        getFileUnit = fu
    end function getFileUnit
end module fileunitMod
