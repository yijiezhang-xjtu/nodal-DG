module progressbarMod

    implicit none

    contains

    subroutine progress(j,total)
        !input
        integer :: j
        integer :: total
        !local
        integer :: k
        character(len=80) :: bar

        bar="|  Progress: [                                                         ] ???%  |"

        write(unit=bar(74:76),fmt="(i3)") 100*j/total

        do k = 1, j*57/total
            bar(14+k:14+k)="#"
        enddo
        ! print the progress bar.
        write(unit=6,fmt="(a1,a80)",advance="no") char(13), bar
        if (j /= total) then
            flush(unit=6)
        else
            write(unit=6, fmt=*)
        endif
        return
    end subroutine progress

end module progressbarMod
