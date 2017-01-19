        program whiten
        use sacio
        implicit none
        integer :: i, flag
        character(len=80) :: filename
        real,allocatable,dimension(:) :: data
        type(sachead) :: head
        real*4    seis_in(10000000),seis_out(10000000)
        real*4   seis_outamp(10000000), seis_outph(10000000)
        real*8   f1,f2,f3,f4,npow,dt,ns,dom
        filename="../IC.BJT.00.BHZ.norm"

        call sacio_readsac(filename, head, data, flag)
        
        dt = 1
        f1 = 0.01
        f2 = 0.02
        f3 = 0.067
        f4 = 0.08
        npow = 1
        do i=1,head%npts
            seis_in(i) = data(i)        
        enddo
        call fllter4(f1,f2,f3,f4,npow,dt,head%npts,seis_in,
     1  seis_out,seis_outamp,seis_outph,ns,dom)
        
        end program


