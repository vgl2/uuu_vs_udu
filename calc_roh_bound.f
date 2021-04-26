        subroutine calc_roh_bound(n,psips,key,r_oh,r_coord,nh,nox)
        implicit real*8(a-h,o-z)
        parameter (nmax=100000)
        parameter (ndim=27)
        dimension psips(ndim,nmax),sec_roh1(nmax),sec_roh2(nmax),
     1  sec_roh3(nmax),sec_roh4(nmax),sec_roh_1(nmax),sec_roh_2(nmax),
     1  r_oh(nmax),r_coord(3,nmax),nh(nmax),nox(nmax)
C		Classifies in a water molecule, which hydrogen is hydrogen bonding.
C		Inputs:
C		n = number of walkers
C		psips = Cartesian coordinates of walkers in units of bohr
C		key = identification of which water to look for hydrogen bound OH.
C		Outputs
C		r_oh =  Hydrogen bound rOH distance in bohr
C		r_coord = Cartesian coordinates of the rOH distance calculated in bohr
C		nh = classification of hydrogen that is hydrogen bound (either 1 or 2)
C		nox = classification of oxygen that the hydrogen is bound to (either 1,2 or 3)
        do i = 1,n
            if (key.eq.1) then
                sec_roh1(i) = sqrt((psips(4,i)-psips(10,i))**2+
     1          (psips(5,i)-psips(11,i))**2+(psips(6,i)-psips(12,i))**2)
                sec_roh2(i) = sqrt((psips(4,i)-psips(19,i))**2+(
     1          psips(5,i)-psips(20,i))**2+(psips(6,i)-
     1          psips(21,i))**2)
                sec_roh3(i) = sqrt((psips(7,i)-psips(10,i))**2+(
     1          psips(8,i)-psips(11,i))**2+(psips(9,i)-
     1          psips(12,i))**2)
                sec_roh4(i) = sqrt((psips(7,i)-psips(19,i))**2+(
     1          psips(8,i)-psips(20,i))**2+(psips(9,i)-
     1          psips(21,i))**2)
            else if (key.eq.2) then
                sec_roh1(i) = sqrt((psips(13,i)-psips(1,i))**2+(
     1          psips(14,i)-psips(2,i))**2+(psips(15,i)-
     1          psips(3,i))**2)
                sec_roh2(i) = sqrt((psips(13,i)-psips(19,i))**2+(
     1          psips(14,i)-psips(20,i))**2+(psips(15,i)-
     1          psips(21,i))**2)
                sec_roh3(i) = sqrt((psips(16,i)-psips(1,i))**2+(
     1          psips(17,i)-psips(2,i))**2+(psips(18,i)-
     1          psips(3,i))**2)
                sec_roh4(i) = sqrt((psips(16,i)-psips(19,i))**2+(
     1          psips(17,i)-psips(20,i))**2+(psips(18,i)-
     1          psips(21,i))**2)
            else if (key.eq.3) then
                sec_roh1(i) = sqrt((psips(22,i)-psips(1,i))**2+(
     1          psips(23,i)-psips(2,i))**2+(psips(24,i)-
     1          psips(3,i))**2)
                sec_roh2(i) = sqrt((psips(22,i)-psips(10,i))**2+(
     1          psips(23,i)-psips(11,i))**2+(psips(24,i)-
     1          psips(12,i))**2)
                sec_roh3(i) = sqrt((psips(25,i)-psips(1,i))**2+(
     1          psips(26,i)-psips(2,i))**2+(psips(27,i)-
     1          psips(3,i))**2)
                sec_roh4(i) = sqrt((psips(25,i)-psips(10,i))**2+(
     1          psips(26,i)-psips(11,i))**2+(psips(27,i)-
     1          psips(12,i))**2)
            endif
        enddo
        do i = 1,n
            if (sec_roh1(i).lt.sec_roh2(i)) then
                sec_roh_1(i) = sec_roh1(i)
            else
                sec_roh_1(i) = sec_roh2(i)
            endif
            if (sec_roh3(i).lt.sec_roh4(i)) then
                sec_roh_2(i) = sec_roh3(i)
            else
                sec_roh_2(i) = sec_roh4(i)
            endif
            if (key.eq.1) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(4,i)-psips(1,i))**2+(
     1              psips(5,i)-psips(2,i))**2+(psips(6,i)-
     1              psips(3,i))**2)
                    r_coord(1,i) = psips(4,i)-psips(1,i)
                    r_coord(2,i) = psips(5,i)-psips(2,i)
                    r_coord(3,i) = psips(6,i)-psips(3,i)
                    nh(i) = 1
                    if (sec_roh1(i).lt.sec_roh2(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                else
                    r_oh(i) = sqrt((psips(7,i)-psips(1,i))**2+(
     1              psips(8,i)-psips(2,i))**2+(psips(9,i)-
     1              psips(3,i))**2)
                    r_coord(1,i) = psips(7,i)-psips(1,i)
                    r_coord(2,i) = psips(8,i)-psips(2,i)
                    r_coord(3,i) = psips(9,i)-psips(3,i)
                    nh(i) = 2
                    if (sec_roh3(i).lt.sec_roh4(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                endif
            else if (key.eq.2) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(13,i)-psips(10,i))**2+
     1              (psips(14,i)-psips(11,i))**2+(psips(15,i)-
     1              psips(12,i))**2)
                    r_coord(1,i) = psips(13,i)-psips(10,i)
                    r_coord(2,i) = psips(14,i)-psips(11,i)
                    r_coord(3,i) = psips(15,i)-psips(12,i)
                    nh(i) = 1
                    if (sec_roh1(i).lt.sec_roh2(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                else
                    r_oh(i) = sqrt((psips(16,i)-psips(10,i))**2+
     1              (psips(17,i)-psips(11,i))**2+(psips(18,i)-
     1              psips(12,i))**2)
                    r_coord(1,i) = psips(16,i)-psips(10,i)
                    r_coord(2,i) = psips(17,i)-psips(11,i)
                    r_coord(3,i) = psips(18,i)-psips(12,i)
                    nh(i) = 2
                    if (sec_roh3(i).lt.sec_roh4(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                endif
            else if (key.eq.3) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(22,i)-psips(19,i))**2+
     1              (psips(23,i)-psips(20,i))**2+(psips(24,i)-
     1              psips(21,i))**2)
                    r_coord(1,i) = psips(22,i)-psips(19,i)
                    r_coord(2,i) = psips(23,i)-psips(20,i)
                    r_coord(3,i) = psips(24,i)-psips(21,i)
                    nh(i) = 1
                    if (sec_roh1(i).lt.sec_roh2(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                else
                    r_oh(i) = sqrt((psips(25,i)-psips(19,i))**2+
     1              (psips(26,i)-psips(20,i))**2+(psips(27,i)-
     1              psips(21,i))**2)
                    r_coord(1,i) = psips(25,i)-psips(19,i)
                    r_coord(2,i) = psips(26,i)-psips(20,i)
                    r_coord(3,i) = psips(27,i)-psips(21,i)
                    nh(i) = 2
                    if (sec_roh3(i).lt.sec_roh4(i)) then
                        nox(i) = 1
                    else
                        nox(i) = 2
                    endif
                endif
            endif
        enddo
        return
        end subroutine


        subroutine calc_roh_free(n,psips,key,r_oh,r_coord,nh)
        implicit real*8(a-h,o-z)
        parameter (nmax=100000)
        parameter (ndim=27)
        dimension psips(ndim,nmax),sec_roh1(nmax),sec_roh2(nmax),
     1  sec_roh3(nmax),sec_roh4(nmax),sec_roh_1(nmax),sec_roh_2(nmax),
     1  r_oh(nmax),r_coord(3,nmax),nh(nmax),test(4,nmax)
C		Classifies in a water molecule, which hydrogen is part of the free OH.
C		Inputs:
C		n = number of walkers
C		psips = Cartesian coordinates of walkers in units of bohr
C		key = identification of which water to look for free OH.
C		Outputs
C		r_oh =  free rOH distance in bohr
C		r_coord = Cartesian coordinates of the rOH distance calculated in bohr
C		nh = classification of hydrogen that is hydrogen bound (either 1 or 2)
        do i = 1,n
            if (key.eq.1) then
                sec_roh1(i) = sqrt((psips(4,i)-psips(10,i))**2+
     1          (psips(5,i)-psips(11,i))**2+(psips(6,i)-psips(12,i))**2)
                sec_roh2(i) = sqrt((psips(4,i)-psips(19,i))**2+(
     1          psips(5,i)-psips(20,i))**2+(psips(6,i)-
     1          psips(21,i))**2)
                sec_roh3(i) = sqrt((psips(7,i)-psips(10,i))**2+(
     1          psips(8,i)-psips(11,i))**2+(psips(9,i)-
     1          psips(12,i))**2)
                sec_roh4(i) = sqrt((psips(7,i)-psips(19,i))**2+(
     1          psips(8,i)-psips(20,i))**2+(psips(9,i)-
     1          psips(21,i))**2)
            else if (key.eq.2) then
                sec_roh1(i) = sqrt((psips(13,i)-psips(1,i))**2+(
     1          psips(14,i)-psips(2,i))**2+(psips(15,i)-
     1          psips(3,i))**2)
                sec_roh2(i) = sqrt((psips(13,i)-psips(19,i))**2+(
     1          psips(14,i)-psips(20,i))**2+(psips(15,i)-
     1          psips(21,i))**2)
                sec_roh3(i) = sqrt((psips(16,i)-psips(1,i))**2+(
     1          psips(17,i)-psips(2,i))**2+(psips(18,i)-
     1          psips(3,i))**2)
                sec_roh4(i) = sqrt((psips(16,i)-psips(19,i))**2+(
     1          psips(17,i)-psips(20,i))**2+(psips(18,i)-
     1          psips(21,i))**2)
            else if (key.eq.3) then
                sec_roh1(i) = sqrt((psips(22,i)-psips(1,i))**2+(
     1          psips(23,i)-psips(2,i))**2+(psips(24,i)-
     1          psips(3,i))**2)
                sec_roh2(i) = sqrt((psips(22,i)-psips(10,i))**2+(
     1          psips(23,i)-psips(11,i))**2+(psips(24,i)-
     1          psips(12,i))**2)
                sec_roh3(i) = sqrt((psips(25,i)-psips(1,i))**2+(
     1          psips(26,i)-psips(2,i))**2+(psips(27,i)-
     1          psips(3,i))**2)
                sec_roh4(i) = sqrt((psips(25,i)-psips(10,i))**2+(
     1          psips(26,i)-psips(11,i))**2+(psips(27,i)-
     1          psips(12,i))**2)
            endif
        enddo
        do i = 1,n
            if (sec_roh1(i).lt.sec_roh2(i)) then
                sec_roh_1(i) = sec_roh1(i)
            else
                sec_roh_1(i) = sec_roh2(i)
            endif
            if (sec_roh3(i).lt.sec_roh4(i)) then
                sec_roh_2(i) = sec_roh3(i)
            else
                sec_roh_2(i) = sec_roh4(i)
            endif
            if (key.eq.1) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(7,i)-psips(1,i))**2+(
     1              psips(8,i)-psips(2,i))**2+(psips(9,i)-
     1              psips(3,i))**2)
                    r_coord(1,i) = psips(7,i)-psips(1,i)
                    r_coord(2,i) = psips(8,i)-psips(2,i)
                    r_coord(3,i) = psips(9,i)-psips(3,i)
                    nh(i) = 2
                else
                    r_oh(i) = sqrt((psips(5,i)-psips(1,i))**2+(
     1              psips(5,i)-psips(2,i))**2+(psips(6,i)-
     1              psips(3,i))**2)
                    r_coord(1,i) = psips(4,i)-psips(1,i)
                    r_coord(2,i) = psips(5,i)-psips(2,i)
                    r_coord(3,i) = psips(6,i)-psips(3,i)
                    nh(i) = 1
                endif
            else if (key.eq.2) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(16,i)-psips(10,i))**2+
     1              (psips(17,i)-psips(11,i))**2+(psips(18,i)-
     1              psips(12,i))**2)
                    r_coord(1,i) = psips(16,i)-psips(10,i)
                    r_coord(2,i) = psips(17,i)-psips(11,i)
                    r_coord(3,i) = psips(18,i)-psips(12,i)
                    nh(i) = 2
                else
                    r_oh(i) = sqrt((psips(13,i)-psips(10,i))**2+
     1              (psips(14,i)-psips(11,i))**2+(psips(15,i)-
     1              psips(12,i))**2)
                    r_coord(1,i) = psips(13,i)-psips(10,i)
                    r_coord(2,i) = psips(14,i)-psips(11,i)
                    r_coord(3,i) = psips(15,i)-psips(12,i)
                    nh(i) = 1
                endif
            else if (key.eq.3) then
                if (sec_roh_1(i).lt.sec_roh_2(i)) then
                    r_oh(i) = sqrt((psips(25,i)-psips(19,i))**2+
     1              (psips(26,i)-psips(20,i))**2+(psips(27,i)-
     1              psips(21,i))**2)
                    r_coord(1,i) = psips(25,i)-psips(19,i)
                    r_coord(2,i) = psips(26,i)-psips(20,i)
                    r_coord(3,i) = psips(27,i)-psips(21,i)
                    nh(i) = 2
                else
                    r_oh(i) = sqrt((psips(22,i)-psips(19,i))**2+
     1              (psips(23,i)-psips(20,i))**2+(psips(24,i)-
     1              psips(21,i))**2)
                    r_coord(1,i) = psips(22,i)-psips(19,i)
                    r_coord(2,i) = psips(23,i)-psips(20,i)
                    r_coord(3,i) = psips(24,i)-psips(21,i)
                    nh(i) = 1
                endif
            endif
        enddo
        return
        end subroutine

        subroutine calc_dot_product(n,vec1,vec2,dot)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        dimension vec1(3,nmax),vec2(3,nmax),dot(nmax)
C		Calculates the dot product between 2 vectors
C		Inputs:
C		n = number of walkers
C		vec1 = coordinates of the first vector
C		vec2 = coordinates of the second vector
C		Outputs:
C		dot = dot product of the 2 vectors
        do i = 1,n
            dot(i) = 0.d0
        enddo
        do i = 1,n
            do j = 1,3
                dot(i) = dot(i) + vec1(j,i)*vec2(j,i)
            enddo
        enddo
        return
        end subroutine
