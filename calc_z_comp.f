        program calc_z
        implicit real*8(a-h,o-z)
        parameter (nmax=100000)
        parameter (nwavetot=100)
        parameter (ndim=27)
        parameter (natoms=ndim/3)
        parameter (nmaxtot=nmax*nwavetot)
        parameter (nbin=200)
        parameter (nmono=3)
        dimension n(nwavetot),n0(nwavetot),simtime(nwavetot),
     1  psips(ndim,nmax),x(natoms,3,nmax),o_center(3,nmax),dot1(nmax),
     1  roox1(nmax),rooy1(nmax),rooz1(nmax),roox2(nmax),rooy2(nmax),
     1  rooz2(nmax),roox3(nmax),rooy3(nmax),rooz3(nmax),cross(3,nmax),
     1  cross_norm(nmax),r_coord1(3,nmax),r_coord2(3,nmax),dot2(nmax),
     1  r_coord3(3,nmax),check1(nmax),check2(nmax),check3(nmax),
     1  uud(nmax),wt(nmax,nwavetot),nwt(nwavetot),tot_wt(nwavetot),
     1  tot_uud(nwavetot),nh1(nmax),nh2(nmax),nh3(nmax),roh1(nmax),
     1  roh2(nmax),roh3(nmax),dot3(nmax),rod1(nmax),rod2(nmax),
     1  r_trans(3,nmax),n2(nwavetot),n01(nwavetot),simtime1(nwavetot),
     1  psips_noopt(ndim,nmax),psips_check(ndim,nmax),rsec(nmax),
     1  psips_check2(ndim,nmax),rsec2(nmax),nhcheck(nmax),rsec3(nmax),
     1  rsec4(nmax),roh3_bound(nmax),r_coord3_bound(3,nmax),nox3(nmax),
     1  nh3_bound(nmax),udu(nmax),duu(nmax),ddu(nmax),dud(nmax),
     1  udd(nmax),tot_udu(nmax),tot_duu(nmax),tot_ddu(nmax),uuu(nmax),
     1  tot_dud(nmax),tot_udd(nmax),ddd(nmax),tot_uuu(nmax),
     1  tot_ddd(nmax)
C		Calculates the orientation of the z-displacement of the 
C		outer hydrogen in water trimer in order to classify the position
C		of the outer hydrogen whether it is above the plane (U) or below
C		the plane (D)
C		Inputs:
C		coordinates of walkers after performing geometry optimization.
C		descendant weights of the walkers
C		Outputs:
C		percentages of the free hydrogens in the relative orentations of 
C		whether or not the free hydrogen is above (U) or below the plane. 
        open(unit=8,file='wf-all-d-opt.dat',status='old')
        open(unit=9,file='tot_weights_all-d.dat',status='old')
        pi = dacos(-1.d0)
        avg_uud = 0.d0
        dev_uud = 0.d0
        avg_udu = 0.d0
        dev_udu = 0.d0
        avg_ddu = 0.d0
        dev_ddu = 0.d0
        avg_duu = 0.d0
        dev_duu = 0.d0
        avg_udd = 0.d0
        dev_udd = 0.d0
        avg_dud = 0.d0
        dev_dud = 0.d0
        avg_uuu = 0.d0
        dev_ddd = 0.d0
        read(8,*) nwave
        do k = 1,nwave
            tot_wt(k) = 0.d0
            read(9,*) nwt(k)
            do i = 1,nwt(k)
                read(9,*) wt(i,k)
                tot_wt(k) = tot_wt(k) + wt(i,k)
            enddo
        enddo
        ip = 0
        do k = 1,nwave
            icount = 0
            read(8,*) n(k),n0(k),simtime(k)
            do i = 1,n(k)
                read(8,*) (psips(j,i),j=1,ndim)
            enddo
            do i = 1,n(k)
                roox1(i) = psips(10,i)-psips(1,i)
                rooy1(i) = psips(11,i)-psips(2,i)
                rooz1(i) = psips(12,i)-psips(3,i)
                roox2(i) = psips(19,i)-psips(10,i)
                rooy2(i) = psips(20,i)-psips(11,i)
                rooz2(i) = psips(21,i)-psips(12,i)
                roox3(i) = psips(19,i)-psips(1,i)
                rooy3(i) = psips(20,i)-psips(2,i)
                rooz3(i) = psips(21,i)-psips(3,i)
            enddo
            do i = 1,n(k)
                cross_norm(i) = 0.d0
            enddo
            do i = 1,n(k)
                cross(1,i) = rooy1(i)*rooz2(i)-(rooz1(i)*rooy2(i))
                cross(2,i)=-(roox1(i)*rooz2(i)-(rooz1(i)*roox2(i)))
                cross(3,i) = roox1(i)*rooy2(i)-(rooy1(i)*roox2(i))
            enddo
            do i = 1,n(k)
                do j = 1,3
                    cross_norm(i) = cross_norm(i)+(cross(j,i)**2)
                enddo
            enddo
            do i = 1,n(k)
                do j = 1,3
                    cross(j,i) = cross(j,i)/sqrt(cross_norm(i))
                enddo
            enddo
            call calc_roh_free(n(k),psips(:,:),1,roh1(:),r_coord1(:,:),
     1      nh1)
            call calc_roh_free(n(k),psips(:,:),2,roh2(:),r_coord2(:,:),
     1      nh2)
            call calc_roh_free(n(k),psips(:,:),3,roh3(:),r_coord3(:,:),
     1      nh3)
            call calc_roh_bound(n(k),psips(:,:),3,roh3_bound(:),
     1      r_coord3_bound(:,:),nh3_bound(:),nox3(:))
            do i = 1,n(k)
                if (nox3(i).ne.1) then
                    r_trans(:,i) = r_coord3(:,i)
                    r_coord3(:,i) = r_coord2(:,i)
                    r_coord2(:,i) = r_trans(:,i)
                endif
            enddo
c           need to check if the cross product is actually perpendicular
c           calculate secondary roh to determine which is h-bonded
            call calc_dot_product(n(k),cross(:,:),r_coord1(:,:),dot1(:))
            call calc_dot_product(n(k),cross(:,:),r_coord2(:,:),dot2(:))
            call calc_dot_product(n(k),cross(:,:),r_coord3(:,:),dot3(:))
            do i = 1,n(k)
                if (dot1(i).gt.0.d0) then
                    check1(i) =1.d0
                else
                    check1(i) = -1.d0
                endif
                if (dot2(i).gt.0.d0) then
                    check2(i) =1.d0
                else
                    check2(i) = -1.d0
                endif
                if (dot3(i).gt.0.d0) then
                    check3(i) =1.d0
                else
                    check3(i) = -1.d0
                endif
            enddo
            do i = 1,n(k)
                if ((check1(i).eq.1.d0).and.(check2(i).eq.1.d0).and.
     1          (check3(i).eq.-1.d0)) then
                    uud(i) = 1.d0
                else
                    uud(i) = 0.d0
                endif
                if ((check1(i).eq.1.d0).and.(check2(i).eq.-1.d0).and.
     1          (check3(i).eq.1.d0)) then
                    udu(i) = 1.d0
                else
                    udu(i) = 0.d0
                endif
                if ((check1(i).eq.-1.d0).and.(check2(i).eq.1.d0).and.
     1          (check3(i).eq.1.d0)) then
                    duu(i) = 1.d0
                else
                    duu(i) = 0.d0
                endif
                if ((check1(i).eq.-1.d0).and.(check2(i).eq.-1.d0).and.
     1          (check3(i).eq.1.d0)) then
                    ddu(i) = 1.d0
                else
                    ddu(i) = 0.d0
                endif
                if ((check1(i).eq.-1.d0).and.(check2(i).eq.1.d0).and.
     1          (check3(i).eq.-1.d0)) then
                    dud(i) = 1.d0
                else
                    dud(i) = 0.d0
                endif
                if ((check1(i).eq.1.d0).and.(check2(i).eq.-1.d0).and.
     1          (check3(i).eq.-1.d0)) then
                    udd(i) = 1.d0
                else
                    udd(i) = 0.d0
                endif
                if ((check1(i).eq.1.d0).and.(check2(i).eq.1.d0).and.
     1          (check3(i).eq.1.d0)) then
                    uuu(i) = 1.d0
                else
                    uuu(i) = 0.d0
                endif
                if ((check1(i).eq.-1.d0).and.(check2(i).eq.-1.d0).and.
     1          (check3(i).eq.-1.d0)) then
                    ddd(i) = 1.d0
                else
                    ddd(i) = 0.d0
                endif
            enddo
            tot_uud(k) = 0.d0
            tot_udu(k) = 0.d0
            tot_duu(k) = 0.d0
            tot_ddu(k) = 0.d0
            tot_dud(k) = 0.d0
            tot_udd(k) = 0.d0
            tot_uuu(k) = 0.d0
            tot_ddd(k) = 0.d0
            do i = 1,n(k)
                tot_uud(k) = tot_uud(k)+ wt(i,k)*uud(i)
                tot_udu(k) = tot_udu(k)+ wt(i,k)*udu(i)
                tot_duu(k) = tot_duu(k)+ wt(i,k)*duu(i)
                tot_ddu(k) = tot_ddu(k)+ wt(i,k)*ddu(i)
                tot_dud(k) = tot_dud(k)+ wt(i,k)*dud(i)
                tot_udd(k) = tot_udd(k)+ wt(i,k)*udd(i)
                tot_uuu(k) = tot_uuu(k)+ wt(i,k)*uuu(i)
                tot_ddd(k) = tot_ddd(k)+ wt(i,k)*ddd(i)
            enddo
            tot_uud(k) = tot_uud(k)/tot_wt(k)
            avg_uud = avg_uud + tot_uud(k)
            tot_udu(k) = tot_udu(k)/tot_wt(k)
            avg_udu = avg_udu + tot_udu(k)
            tot_duu(k) = tot_duu(k)/tot_wt(k)
            avg_duu = avg_duu + tot_duu(k)
            tot_ddu(k) = tot_ddu(k)/tot_wt(k)
            avg_ddu = avg_ddu + tot_ddu(k)
            tot_dud(k) = tot_dud(k)/tot_wt(k)
            avg_dud = avg_dud + tot_dud(k)
            tot_udd(k) = tot_udd(k)/tot_wt(k)
            avg_udd = avg_udd + tot_udd(k)
            tot_uuu(k) = tot_uuu(k)/tot_wt(k)
            avg_uuu = avg_uuu + tot_uuu(k)
            tot_ddd(k) = tot_ddd(k)/tot_wt(k)
            avg_ddd = avg_ddd + tot_ddd(k)
        enddo
        avg_uud = avg_uud/dfloat(nwave)
        avg_udu = avg_udu/dfloat(nwave)
        avg_duu = avg_duu/dfloat(nwave)
        avg_ddu = avg_ddu/dfloat(nwave)
        avg_dud = avg_dud/dfloat(nwave)
        avg_udd = avg_udd/dfloat(nwave)
        avg_uuu = avg_uuu/dfloat(nwave)
        avg_ddd = avg_ddd/dfloat(nwave)
        do k = 1,nwave
            dev_uud = dev_uud +(tot_uud(k)-avg_uud)**2
            dev_udu = dev_udu +(tot_udu(k)-avg_udu)**2
            dev_duu = dev_duu +(tot_duu(k)-avg_duu)**2
            dev_ddu = dev_ddu +(tot_ddu(k)-avg_ddu)**2
            dev_dud = dev_dud +(tot_dud(k)-avg_dud)**2
            dev_udd = dev_udd +(tot_udd(k)-avg_udd)**2
            dev_uuu = dev_uuu +(tot_uuu(k)-avg_uuu)**2
            dev_ddd = dev_ddd +(tot_ddd(k)-avg_ddd)**2
        enddo
        dev_uud = sqrt(dev_uud/dfloat(nwave))
        dev_udu = sqrt(dev_udu/dfloat(nwave))
        dev_duu = sqrt(dev_duu/dfloat(nwave))
        dev_ddu = sqrt(dev_ddu/dfloat(nwave))
        dev_dud = sqrt(dev_dud/dfloat(nwave))
        dev_udd = sqrt(dev_udd/dfloat(nwave))
        dev_uuu = sqrt(dev_uuu/dfloat(nwave))
        dev_udd = sqrt(dev_udd/dfloat(nwave))
        print *, 'uud,udu,duu,ddu,dud,udd,uuu,ddd'
        print *, avg_uud,dev_uud,avg_udu,dev_udu,avg_duu,dev_duu,
     1  avg_ddu,dev_ddu,avg_dud,dev_dud,avg_udd,dev_udd,avg_uuu,dev_uuu,
     1  avg_ddd,dev_ddd
        end program
           
                
            
