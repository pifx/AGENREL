program sargo_inb
! compute inbreeding as in sargolzaei et al. 2005
! works from a reordered and renumbered pedigree
! computes inbreeding coefficient for all animals
implicit none
integer, parameter:: unped=15
integer:: n, i, k, maxgen, npar, c, nb_anim
integer, allocatable:: ped(:,:), psigen(:), ispar(:), bgen(:)
integer, allocatable:: xrf(:), ped2(:,:), isex(:), ispar2(:)
real:: tot_inb, t1, t2, t11, t22
real, allocatable:: dd(:), inb(:), A22(:,:), dd2(:)
character*40:: fped

print *, 'name of pedigree?'
read (*,*) fped
call cpu_time(t1)
fped=trim(adjustl(fped))
open(unped, file=fped)
n=count_lines(unped)
allocate(ped(n,3),psigen(n),dd(n),inb(n),ispar(n),isex(n))
print *, '  reading pedigree and computing pseudo-generations ...'
do i=1,n
        read(unped,*) ped(i,1), ped(i,2), ped(i,3)
enddo
print *, 'x'
psigen=pseudo_gen(ped,n)
print *, '             ... done'
maxgen=maxval(psigen)+1
print '(a,i9)', ' # of animals in pedigree:', n
print '(a,i3)', ' # of pseudo-generations (including founders):', maxgen
allocate(bgen(maxgen))
bgen=bounds_gen(psigen,n,maxgen)
dd=1.d0; inb=0.d0
print *, ' working on first generations ...'
! animals from generation 0 : inb = 0
! animals from generation 1 : inb = 0, but dd can be <1
do i=1, n
        if (psigen(i)==0) cycle
        if (psigen(i)>1) exit
        if (ped(i,2)>0.and.ped(i,3)>0) dd(i)=0.5
        if (ped(i,2)>0.and.ped(i,3)==0) dd(i)=0.75
        if (ped(i,2)==0.and.ped(i,3)>0) dd(i)=0.75
enddo
print *, '             ... done'
print *, '---------------------------------------------------------'
print '(a)', ' PSEUDO-GEN   #ANIMALS   #PARENTS AVG_INBREEDING TIME[MIN]'
! animals from generation 2 to max: inb can be /= 0, dd can be <1
do k=3, maxgen
        ! pseudogeneration (=k-1)/ nb parents / nb offspring/ time /
        ! avg(inb)
        ! find the parents of generation k-1
        ! (those parents are in generation 0 (k=1)) to k-2
        call cpu_time(t11)
        allocate(xrf(bgen(k-1)))
        ispar=0; tot_inb=0.d0; nb_anim=0
        do i=1, n
                if (psigen(i)>k-1) exit
                if (psigen(i)<k-1) cycle
                if (ped(i,2)>0) ispar(ped(i,2))=1
                if (ped(i,3)>0) ispar(ped(i,3))=1
        enddo
        npar=sum(ispar)
        ! get cross-references for the parents
        c=1
        do i=1, bgen(k-1)
                if (ispar(i)==1) then
                        xrf(i)=c
                        c=c+1
                endif
        enddo
        ! extract the pedigree
        isex=ispar
        call extraction1(ped,isex,n)
        nex=sum(isex)
        allocate(ped2(nex,3),ispar2(nex))
        call extraction2(ped,isex,ispar,n,ped2,ispar2,nex)
        ! compute A22 between parents
        allocate(dd2(nex),A22(npar,npar))
        dd2=pack(dd,isex==1)
        call colleau3(ped2,ispar2,dd2,nex,A22,npar)
        ! get the inbreeding and dd for animals of generation k-1
        do i=1, n
                if (psigen(i)<k-1) cycle
                if (psigen(i)>k-1) exit
                if (ped(i,2)>0.and.ped(i,3)>0) then
                        inb(i)=.5*A22(xrf(ped(i,2)),xrf(ped(i,3)))
                        dd(i)=.5-.25*(inb(ped(i,2))+inb(ped(i,3)))
                elseif (ped(i,2)>0.and.ped(i,3)==0) then
                        dd(i)=.75-.25*inb(ped(i,2))
                elseif (ped(i,2)==0.and.ped(i,3)>0) then
                        dd(i)=.75-.25*inb(ped(i,3))
                endif
                tot_inb=tot_inb+inb(i)
                nb_anim=nb_anim+1
        enddo
        deallocate(ispar2,xrf,ped2,A22,dd2)
        call cpu_time(t22)
        print '(3i11,f12.4,f12.2)', k-1, nb_anim,  npar, tot_inb/real(nb_anim), (t22-t11)/60
enddo
print *, '----------------------------------------------------------'
print *, ' saving results ...'
! save results
open(30,file=fped(:len_trim(fped))//'.inb')
do i=1, n
        write(30,'(3i12,2f9.4,i4)') ped(i,1), ped(i,2), ped(i,3), inb(i), dd(i), psigen(i)
enddo
print *, '             ... done'
call cpu_time(t2)
print '(a,f10.2,a)', ' computations achieved in ', (t2-t1)/60,' min.'
print *, ' file *.inb contains pedigree, inbreeding, D-values & pseudo-generations'
close(30)

contains


function count_lines(un_file) result(n)
implicit none
integer:: un_file, n, k
character*1:: junk
rewind(un_file); n=0
do
        read(un_file,*,iostat=k) junk
        if (k/=0) exit
        n=n+1
enddo
rewind(un_file)
end function count_lines

function pseudo_gen(ped,n) result(psi)
implicit none
integer:: n, ped(n,3), psi(n), psip, psid
do i=1, n
        if (ped(i,2)==0.and.ped(i,3)==0) then
                psi(i)=0
        else
                psip=0; psid=0
                if (ped(i,2)>0) psip=psi(ped(i,2))
                if (ped(i,3)>0) psid=psi(ped(i,3))
                psi(i)=maxval((/psip, psid/))+1
        endif
enddo
end function


function bounds_gen(psigen,n,m) result(bgen)
! m = maxgen
implicit none
integer:: n, m, psigen(n), bgen(m), i
do i=2,n
        if (psigen(i-1)/=psigen(i)) bgen(psigen(i-1)+1)=i-1
enddo
bgen(m)=n
end function bounds_gen

subroutine extraction1(ped,isex,n)
implicit none
integer:: i, n, ped(n,3), isex(n)
do i=n,1,-1
        if (isex(i)==0) cycle
        if (ped(i,2)>0) isex(ped(i,2))=1
        if (ped(i,3)>0) isex(ped(i,3))=1
enddo
end subroutine extraction1

subroutine extraction2(ped,isex,ispar,n,ped2,ispar2,nex)
implicit none
integer:: n, nex ped(n,3), isex(n), ispar(n), ped2(nex,3), ispar2(nex)
integer:: i, c
integer, allocatable:: xrf(:)
allocate(xrf(n))
ispar2=0; ped2=0; xrf=0; c=0
do i=1,n
        if (isex(i)==0) cycle
        c=c+1
        xrf(i)=c
        ped2(c)=c
        if (ped(i,2)>0) ped2(c,2)=xrf(ped(i,2))
        if (ped(i,3)>0) ped2(c,3)=xrf(ped(i,3))
        if (ispar(i)==1) ispar2(c)=1
enddo
end subroutine extraction2

subroutine colleau3(ped,isgen,d,n,A22,ng)
implicit none
integer:: i, k, n, ng, c, isgen(n), ped(n,3)
integer, allocatable:: inds(:)
real:: vs, vd, x, d(n), A22(ng,ng)
real, allocatable:: q(:), r(:), u(:), v(:)
allocate(inds(ng), u(n), q(n), r(n), v(n))
c=1
do i=1,n
        if (isgen(i)==1) then
                inds(c)=i
                c=c+1
        endif
enddo
! computation of A22
do k=1,ng
        q=0; v=0; r=0
        r(inds(k))=1
        do i=n, 1, -1
                q(i)=q(i)+r(i)
                if (ped(i,2)>0) q(ped(i,2))=q(ped(i,2))+q(i)/2
                if (ped(i,3)>0) q(ped(i,3))=q(ped(i,3))+q(i)/2
        enddo
        do i=1,n
                vs=0.0;vd=0.0
                if (ped(i,2)>0) vs=v(ped(i,2))
                if (ped(i,3)>0) vd=v(ped(i,3))
                v(i)=q(i)*d(i)+(vs+vd)/2
        enddo
        A22(:,k)=v(inds)
enddo
!A22=v(inds,:)
deallocate(v)

end subroutine colleau3


end program
