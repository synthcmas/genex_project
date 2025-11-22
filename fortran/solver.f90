module solve_diff
	use system_of_equations
	implicit none
contains
	subroutine solve(these_initials,these_finals,snapshots)
	    USE RKSUITE_90_PREC
   		USE RKSUITE_90
    	Use system_of_equations
		implicit none
		REAL(WP),DIMENSION(:,:,:),Intent(IN)::these_initials    
   	 	REAL(WP),DIMENSION(:,:,:,:),Intent(OUT)::these_finals
    	REAL(WP),DIMENSION(:),Intent(IN)::snapshots    
    	REAL(WP),DIMENSION(:,:,:),allocatable::DUMMY_VECTOR
    	LOGICAL,DIMENSION(:,:,:),allocatable::THIS_MASK    
	    Real(WP),dimension(:),allocatable::Y_GOT,YDERIV_GOT,THRESHOLDS,initials
    	Real(WP),dimension(:,:),allocatable::express_OUT
        REAL(WP)::TOLERANCE,this_time_end,time_int_success
   		TYPE(RK_COMM_REAL_1D) :: COMMB
    	integer::i,flag
   		allocate(initials(size(these_initials)))

        
	    initials=pack(these_initials,.True.)        



	    allocate(DUMMY_VECTOR,mold=these_initials)
	    allocate(THIS_MASK(size(these_initials,1),size(these_initials,2),size(these_initials,3)))    
    	! Storing the values after unpacking 
	    THIS_MASK=.True.


	    ALLOCATE(Y_GOT(size(initials)),YDERIV_GOT(size(initials)),THRESHOLDS(size(initials)))
    	ALLOCATE(express_OUT(size(snapshots),size(initials)))

    	express_OUT(1,:)=initials

    	TOLERANCE=1.0D-3
    	THRESHOLDS=1.0D-4

    	adaptive_expression:do i=1,size(snapshots)-1

        CALL SETUP(COMMB,snapshots(I),initials,snapshots(I+1),TOLERANCE,THRESHOLDS,method="H",MESSAGE=.TRUE.)
    
    	print*,'step',i,'done'         
        this_time_end=snapshots(I+1)
        
        CALL range_integrate(COMMB,system_eq,this_time_end,time_int_success,Y_GOT,YDERIV_GOT,FLAG=FLAG)
        IF(FLAG/=1) THEN
            PRINT*,"INTEGRATION NOT SUCCESSFUL"
            PRINT*,'COULD EVALUATE TILL day=',time_int_success,"INSTEAD N=",this_time_end
            IF(ABS(time_int_success-this_time_end)>1.0d-3) STOP "EXITING"
        ENDIF
        !==================== STORE INTERMEDIATE VALUES TILL NEND ==================
        express_OUT(I+1,:)=Y_GOT(:)    

        initials=Y_GOT

        CALL COLLECT_GARBAGE(COMMB)

    	end do adaptive_expression

      	do i=1,size(snapshots)
        	these_finals(i,:,:,:)=unpack(express_OUT(i,:),THIS_MASK,dummy_vector)
 	   	end do

	end subroutine solve	
end module solve_diff

