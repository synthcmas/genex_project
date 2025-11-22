program main
	USE parameters
	USE RKSUITE_90_PREC
	use solve_diff
 	USE INIFILE    
 	use utils 	
	implicit none
	real(wp),allocatable,dimension(:,:,:)::init_values	
	real(wp),allocatable,dimension(:,:,:,:)::final_values
	real(wp),allocatable,dimension(:)::timesteps
	integer::i,m,n,writeUnit
	character(LEN=128)::foldername='plot_data',write_filename
	character(LEN=128)::slice_num
	character(LEN=Ini_max_string_len)::InputFile
	logical::bad


	InputFile = ''
	if (iargc() /= 0)  call getarg(1,InputFile)
 	if (InputFile == '') stop 'No parameter input file'

 	call Ini_Open(InputFile, 1, bad, .false.)
  	if (bad) stop 'Error opening parameter file'
 	Ini_fail_on_not_found = .false.


    num_mRNA=Ini_Read_Int('Num_mRNA',100)
    num_Protein=Ini_Read_Int('Num_Protein',200)

    aa=Ini_Read_Double('aa')
    bb=Ini_Read_Double('bb')

    gamma=Ini_Read_Double('gamma',0.005d0)

    kappa_0=Ini_Read_Double('kappa0',0.0005d0)
    kappa_1=Ini_Read_Double('kappa1',0.005d0)


    num_steps=Ini_Read_Int('SnapShots',100)
	t_i=Ini_Read_Double('Initial_Time',1d0)
	t_f=7d0/gamma

	print*,'protein mean:',(aa+bb*kappa_0)/(kappa_0+kappa_1)
	!Ini_Read_Double('Final_Time',20d0)

	allocate(timesteps(num_steps))

	forall(i=1:num_steps) timesteps(i)=t_i+(t_f-t_i)*dble(i-1)/dble(num_steps-1)

    allocate(init_values(num_DNA_states,num_mRNA,num_Protein))
    allocate(final_values(size(timesteps),num_DNA_states,num_mRNA,num_Protein)) 

    ! Initialize all probilities same 1/(num_mRNA*num_Protein)
    !init_values=1d0/dble(num_mRNA*num_Protein)
    init_values=0d0

    
	init_values(1,2,1)=1d0
 

	call solve(init_values,final_values,timesteps)


	do i=1,size(timesteps)	
		write(slice_num,'(i10)') i
		write_filename=trim(adjustl(foldername))//"/"//trim("Genex3_OUT_Protein")//&
		&trim("_timeslice_")//trim(adjustl(slice_num))//".txt"
		OPEN(NEWUNIT=writeUnit,FILE=write_filename,action='write',status='REPLACE')   			
		write(writeUnit,*)'Protein',' ','Probability' 			
        do n=1,num_Protein
			write(writeUnit,*)n-1,sum(final_values(i,:,:,n))
		end do	
		close(writeUnit)
	end do


	do i=1,size(timesteps)	
		write(slice_num,'(i10)') i
		write_filename=trim(adjustl(foldername))//"/"//trim("Genex3_OUT_mRNA")//&
		&trim("_timeslice_")//trim(adjustl(slice_num))//".txt"
		OPEN(NEWUNIT=writeUnit,FILE=write_filename,action='write',status='REPLACE') 
		write(writeUnit,*)'mRNA',' ','Probability'		
        do m=1,num_mRNA
			write(writeUnit,*)m-1,sum(final_values(i,:,m,:))
		end do	
		close(writeUnit)
	end do

	write_filename=trim(adjustl(foldername))//"/"//trim("Genex3_OUT_timeslices.txt")		
	OPEN(NEWUNIT=writeUnit,FILE=write_filename,action='write',status='REPLACE')   			
	do i=1,size(timesteps)	
		write(writeUnit,*)timesteps(i)
	end do
	close(writeUnit)
	print*,'Sum of probabilities:',sum(final_values(size(timesteps),:,:,:))

end program main