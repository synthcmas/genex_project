module system_of_equations
implicit none
contains
function system_eq(time,vals)
    USE RKSUITE_90_PREC
    use parameters
    implicit none
    REAL(WP), INTENT(IN) :: time
    REAL(WP), DIMENSION(:), INTENT(IN) :: vals
    REAL(KIND=WP), DIMENSION(SIZE(vals)) :: system_eq
    REAL(WP),DIMENSION(num_DNA_states,num_mRNA,num_Protein)::P_Matrix,dP_Matrix,THIS_DUMMY_VECTOR
    LOGICAL,DIMENSION(num_DNA_states,num_mRNA,num_Protein)::THIS_MASK    
    integer::m,n,mm,nn

    THIS_MASK=.True.
    THIS_DUMMY_VECTOR=0d0

    P_Matrix=UNPACK(vals,THIS_MASK,THIS_DUMMY_VECTOR)

    ! Not considering probabilities of 0 mRNA and 0 Protein

    do mm=1,num_mRNA
        do nn=1,num_Protein
            m=mm-1
            n=nn-1
            dP_Matrix(inactive,mm,nn)=-(kappa_0+gamma*bb*dble(m)+gamma*dble(m)+dble(n))*P_Matrix(inactive,mm,nn)+kappa_1*P_Matrix(active,mm,nn)
            if(nn-1>=1) dP_Matrix(inactive,mm,nn)=dP_Matrix(inactive,mm,nn)+gamma*bb*dble(m)*P_Matrix(inactive,mm,nn-1)
            if(mm+1<=num_mRNA) dP_Matrix(inactive,mm,nn)=dP_Matrix(inactive,mm,nn)+gamma*dble(m+1)*P_Matrix(inactive,mm+1,nn)
            if(nn+1<=num_Protein) dP_Matrix(inactive,mm,nn)=dP_Matrix(inactive,mm,nn)+dble(n+1)*P_Matrix(inactive,mm,nn+1)

            dP_Matrix(active,mm,nn)=-(kappa_1+aa+gamma*dble(m)+gamma*bb*dble(m)+dble(n))*P_Matrix(active,mm,nn)+kappa_0*P_Matrix(inactive,mm,nn)
            if(mm-1>=1) dP_Matrix(active,mm,nn)=dP_Matrix(active,mm,nn)+aa*P_Matrix(active,mm-1,nn)
            if(nn-1>=1) dP_Matrix(active,mm,nn)=dP_Matrix(active,mm,nn)+gamma*bb*dble(m)*P_Matrix(active,mm,nn-1)
            if(mm+1<=num_mRNA) dP_Matrix(active,mm,nn)=dP_Matrix(active,mm,nn)+gamma*dble(m+1)*P_Matrix(active,mm+1,nn)
            if(nn+1<=num_Protein) dP_Matrix(active,mm,nn)=dP_Matrix(active,mm,nn)+dble(n+1)*P_Matrix(active,mm,nn+1)

        end do        
    end do

    system_eq=pack(dP_Matrix,.True.)
end function system_eq  

end module system_of_equations
