module parameters
	USE RKSUITE_90_PREC
	implicit none
	real(wp)::aa,bb,gamma,kappa_0,kappa_1
	integer::num_mRNA=100,num_Protein=200,num_steps=100
	real(wp)::t_i=1d0,t_f=10d0
	integer,parameter::num_DNA_states=2,inactive=1,active=2
end module parameters	