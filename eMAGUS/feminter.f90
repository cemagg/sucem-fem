! Last changed: 
! 25 April 2003 DBD: CRS_TO_SKYLINE interface added.
! 27 March 2003 DBD : GAUSS_QUAD_VBF interface extended.  
! 01 May 2002 DBD:  Interface added to subroutine GW_OUTPUT_FIELDS
! and interface to GW_OUTPUT_S_PARAMS corrected.

MODULE feminterface 
  IMPLICIT NONE


INTERFACE
  SUBROUTINE ASSIGN_ELEMENT_ORDERS
  USE adaptive_fem
  USE geometry
  USE nrtype
  END SUBROUTINE ASSIGN_ELEMENT_ORDERS
END INTERFACE


INTERFACE
  SUBROUTINE ASSIGN_ERROR_LEVEL_FLAGS
  USE adaptive_fem
  USE nrtype
  END SUBROUTINE ASSIGN_ERROR_LEVEL_FLAGS
END INTERFACE


INTERFACE
  SUBROUTINE ASSIGN_FACE_AND_EDGE_ORDERS(init_flag)
  USE geometry
  USE problem_info
  USE nrtype
  LOGICAL(LGT), INTENT(IN) :: init_flag
  END SUBROUTINE ASSIGN_FACE_AND_EDGE_ORDERS
END INTERFACE


INTERFACE
  FUNCTION VBF(vbf_type,lambda,grad_lambda,node1,node2,node3)
  USE nrtype
  USE problem_info
  USE unit_numbers
  INTEGER(I4B), INTENT(IN) ::   vbf_type 
  REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
  REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda 
  INTEGER(I4B), INTENT(IN) ::   node1,node2
  INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 
  REAL(SP), DIMENSION (3) :: VBF                   
  END FUNCTION VBF
END INTERFACE


INTERFACE
  FUNCTION VBF_S(vbf_type,lambda,grad_lambda,normal,node1,node2,node3)
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) ::   vbf_type 
  REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
  REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda 
  REAL(SP), DIMENSION(3), INTENT(IN) :: normal
  INTEGER(I4B), INTENT(IN) ::   node1,node2
  INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 
  REAL(SP), DIMENSION (3) :: VBF_S
  END FUNCTION VBF_S
END INTERFACE


INTERFACE
  SUBROUTINE GW_OUTPUT_FIELDS(freqcount)
  USE boundary_conditions
  USE frequency_data
  USE gw_data
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: freqcount
  END SUBROUTINE GW_OUTPUT_FIELDS
END INTERFACE


INTERFACE
  SUBROUTINE GW_OUTPUT_S_PARAMS(freqcount)
  USE boundary_conditions
  USE math_tools, ONLY: PHASE ! Corrected DBD 01 May 02
  USE frequency_data
  USE geometry
  USE gw_data
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: freqcount
  END SUBROUTINE GW_OUTPUT_S_PARAMS
END INTERFACE



INTERFACE
SUBROUTINE BAND_EIGENVALUES
  USE nrtype
  USE geometry
  USE bandwidth
  USE unit_numbers
  USE matrix
  USE problem_info
  END SUBROUTINE BAND_EIGENVALUES
END INTERFACE




INTERFACE 
  SUBROUTINE MAKE_FE_AMATRIX(k0)
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE matrix
  USE nrtype
  USE problem_info
  REAL(SP), INTENT(IN) :: k0
  END SUBROUTINE MAKE_FE_AMATRIX
END INTERFACE




INTERFACE 
  SUBROUTINE CRS_TO_SKYLINE
  USE nrtype
  USE unit_numbers
  USE matrix
  USE problem_info
  END SUBROUTINE CRS_TO_SKYLINE
END INTERFACE


INTERFACE
  SUBROUTINE CONTROL_ADAPTIVE
  USE nrtype
  END SUBROUTINE CONTROL_ADAPTIVE
END INTERFACE


INTERFACE
 SUBROUTINE COUNT_NONZEROS(listsize,list,nonzeros)
 USE nrtype
 USE problem_info
 USE unit_numbers
 INTEGER(I4B),INTENT(IN) :: listsize
 INTEGER(I4B), DIMENSION(listsize),INTENT(IN) :: list
 INTEGER(I4B), INTENT(OUT) :: nonzeros
 END SUBROUTINE COUNT_NONZEROS
END INTERFACE


INTERFACE
  SUBROUTINE DEALLOCATE_ADAPTIVE
  USE adaptive_fem
  USE nrtype
  END SUBROUTINE DEALLOCATE_ADAPTIVE
END INTERFACE


INTERFACE
  SUBROUTINE DIRECT_SOLVE(flag,time_taken)
  USE matrix
  USE unit_numbers
  INTEGER(I4B), INTENT(IN) :: flag
  REAL(SP) :: time_taken
  END SUBROUTINE DIRECT_SOLVE
END INTERFACE


INTERFACE
  SUBROUTINE EIGEN_SYSMAT(timer2,timer3)
  USE nrtype
  USE geometry
  USE unit_numbers   
  USE problem_info
  USE matrix
  REAL(SP), INTENT(OUT) :: timer2, timer3       
  END SUBROUTINE EIGEN_SYSMAT
END INTERFACE


INTERFACE
  SUBROUTINE EIG_MAKE_FE_STMATRICES(bandwidth_determine)
  USE bandwidth
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  LOGICAL(LGT), INTENT(IN) :: bandwidth_determine
  END SUBROUTINE EIG_MAKE_FE_STMATRICES
END INTERFACE


INTERFACE
  SUBROUTINE ENLARGE_SPARSE_MATRICES(estim_nnz,rowstart)
  USE matrix
  INTEGER(I4B), INTENT(INOUT) :: estim_nnz
  INTEGER(I4B), INTENT(IN)  :: rowstart 
  END SUBROUTINE ENLARGE_SPARSE_MATRICES
END INTERFACE


INTERFACE
  SUBROUTINE ERM_INTER_ELEMENT_BC_MAKE(elem1,local_face1,Nvec)
  USE basis_function
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: elem1,local_face1
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE), INTENT(OUT) :: Nvec
  END SUBROUTINE ERM_INTER_ELEMENT_BC_MAKE
END INTERFACE


INTERFACE
  FUNCTION EVALUATE_ELECTRIC_FIELD(elem,s_vec,curl_power,elem_dof_values)
  USE basis_function
  USE geometry !, ONLY: SIMPLEX_COEFFICIENTS
  USE matrix
  USE nrtype
  USE problem_info
  INTEGER(I4B), INTENT(IN) :: elem
  REAL(SP), DIMENSION(4), INTENT(IN) :: s_vec
  INTEGER(I4B), INTENT(IN) :: curl_power
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), OPTIONAL, INTENT(IN) :: elem_dof_values
  COMPLEX(SPC), DIMENSION(3) :: EVALUATE_ELECTRIC_FIELD
  END FUNCTION EVALUATE_ELECTRIC_FIELD
END INTERFACE


INTERFACE
  SUBROUTINE EVALUATE_RESIDUAL_INDICATORS
  USE adaptive_fem
  USE geometry
  USE nrtype
  END SUBROUTINE EVALUATE_RESIDUAL_INDICATORS
END INTERFACE


INTERFACE
  SUBROUTINE EXCHANGE_EDGES(new_edge,old_edge)
  USE problem_info
  USE geometry
  USE unit_numbers
  USE matriX
  INTEGER(I4B), INTENT(IN) :: new_edge, old_edge
  END SUBROUTINE EXCHANGE_EDGES
END INTERFACE 


INTERFACE
  SUBROUTINE exchange_layeredges(current_layer,next_layer,in_this_layer,edge_counter,layer_number)
  USE geometry
  USE unit_numbers
  USE problem_info
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: current_layer
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: next_layer
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: in_this_layer
  INTEGER(I4B),  INTENT(INOUT) :: edge_counter
  INTEGER(I4B),  INTENT(IN) :: layer_number
  END SUBROUTINE exchange_layeredges
END INTERFACE


! Added DBD 07 Aug 04

INTERFACE
  SUBROUTINE FD_SCAT_BVECTOR_SCAT(k0)
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  REAL(SP), INTENT(IN) :: k0    ! wavenumber at which assembly must take place
  END SUBROUTINE FD_SCAT_BVECTOR_SCAT
END INTERFACE

INTERFACE
  SUBROUTINE FD_SCAT_BVECTOR_TOT(k)
  USE boundary_conditions
  USE fd_scat_source
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  REAL(SP), INTENT(IN) :: k  ! wavenumber in external medium 
  END SUBROUTINE FD_SCAT_BVECTOR_TOT
END INTERFACE

INTERFACE
  SUBROUTINE FD_SCAT_INC_SCAT(freq)
  USE basis_function
  USE boundary_conditions
  USE fd_scat_source
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  END SUBROUTINE FD_SCAT_INC_SCAT
END INTERFACE

INTERFACE
  SUBROUTINE FD_SCAT_MATRIX_ALLOCATE
  USE matrix
  USE problem_info
  END SUBROUTINE FD_SCAT_MATRIX_ALLOCATE
END INTERFACE

INTERFACE
  SUBROUTINE FD_SCAT_PRE_X_VEC
  USE basis_function
  USE boundary_conditions
  USE fd_scat_source
  USE matrix
  USE problem_info
  END SUBROUTINE FD_SCAT_PRE_X_VEC
END INTERFACE
! End added DBD 07 Aug 04


INTERFACE
  SUBROUTINE FD_SCAT_SYSMAT
  USE frequency_data
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  END SUBROUTINE FD_SCAT_SYSMAT
END INTERFACE


INTERFACE
  SUBROUTINE FEM_FIELDCALC(xob,yob,zob,E_xyz,H_xyz,eigen_mode)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  REAL(SP), INTENT(IN) :: xob,yob,zob
  COMPLEX(SPC), DIMENSION(3), INTENT(OUT) :: E_xyz,H_xyz
  INTEGER(I4B), INTENT(IN), OPTIONAL :: eigen_mode
  END SUBROUTINE FEM_FIELDCALC
END INTERFACE


INTERFACE
! Changed DBD 29 March 2001
  FUNCTION GAUSS_QUAD_VBF(port_num,element_num,local_face_num,&
         VBF_type,i1,i2,normal,i3,ell)
  USE basis_function, ONLY: VBF
  USE boundary_conditions
  USE geometry
  USE gw_data
  USE math_tools, ONLY: CROSS_PRODUCT
  USE problem_info
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: port_num,element_num,& 
                              local_face_num,VBF_type
!  REAL(SP), INTENT(IN) :: x0
! End DBD changes 29 March 2001
  INTEGER(I4B), INTENT(IN):: i1, i2
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i3
  REAL(SP), INTENT(IN), OPTIONAL :: ell
  REAL(SP) :: GAUSS_QUAD_VBF
  END FUNCTION GAUSS_QUAD_VBF
END INTERFACE

INTERFACE
! Added DBD 16 Dec 2003.
  FUNCTION GAUSS_QUAD_VBF_SPC(ABC_num,element_num,local_face_num,&
                        VBF_type,i1,i2,normal,i3,ell)
  USE basis_function, ONLY: VBF
  USE boundary_conditions
  USE geometry
  USE gw_data
  USE math_tools, ONLY: CROSS_PRODUCT
  USE problem_info
  USE quad_tables
  USE FD_scat_source
  INTEGER(I4B), INTENT(IN) :: ABC_num,element_num,& 
                              local_face_num,VBF_type 
  INTEGER(I4B), INTENT(IN):: i1, i2                                                           
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal 
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i3     
  REAL(SP), INTENT(IN), OPTIONAL :: ell        
  COMPLEX(SPC) :: GAUSS_QUAD_VBF_SPC
  END FUNCTION GAUSS_QUAD_VBF_SPC
END INTERFACE


! Added DBD 04 June 2003 
!INTERFACE
!  FUNCTION GAUSS_QUAD_LINE_VBF(element_num,local_edge_num,&
!                             VBF_type,time_derivative_order,ell)  
!  USE basis_function, ONLY: VBF
!  USE boundary_conditions
!  USE geometry
!  USE problem_info
!  USE quad_tables
!  USE TD_source, ONLY: TD_INC_FIELD,TD_D_BY_DT_INC_FIELD,TD_D2_BY_DT2_INC_FIELD
!  INTEGER(I4B), INTENT(IN) :: element_num,local_edge_num,VBF_type,time_derivative_order
!  REAL(SP), INTENT(IN) :: ell        
!  REAL(SP) :: GAUSS_QUAD_LINE_VBF 
!  END FUNCTION GAUSS_QUAD_LINE_VBF
!END INTERFACE
! End added DBD 04 June 2003 

INTERFACE
! Changed DBD 2 March 2002
FUNCTION GAUSS_QUAD_CFIELD(port_num,element_num,local_face_num,&
                           e1,e2,e3,f1,f2,f3,ell)
  USE boundary_conditions
  USE geometry, ONLY: XYZ_COORDINATES, GRADIENT_LAMBDA, FACE_AREA, &
                      LOCAL_FACENODES
  USE gw_data
  USE problem_info
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: port_num, element_num, & 
                              local_face_num
  COMPLEX(SPC), DIMENSION(3), INTENT(IN) :: e1
  COMPLEX(SPC), DIMENSION(3), INTENT(IN), OPTIONAL :: e2,e3
  COMPLEX(SPC), INTENT(IN), OPTIONAL :: f1,f2,f3
  REAL(SP), DIMENSION(3), INTENT(IN),OPTIONAL  :: ell
  COMPLEX(SPC) :: GAUSS_QUAD_CFIELD
  END FUNCTION GAUSS_QUAD_CFIELD
END INTERFACE


INTERFACE
  SUBROUTINE INPUT_ELEMENT_ORDERS
  USE geometry
  USE nrtype
  END SUBROUTINE INPUT_ELEMENT_ORDERS
END INTERFACE


INTERFACE
  SUBROUTINE ITER_SOLVE(time_taken)
  USE CBAA_data
  USE matrix
  USE problem_info
  USE unit_numbers
  REAL(SP), INTENT(OUT) :: time_taken
  END SUBROUTINE ITER_SOLVE
END INTERFACE


INTERFACE
  SUBROUTINE ITER_SOLVE_DP(time_taken,smallest_its,largest_its,average_its)
  USE CBAA_data
  USE matrix
  USE problem_info
  USE unit_numbers
  REAL(SP), INTENT(OUT) :: time_taken
  INTEGER(I4B), INTENT(OUT), OPTIONAL :: smallest_its,largest_its,average_its
  END SUBROUTINE ITER_SOLVE_DP
END INTERFACE


INTERFACE
  SUBROUTINE LAYER_THE_LAYER(this_layer_size,this_layer,updated_this_layer,layer_number)
  USE geometry
  USE unit_numbers
  USE problem_info
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: this_layer_size
  INTEGER(I4B), DIMENSION(this_layer_size), INTENT(IN) :: this_layer
  INTEGER(I4B) , DIMENSION(this_layer_size), INTENT(OUT) :: updated_this_layer
  INTEGER(I4B), INTENT(IN) :: layer_number
  END SUBROUTINE layer_the_layer
END INTERFACE    


INTERFACE    
  FUNCTION LOCAL_TO_GLOBAL_INDEX_TET(elem,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  INTEGER(I4B), INTENT(IN) :: elem,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_TET
  END FUNCTION LOCAL_TO_GLOBAL_INDEX_TET
END INTERFACE    


INTERFACE    
  FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TET(elem,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  INTEGER(I4B), INTENT(IN) :: elem,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_PRE_TET
  END FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TET
END INTERFACE    


INTERFACE    
  FUNCTION LOCAL_TO_GLOBAL_INDEX_TRI(elem,local_face,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  INTEGER(I4B), INTENT(IN) :: elem,local_face,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_TRI
  END FUNCTION LOCAL_TO_GLOBAL_INDEX_TRI
END INTERFACE    


INTERFACE    
  FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TRI(elem,local_face,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  INTEGER(I4B), INTENT(IN) :: elem,local_face,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_PRE_TRI
  END FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TRI
END INTERFACE    

INTERFACE
  FUNCTION L_SQ_NORM_ELEMENTAL(elem,elem_dof_vec)
  USE basis_function
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: elem
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), OPTIONAL, INTENT(IN) :: elem_dof_vec
  REAL(SP) :: L_SQ_NORM_ELEMENTAL 
  END FUNCTION L_SQ_NORM_ELEMENTAL
END INTERFACE


INTERFACE
  SUBROUTINE MESH_INFO_WRITE 
  USE geometry
  USE matrix
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
  END SUBROUTINE MESH_INFO_WRITE 
END INTERFACE


INTERFACE  
  SUBROUTINE NUMBER_DOF
  USE eigen_analysis_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  END SUBROUTINE NUMBER_DOF
END INTERFACE  


INTERFACE  
  SUBROUTINE NUMBER_PRE_DOF
  USE eigen_analysis_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  END SUBROUTINE NUMBER_PRE_DOF
END INTERFACE  


INTERFACE  
  SUBROUTINE OUTPUT_ELEMENT_ORDERS
  USE adaptive_fem
  USE geometry
  USE nrtype
  END SUBROUTINE OUTPUT_ELEMENT_ORDERS
END INTERFACE  


INTERFACE  
  SUBROUTINE OUTPUT_ERROR_DISTRIBUTION
  USE adaptive_fem
  USE geometry
  USE nrtype
  END SUBROUTINE OUTPUT_ERROR_DISTRIBUTION
END INTERFACE  


INTERFACE  
  SUBROUTINE OUTPUT_FE_RESULTS(feltyp,n_x,n_y,n_z,felkor,x0,y0,z0,delta_x,delta_y,delta_z,i_mode)
  USE problem_info
  USE unit_numbers
  INTEGER(I4B), INTENT(IN) :: feltyp,felkor
  INTEGER(I4B), INTENT(IN) :: n_x,n_y,n_z
  REAL(SP), INTENT(IN) :: x0,y0,z0,delta_x,delta_y,delta_z
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i_mode
  END SUBROUTINE OUTPUT_FE_RESULTS
END INTERFACE  


INTERFACE  
  SUBROUTINE OUTPUT_FF_RESULTS(ntheta,nphi,theta0,phi0,dtheta,dphi,data_type,Einc_abs)
  USE problem_info
  USE unit_numbers
  INTEGER(I4B), INTENT(IN) :: ntheta,nphi
  REAL(SP), INTENT(IN) :: theta0,phi0,dtheta,dphi
  INTEGER(I4B), INTENT(IN) :: data_type
  REAL(SP), INTENT(IN), OPTIONAL :: Einc_abs
  END SUBROUTINE OUTPUT_FF_RESULTS
END INTERFACE  


INTERFACE
  SUBROUTINE REORDER
  USE problem_info
  USE geometry
  END SUBROUTINE REORDER
END INTERFACE   


INTERFACE
 SUBROUTINE remove_multiple(listsize, list, newlist)
 USE unit_numbers
 USE problem_info
 INTEGER(I4B), INTENT(IN) :: listsize
 INTEGER(I4B), DIMENSION(listsize), INTENT(IN)  :: list 
 INTEGER(I4B), DIMENSION(listsize), INTENT(OUT) :: newlist
 END SUBROUTINE remove_multiple 
END INTERFACE


!AddedMMB
INTERFACE
  FUNCTION RESIDUAL_ELEMENT_HOMOGENEOUS(elem)
  USE frequency_data
  USE geometry
  USE nrtype
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: elem
  REAL(SP) :: RESIDUAL_ELEMENT_HOMOGENEOUS
  END FUNCTION RESIDUAL_ELEMENT_HOMOGENEOUS
END INTERFACE


INTERFACE 
  FUNCTION RESIDUAL_FACE(face_num)
  USE frequency_data
  USE geometry
  USE gw_data
  USE math_tools, ONLY: FIND_IN_LIST,CROSS_PRODUCT
  USE nrtype
  USE quad_tables
  INTEGER(I4B), INTENT(IN) :: face_num
  REAL(SP) :: RESIDUAL_FACE
  END FUNCTION RESIDUAL_FACE
END INTERFACE 
!AddedMMB

INTERFACE
   SUBROUTINE MATRIX_SPARSE_ALLOCATE
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   END SUBROUTINE MATRIX_SPARSE_ALLOCATE
END INTERFACE


INTERFACE
   SUBROUTINE find_all_in_list(listsize,list,element,poslist,dof)
   USE unit_numbers
   USE problem_info
   USE nrtype
   INTEGER(I4B), INTENT(IN) :: listsize,dof
   INTEGER(I4B),  INTENT(IN), DIMENSION(listsize):: list
   INTEGER(I4B), INTENT(IN)  :: element
   INTEGER(I4B),  INTENT(OUT), DIMENSION(dof) :: poslist   
   END SUBROUTINE find_all_in_list
END INTERFACE


INTERFACE
   SUBROUTINE fill_ones(listsize,fillsize,fillist,list)
   USE unit_numbers
   USE problem_info
   USE nrtype
   INTEGER(I4B),INTENT(IN) :: listsize, fillsize
   INTEGER(I4B),DIMENSION(fillsize), INTENT(IN)  :: fillist
   INTEGER(I4B), DIMENSION(listsize), INTENT(INOUT)  :: list
   END SUBROUTINE fill_ones
END INTERFACE


INTERFACE
   SUBROUTINE find_greater_than(rowsize,rowind,elemk,nextrow)
   USE unit_numbers
   USE problem_info
   USE nrtype
   INTEGER(I4B),INTENT(IN) :: elemk, rowsize
   INTEGER(I4B), DIMENSION(rowsize), INTENT(IN) :: rowind
   INTEGER(I4B), INTENT(OUT) :: nextrow
   END SUBROUTINE find_greater_than
END INTERFACE


INTERFACE
   FUNCTION convertfillcor(row,col)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: row,col
   INTEGER(I4B) convertfillcor
   END FUNCTION convertfillcor 
END INTERFACE


INTERFACE
   FUNCTION convertcor(row,col,flag)
   USE CBAA_data
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: row,col
   INTEGER(I4B), INTENT(IN), OPTIONAL :: flag
   INTEGER(I4B) convertcor
   END FUNCTION convertcor 
END INTERFACE


INTERFACE MATVECPROD
   SUBROUTINE MATVECPROD_REAL(flag,x,y,n,sig)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: flag,n
! Error corrected 1 May 03
!   REAL(I4B), DIMENSION(n), INTENT(IN) :: x
!   REAL(I4B), DIMENSION(n), INTENT(OUT) :: y
   REAL(SP), DIMENSION(n), INTENT(IN) :: x
   REAL(SP), DIMENSION(n), INTENT(OUT) :: y
   REAL(SP), INTENT(IN),OPTIONAL  :: sig
   END SUBROUTINE MATVECPROD_REAL

   SUBROUTINE MATVECPROD_DP(flag,x_DP,y_DP,n,sig)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: flag,n
   REAL(DP), DIMENSION(n), INTENT(IN) :: x_DP
   REAL(DP), DIMENSION(n), INTENT(OUT) :: y_DP
   REAL(DP), INTENT(IN),OPTIONAL  :: sig
   END SUBROUTINE MATVECPROD_DP

   SUBROUTINE MATVECPROD_COMPLEX(flag,x_c,y_c,n,sig)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: flag,n
   COMPLEX(SPC), DIMENSION(n), INTENT(IN) :: x_c
   COMPLEX(SPC), DIMENSION(n), INTENT(OUT) :: y_c
   REAL(SP), INTENT(IN),OPTIONAL  :: sig
   END SUBROUTINE MATVECPROD_COMPLEX
END INTERFACE 


INTERFACE
   SUBROUTINE uppersolve(n,x,y)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: n
   REAL(SP), DIMENSION(n), INTENT(IN) :: x
   REAL(SP), DIMENSION(n), INTENT(OUT) :: y
   END SUBROUTINE uppersolve
END INTERFACE


INTERFACE
   SUBROUTINE lowersolve(n,m,x)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: n
   REAL(SP), DIMENSION(n), INTENT(IN) :: m
   REAL(SP), DIMENSION(n), INTENT(OUT) :: x
   END SUBROUTINE lowersolve
END INTERFACE


INTERFACE
   SUBROUTINE SYMBOLIC_FILL
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   END SUBROUTINE SYMBOLIC_FILL
END INTERFACE


INTERFACE
   SUBROUTINE uppersolve2(n,x,y)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: n
   REAL(SP), DIMENSION(n), INTENT(IN) :: x
   REAL(SP), DIMENSION(n), INTENT(OUT) :: y
   END SUBROUTINE uppersolve2
END INTERFACE


INTERFACE
   SUBROUTINE lowersolve2(n,m,x)
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   INTEGER(I4B), INTENT(IN) :: n
   REAL(SP), DIMENSION(n), INTENT(IN) :: m
   REAL(SP), DIMENSION(n), INTENT(OUT) :: x
   END SUBROUTINE lowersolve2
END INTERFACE


INTERFACE
   SUBROUTINE DOOLITTLE_LU
   USE matrix
   USE unit_numbers
   USE problem_info
   USE nrtype         
   END SUBROUTINE DOOLITTLE_LU
END INTERFACE


INTERFACE
  SUBROUTINE CCStoCRS
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype        
  END SUBROUTINE CCStoCRS
END INTERFACE


INTERFACE
 SUBROUTINE sparse_ARPACKD(howmny, select, d, z, n, numev, resid, numcv, &
                          v, ldv, iparam, workd, workl, lworkl, iwork, info)
 USE nrtype
 USE matrix
 USE unit_numbers
 USE problem_info
 Character,                 INTENT(IN)      :: howmny
 Logical(LGT),DIMENSION(numcv),INTENT(INOUT):: select
 Integer(I4B),              INTENT(IN)      :: n, numev, &
                         numcv, ldv, lworkl
 INTEGER(I4B) , INTENT(INOUT)               :: info  
 Integer(I4B), DIMENSION(11), INTENT(INOUT) :: iparam
 Integer(I4B), DIMENSION(n), INTENT(INOUT)  :: iwork
 Real(SP), DIMENSION(lworkl),INTENT(INOUT)  :: workl                 
 Real(SP), DIMENSION(numev),  INTENT(INOUT) :: d
 Real(SP), DIMENSION(n) ,    INTENT(INOUT)  :: resid
 Real(SP), DIMENSION(ldv,numcv), INTENT(INOUT) :: v  
 Real(SP), DIMENSION(3*n), INTENT(INOUT)    :: workd   
 Real(SP), DIMENSION(n,numev), INTENT(INOUT) :: z
 END SUBROUTINE SPARSE_ARPACKD
END INTERFACE


INTERFACE
  SUBROUTINE sortcolind(collistsize, collist,  sortedlistsize, sortedlist)
  USE nrtype
  INTEGER(I4B), INTENT(IN) :: collistsize
  INTEGER(I4B), DIMENSION(collistsize), INTENT(IN) :: collist
  INTEGER(I4B),  INTENT(OUT) :: sortedlistsize
  INTEGER(I4B), DIMENSION(collistsize), INTENT(OUT) :: sortedlist
  END SUBROUTINE sortcolind
END INTERFACE


INTERFACE
  SUBROUTINE symtest(passedtest)
  USE nrtype
  USE matrix
  LOGICAL(LGT) , INTENT(OUT) :: passedtest
  END SUBROUTINE symtest
END INTERFACE


INTERFACE
  SUBROUTINE remove_zeros(listsize, list, newlist)
  USE unit_numbers
  USE problem_info
  INTEGER(I4B) :: listsize
  INTEGER(I4B), DIMENSION(listsize), INTENT(IN)  :: list 
  INTEGER(I4B), DIMENSION(listsize), INTENT(OUT) :: newlist
  END SUBROUTINE remove_zeros 
END INTERFACE


! Added DBD 03 April 2003


INTERFACE
  SUBROUTINE TD_MATRIX_ALLOCATE
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
  END SUBROUTINE TD_MATRIX_ALLOCATE
END INTERFACE

INTERFACE
  SUBROUTINE TD_SYSMAT
  USE bandwidth
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE td_source
  USE unit_numbers
  END SUBROUTINE TD_SYSMAT
END INTERFACE

INTERFACE
  SUBROUTINE TD_TIMESTEP
  USE matrix
  USE problem_info
  USE td_source
  USE output_error, ONLY: ERROR_FEMFEKO
  USE unit_numbers
  END SUBROUTINE TD_TIMESTEP
END INTERFACE

! End added DBD 03 April 2003



INTERFACE
! Changed DBD 31 March 2003
  SUBROUTINE U_MAKE_HIERARCHAL(element_num,local_face_num,element_order,port_or_ABC_num,&
                               normal,Us)
! End changed DBD 31 March 2003

  USE geometry
  USE nrtype
  USE math_tools, ONLY: CROSS_PRODUCT
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers

! Changed DBD 31 March 2003
  INTEGER(I4B), INTENT(IN) :: element_num, element_order, & 
                              local_face_num, port_or_ABC_num
!  REAL(SP), INTENT(IN) ::  k0 
! End changed DBD 31 March 2003
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal
  COMPLEX(SPC), INTENT(OUT), DIMENSION(8) :: Us
  END SUBROUTINE U_MAKE_HIERARCHAL
END INTERFACE


END MODULE feminterface



