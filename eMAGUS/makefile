PROG =	femfeko

SRCS =	bandeig.f90 basisfun.f90 blockinput_tools.f90 cbaa_sys.f90	      \
	 datastrc.f90 debugvar.f90			      \
	eig_sys.f90 eigen_analysis_data.f90 element_matrices.f90	      \
	fd_scat_source.f90 fd_sysmat.f90 femfeko.f90 femfeko_main.f90	      \
	femfuncs.f90 feminter.f90 frequency_data.f90 geometry.f90 gw_sys.f90  \
	input_aux_inputfile.f90 input_aux_materials.f90 input_eigen.f90	      \
	input_frequency.f90 input_mesh.f90 input_problem_info.f90	      \
	input_scattering.f90 integrate.f90 material_properties.f90	      \
	mathtool.f90 meshin.f90 meshin_old.f90 namelist_tools.f90 nrtype.f90  \
	output_data.f90 output_error.f90 output_mesh.f90 parseinput.f90	      \
	post_pro.f90 problem_info.f90 quad_tables.f90 renumber.f90	      \
	reorder.f90 scattering_analysis_data.f90 sparse.f90 sprsrtns.f90      \
	td_source.f90 td_sysmat.f90 twoform_vbf.f90 write_matrix.f90 bessel.f \
	ssband.f umf4_f77wrapper.c geomwrap.f90

OBJS =	bandeig.o basisfun.o blockinput_tools.o cbaa_sys.o 	      \
	 datastrc.o debugvar.o eig_sys.o				      \
	eigen_analysis_data.o element_matrices.o fd_scat_source.o fd_sysmat.o \
	femfeko.o femfeko_main.o femfuncs.o feminter.o frequency_data.o	      \
	geometry.o gw_sys.o input_aux_inputfile.o input_aux_materials.o	      \
	input_eigen.o input_frequency.o input_mesh.o input_problem_info.o     \
	input_scattering.o integrate.o material_properties.o mathtool.o	      \
	meshin.o meshin_old.o namelist_tools.o nrtype.o output_data.o	      \
	output_error.o output_mesh.o parseinput.o post_pro.o problem_info.o   \
	quad_tables.o renumber.o reorder.o scattering_analysis_data.o	      \
	sparse.o sprsrtns.o td_source.o td_sysmat.o twoform_vbf.o	      \
	write_matrix.o bessel.o ssband.o umf4_f77wrapper.o geomwrap.o

INCLUDE = /usr/include/suitesparse
LIBS = -lminpack -larpack -llapack -lf77blas -lumfpack -lamd -lm	

CC = gcc
CFLAGS = -O2 -fPIC
FC = gfortran
FFLAGS = -O2 -fPIC
F90 = gfortran
F90FLAGS = -O2 -fPIC -ffixed-line-length-none 
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $<

%.o : %.mod

bandeig.o: datastrc.o geometry.o nrtype.o output_error.o problem_info.o
basisfun.o: basisfun.o datastrc.o debugvar.o fd_scat_source.o \
	frequency_data.o geometry.o material_properties.o mathtool.o nrtype.o \
	problem_info.o quad_tables.o scattering_analysis_data.o td_source.o
blockinput_tools.o: nrtype.o problem_info.o
cbaa_sys.o: basisfun.o datastrc.o element_matrices.o feminter.o \
	frequency_data.o geometry.o mathtool.o nrtype.o output_error.o \
	problem_info.o quad_tables.o
datastrc.o: nrtype.o
debugvar.o: nrtype.o
eig_sys.o: datastrc.o eigen_analysis_data.o element_matrices.o feminter.o \
	geometry.o nrtype.o output_error.o problem_info.o
eigen_analysis_data.o: nrtype.o problem_info.o
element_matrices.o: basisfun.o datastrc.o element_matrices.o geometry.o \
	integrate.o material_properties.o mathtool.o nrtype.o problem_info.o \
	quad_tables.o
fd_scat_source.o: datastrc.o frequency_data.o geometry.o nrtype.o \
	output_error.o problem_info.o scattering_analysis_data.o
fd_sysmat.o: basisfun.o datastrc.o element_matrices.o fd_scat_source.o \
	feminter.o frequency_data.o geometry.o gw_sys.o mathtool.o nrtype.o \
	output_error.o problem_info.o quad_tables.o
femfeko.o: cbaa_sys.o datastrc.o eigen_analysis_data.o \
	feminter.o frequency_data.o geometry.o gw_sys.o input_mesh.o \
	material_properties.o mathtool.o nrtype.o output_data.o \
	output_error.o output_mesh.o parseinput.o problem_info.o \
	quad_tables.o scattering_analysis_data.o write_matrix.o
femfeko_main.o: femfeko.o geomwrap.o
femfuncs.o: basisfun.o datastrc.o eigen_analysis_data.o element_matrices.o \
	fd_scat_source.o feminter.o geometry.o gw_sys.o material_properties.o \
	mathtool.o nrtype.o output_error.o problem_info.o quad_tables.o \
	scattering_analysis_data.o td_source.o
feminter.o: basisfun.o datastrc.o eigen_analysis_data.o element_matrices.o \
	fd_scat_source.o frequency_data.o geometry.o mathtool.o nrtype.o \
	output_error.o problem_info.o quad_tables.o td_source.o
frequency_data.o: nrtype.o output_error.o
geometry.o: datastrc.o eigen_analysis_data.o material_properties.o mathtool.o \
	nrtype.o output_error.o problem_info.o scattering_analysis_data.o
gw_sys.o: datastrc.o element_matrices.o feminter.o frequency_data.o \
	geometry.o mathtool.o nrtype.o output_error.o problem_info.o
input_aux_inputfile.o: blockinput_tools.o input_aux_materials.o nrtype.o \
	problem_info.o
input_aux_materials.o: material_properties.o nrtype.o problem_info.o
input_eigen.o: eigen_analysis_data.o namelist_tools.o nrtype.o problem_info.o
input_frequency.o: frequency_data.o namelist_tools.o nrtype.o
input_mesh.o: blockinput_tools.o datastrc.o geometry.o nrtype.o \
	problem_info.o
input_problem_info.o: datastrc.o input_eigen.o input_frequency.o \
	input_scattering.o namelist_tools.o nrtype.o problem_info.o \
	quad_tables.o
input_scattering.o: namelist_tools.o nrtype.o problem_info.o \
	scattering_analysis_data.o
integrate.o: nrtype.o quad_tables.o
material_properties.o: nrtype.o output_error.o problem_info.o
mathtool.o: nrtype.o
meshin.o: datastrc.o eigen_analysis_data.o frequency_data.o geometry.o \
	material_properties.o nrtype.o output_error.o parseinput.o \
	problem_info.o quad_tables.o scattering_analysis_data.o td_source.o
meshin_old.o: datastrc.o eigen_analysis_data.o frequency_data.o geometry.o \
	material_properties.o nrtype.o output_error.o parseinput.o \
	problem_info.o quad_tables.o scattering_analysis_data.o td_source.o
namelist_tools.o: nrtype.o output_error.o
output_data.o: datastrc.o geometry.o material_properties.o nrtype.o \
	output_error.o problem_info.o quad_tables.o
output_error.o: datastrc.o nrtype.o
output_mesh.o: datastrc.o geometry.o nrtype.o
parseinput.o: datastrc.o input_aux_inputfile.o input_problem_info.o \
	namelist_tools.o problem_info.o
post_pro.o: basisfun.o cbaa_sys.o datastrc.o element_matrices.o \
	fd_scat_source.o feminter.o frequency_data.o geometry.o gw_sys.o \
	material_properties.o mathtool.o nrtype.o output_error.o \
	problem_info.o quad_tables.o scattering_analysis_data.o td_source.o
problem_info.o: debugvar.o nrtype.o
quad_tables.o: nrtype.o
renumber.o: datastrc.o feminter.o geometry.o mathtool.o nrtype.o \
	problem_info.o
reorder.o: datastrc.o feminter.o geometry.o problem_info.o
scattering_analysis_data.o: nrtype.o problem_info.o
sparse.o: datastrc.o feminter.o nrtype.o output_error.o problem_info.o
sprsrtns.o: cbaa_sys.o datastrc.o feminter.o geometry.o mathtool.o nrtype.o \
	output_error.o problem_info.o
td_source.o: datastrc.o geometry.o nrtype.o output_error.o problem_info.o \
	scattering_analysis_data.o
td_sysmat.o: basisfun.o datastrc.o element_matrices.o feminter.o \
	frequency_data.o geometry.o gw_sys.o material_properties.o nrtype.o \
	output_error.o problem_info.o scattering_analysis_data.o td_source.o
twoform_vbf.o: mathtool.o nrtype.o
write_matrix.o: nrtype.o problem_info.o
