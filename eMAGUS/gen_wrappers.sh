
	#! /bin/sh

rm tmp_femfeko.pyf 

f2py -m femfeko -h tmp_femfeko.pyf \
    nrtype.f90 \
    problem_info.f90 \
    femfeko.f90 \
    geometry.f90 \
    post_pro.f90\
    datastrc.f90\
    geomwrap.f90\
    twoform_vbf.f90\
    basisfun.f90\
    parseinput.f90\
    only: \
    run_femfeko \
    wrap_dofmap \
    get_n xyz_to_elnum \
    xyz_coordinates \
    gradient_lambda \
    simplex_coordinates \
    fem_fieldcalc\
    init_geom \
    twoform_vbf_face \
    vbf : 

grep -v type\( tmp_femfeko.pyf | sed s/private,public/public/ > femfeko.pyf