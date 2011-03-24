#! /bin/sh

rm femfeko_main.so
make

mv femfeko.so femfeko.so.old
# For pre-numpy f2py
#f2py --fcompiler=intel -lminpack -larpack -llapack -lf77blas -lumfpack -lamd \
#-c femfeko.pyf *.o 
# For new f2py
# f2py --fcompiler=intelem -lminpack -larpack -llapack -lf77blas -lumfpack -lamd \
# -c femfeko.pyf *.o 
# Try with gfortran
f2py --fcompiler=gnu95 -lminpack -larpack -llapack -lf77blas -lumfpack -lamd \
-c femfeko.pyf *.o 
