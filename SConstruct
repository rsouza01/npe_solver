
# Comment lines start with the # symbol
# The following sets up an Compile Environment Object with gfortran as the linker.
env = Environment(tools=['default','gfortran'],F90='gfortran',LINK='gfortran',LINKFLAGS='-g',F90FLAGS='-g')

env.VariantDir('build', 'src', duplicate=0)

# The next line of code is an array of the source files names used in the program.
sources = ['build/cgs_constants.f90','build/math_constants.f90','build/npe_functions.f90','build/npe_solver.f90']

# The next line is the actual code that links the executable. env.Program is generates an executable.
objs = env.Program('build/npe_solver', sources)
