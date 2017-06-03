vars = Variables()
vars.Add('INFO', '', 1)
vars.Add('DEBUG', '', 1)
vars.Add('FULL_DEBUG', '', 0)
vars.Add('NAME', '', 0)
vars.Add('ORG', '', 0)
vars.Add('NOISE', '', 0)
vars.Add('R0', '', 0)
vars.Add('BETA', '', 1.0)

env = Environment(variables = vars)
env.Append(CPPDEFINES={
	'INFO' : '${INFO}',
	'DEBUG' : '${DEBUG}', 
	'FULL_DEBUG' : '${FULL_DEBUG}',
	'NAME' : '${NAME}',
	'ORG' : '${ORG}',
	'NOISE' : '${NOISE}',
	'R0' : '${R0}',
	'BETA' : '${BETA}'})


env.Append(CCFLAGS = ['-Wall', '-O3', '-fopenmp'])
env.Append(FORTRANFLAGS = ['-O2', '-march=native'])

env.Append(LIBS = ['libgfortran', 'blas', 'lapack', 'png', 'gomp'])

env.Program(target='unitTests', source=env.Object(source = ["unitTests.cpp", "dc2d.f"]))

env.Program(target='test', source=env.Object(source = ["test.cpp", "io_png.c", "dc2d.f", "augLagL1.cpp"]))