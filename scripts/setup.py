from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        'MeshExtension',
        sources=['../src/setup.cpp'],
        include_dirs=['/usr/include/petsc', '/usr/include', '/usr/include/python3.10', '/path/to/pybind11/include', '/usr/lib/x86_64-linux-gnu/openmpi/include'],
        library_dirs=['/usr/lib'],
        libraries=['petsc', 'stdc++', 'gmsh', 'python3.10', 'mpi'],
        extra_compile_args=['-fdiagnostics-color=always', '-g']
    ),
]

setup(
    name="MeshExtension",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)