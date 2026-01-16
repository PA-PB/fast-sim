from setuptools import setup, Extension, find_packages
import pybind11
import os

HERE = os.path.abspath(os.path.dirname(__file__))
eigen_path = os.path.join(HERE, "extern", "eigen")

if os.name == 'nt':
    compile_args = ['/O2', '/arch:AVX2', '/std:c++17', '/openmp']
else:
    compile_args = ['-O3', '-march=native', '-std=c++17', '-fopenmp']

ext_modules = [
    Extension(
      
        'fast_spice._core', 
        [
            'src/bindings.cpp', 
            'src/Netlist.cpp', 
            'src/Solver.cpp'
        ],
        include_dirs=[pybind11.get_include(), eigen_path, 'src', '.'],
        language='c++',
        extra_compile_args=compile_args,
        define_macros=[('NDEBUG', None), ('EIGEN_NO_DEBUG', None)],
        py_limited_api=True,
    ),

]

setup(
    name='fast_spice',
    version='1.0',
    author='Pedro Augusto Pappis Bandeira',
    author_email='pedro.pappisbandeira@gmail.com',
    description='Simulador SPICE C++ com wrapper Python.',
    
  
    packages=find_packages(),
    
    # Compila a extensão
    ext_modules=ext_modules,


    package_data={
        'fast_spice': ['*.pyi', 'py.typed']
    },
    include_package_data=True,
    zip_safe=False,
)