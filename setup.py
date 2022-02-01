import glob

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np


common_kw = {
    'extra_link_args': ['--verbose'],
    'include_dirs': [
        'src/passpredict',
        'src/passpredict/sgp4',
        'src/passpredict/sofa',
        np.get_include()
    ],
    'extra_compile_args': ['-O2'],
    # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
}


extensions = [
    Extension("passpredict._time",
        ['src/passpredict/_time.pyx'] + sorted(glob.glob('src/passpredict/sgp4/*.cpp')),
        **common_kw,
    ),
    Extension("passpredict._rotations",
        ['src/passpredict/_rotations.pyx'] + sorted(glob.glob('src/passpredict/sofa/*.c')),
        **common_kw,
    ),
    Extension(
        'passpredict._solar',
        ['src/passpredict/_solar.pyx'],
        **common_kw
    ),
]


setup(
    ext_modules = cythonize(extensions, language_level=3),
)