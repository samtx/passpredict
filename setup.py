import glob

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np


common_kw = {
    'extra_link_args': ['--verbose'],
    'include_dirs': [
        'passpredict',
        'passpredict/sgp4',
        'passpredict/sofa',
        np.get_include()
    ],
    'extra_compile_args': ['-O2'],
    # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
}

source_files = []
source_files += glob.glob('passpredict/sgp4/*.cpp')

extensions = [
    Extension("passpredict._time",
        ['passpredict/_time.pyx'] + sorted(glob.glob('passpredict/sgp4/*.cpp')),
        **common_kw,
    ),
    Extension("passpredict._rotations",
        ['passpredict/_rotations.pyx'] + sorted(glob.glob('passpredict/sofa/*.c')),
        **common_kw,
    ),
    Extension(
        'passpredict._solar',
        ['passpredict/_solar.pyx'],
        **common_kw
    ),
]

with open("README.md") as f:
    long_description = f.read()

setup(
    name="passpredict",
    version="0.1.1",
    packages=['passpredict'],
    python_requires=">=3.9",
    install_requires=[
        'sgp4>=2.12',
        'numpy>=1.22',
        'scipy>=1.7.3',
        'orbit-predictor',
        'pydantic',
        'click',
        'rich',
        'httpx',
        'timezonefinder',
    ],
    extras_require={
        'dev': [
            'Cython',
            'pytest>=6.2.5',
            'pytest-html',
            'pytest-cov',
            'flake8',
        ]
    },
    package_data={
        # If any package contains *.txt, *.rst, *.dat, *.csv files, include them:
        '': ['*.txt', '*.rst', '*.dat', '*.csv'],
    },
    # metadata to display on PyPI
    author="Sam Friedman",
    author_email=None,
    description="Predict upcoming satellite overpasses",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="iss predict orbit sgp4 satellite",
    url="https://www.github.com/samtx/passpredict",   # project home page, if any
    classifiers=[
        'License :: OSI Approved :: Python Software Foundation License'
    ],
    # could also include long_description, download_url, etc.
    zip_safe = False,
    entry_points = {
        'console_scripts': [
            'passpredict = passpredict.cli:main'
        ]
    },
    ext_modules = cythonize(extensions, language_level=3),
)