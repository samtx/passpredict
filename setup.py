import glob
import pathlib

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
meta = dict()

# Get long description
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()
meta['long_description'] = long_description

# Get version number
with open("passpredict/__init__.py", encoding="utf-8") as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith('__version__'):
            exec(line, {}, meta)
            assert '__version__' in meta
            break

setup(
    name="passpredict",
    version=meta['__version__'],
    packages=[
        'passpredict',
        'passpredict.observers',
        'passpredict.satellites',
    ],
    python_requires=">=3.8",
    install_requires=[
        'sgp4>=2.12',
        'numpy>=1.22',
        'scipy>=1.7.3',
        'orbit-predictor',
        'click',
        'rich',
        'httpx',
        'timezonefinder',
        'backports.zoneinfo;python_version<"3.9"'
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
        '': ['*.txt', '*.rst', '*.dat', '*.csv', '*.md'],
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