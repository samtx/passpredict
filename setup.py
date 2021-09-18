import glob

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

sofa_files = [
    'passpredict/ext/sofa/gmst82.c',
    'passpredict/ext/sofa/utctai.c',
    'passpredict/ext/sofa/taitt.c',
    'passpredict/ext/sofa/obl80.c',
    'passpredict/ext/sofa/nut80.c',
    'passpredict/ext/sofa/eqeq94.c',
    'passpredict/ext/sofa/anp.c',
    'passpredict/ext/sofa/jd2cal.c',
    'passpredict/ext/sofa/dat.c',
    'passpredict/ext/sofa/cal2jd.c',
    'passpredict/ext/sofa/anpm.c',
    'passpredict/ext/sofa/dtf2d.c',
    'passpredict/ext/sofa/d2dtf.c',
    'passpredict/ext/sofa/d2tf.c',
]
sgp4_files = glob.glob('passpredict/ext/sgp4/*.cpp')
ast2body_files = glob.glob('passpredict/ext/ast2body/*.cpp')
passpredict_cpp_files = glob.glob("passpredict/ext/*.cpp")

include_dirs = [
    numpy.get_include(),
    'passpredict',
    'passpredict/ext',
    'passpredict/ext/sgp4',
    'passpredict/ext/sofa',
    'passpredict/ext/ast2body',
]

extensions = [
    Extension("passpredict.timefn",
        sgp4_files + ['passpredict/timefn.pyx'],
        include_dirs=include_dirs,
        language="c++",
    ),
    Extension("passpredict.predict",
        sofa_files + sgp4_files + passpredict_cpp_files + ['passpredict/predict.pyx'],
        include_dirs=include_dirs,
        language="c++",
    ),
    Extension("passpredict.propagate",
        ast2body_files + ['passpredict/propagate.pyx'],
        include_dirs=include_dirs,
        language="c++",
    ),
]

with open("README.md") as f:
    long_description = f.read()

setup(
    name="passpredict",
    version="0.0.12",
    packages=['passpredict'],
    python_requires=">=3.8",
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=[
        "numpy",
        "requests",
        "pydantic",
        "click",
    ],
    extras_require={
        'dev': [
            'Cython',
            'pytest == 5.4.3',
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
    ext_modules = cythonize(extensions),
)