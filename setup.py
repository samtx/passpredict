import glob

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

source_files = []
source_files += [
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
]
source_files += glob.glob("passpredict/ext/sgp4/*.cpp")
source_files += glob.glob("passpredict/ext/*.cpp")
source_files = glob.glob("passpredict/*.pyx")

include_dirs = [
    numpy.get_include(),
    'passpredict/ext',
    'passpredict/ext/sgp4',
    'passpredict/ext/sofa',
]

extensions = [
    Extension("passpredict.timefn", source_files,
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