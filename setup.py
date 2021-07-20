from setuptools import Extension, setup
from Cython.Build import cythonize

include_dir = "./include"
library_dir = "./lib"

extensions = [
    Extension("passpredict.timefn", ["passpredict/timefn.pyx"],
        include_dirs=[include_dir],
        libraries=['sgp4','m'],
        library_dirs=[library_dir],
        language="c++"
    ),
]
# from Cython.Build import cythonize
# ext.append(cythonize(['passpredict/timefn_ext.pyx']))

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