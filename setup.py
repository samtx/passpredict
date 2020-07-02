from setuptools import setup, find_packages

ext = []
# from Cython.Build import cythonize
# ext.append(cythonize(['passpredict/timefn_ext.pyx']))

with open("README.md") as f:
    long_description = f.read()

setup(
    name="passpredict",
    version="0.0.5-alpha.1",
    packages=find_packages(
        exclude=[
            'docs', 'tests'
        ]
    ),
    python_requires=">=3.7", 
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=[
        "numpy",        
        "sgp4 >= 2.8",
        "requests",
        "pydantic",
        "click",
    ],
    extras_require={
        'dev': [
            'Cython',
            'pytest',
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
    author_email="samtx@outlook.com",
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
    extensions = ext,
)