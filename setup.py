from setuptools import setup, find_packages

setup(
    name='EvolveX',
    version='1',
    description='EvolveX antibody design pipeline',
    author='Gabriel Cia and Rob Van Der Kant - SwitchLab',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.9',
    install_requires=['pyyaml', 'biopython', 'pandas', 'dask', 'distributed', 'dask-jobqueue'],
    entry_points={
        'console_scripts': [
            'evolvex = evolvex.main:main',
        ]
    },
)