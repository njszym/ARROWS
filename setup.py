from setuptools import setup, find_packages
from shutil import copytree
import os


# Install dependencies
setup(
    name='ARROWS',
    version='0.0.1',
    description='A package designed to guide solid-state synthesis experiments toward maximal yield of a target phase.'
    author='Nathan J. Szymanski',
    author_email='nathan_szymanski@berkeley.edu',
    python_requires='>=3.9.0',
    url='https://github.com/njszym/ARROWS',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'pymatgen', 'scipy', 'mp_api']
)
