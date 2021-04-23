from distutils.core import setup
from setuptools import find_packages

setup(
	name='gdm',
	version='0.2',
	packages=find_packages(),
	author=open('AUTHORS.md').read(),
	long_description=open('README.md').read(),
	license='GPLv3',
	package_data={'': ['*.csv', '*.pgm']},
	include_package_data=True,
)

