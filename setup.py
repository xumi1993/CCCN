#!/usr/bin/env python
from setuptools import find_packages, setup
packages=find_packages()

VERSION = "1.1"
setup(name='cccn',
      version=VERSION,
      author='Mijian Xu',
      author_email='gomijianxu@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      package_dir={'cccn': 'cccn'},
      install_requires=['obspy',
                        'pyasdf'],
      include_package_data=False,
      zip_safe=False
      )
