"""Setup file and install scripts.
Version 1.0.0 
"""
import os

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


DESCRIPTION = "ASEkit: A kit apply for ASE analysis pipeline"
DISTNAME = 'ASEkit'
MAINTAINER = 'huangke'
MAINTAINER_EMAIL = 'huangke@shanghaitech.edu.cn'
LICENSE = 'BSD (3-clause)'
#DOWNLOAD_URL = ''
VERSION = "1.2.1"

ROOT_DIR = os.path.split(os.path.realpath(__file__))[0]

if __name__ == "__main__":
    setup(name=DISTNAME,
          author=MAINTAINER,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=(open(ROOT_DIR + "/README.rst").read()),
          license=LICENSE,
          #url=URL,
          #download_url=DOWNLOAD_URL,
          packages=find_packages(),
          install_requires=[
              'pandas',
              'statsmodels',
              'scipy',
              ],
          version=VERSION,
          include_package_data=True,
          # scripts = [],
          entry_points={
              "console_scripts": [
                  'ASEkit = ASEkit:main'
              ]
          },
          classifiers=[
             'Intended Audience :: Science/Research',
             'Programming Language :: Python :: 3.5',
             'License :: OSI Approved :: BSD License',
             'Topic :: Scientific/Engineering :: Bio-Informatics',
             'Operating System :: POSIX',
             'Operating System :: POSIX :: Linux',
             'Operating System :: MacOS'],
          )
