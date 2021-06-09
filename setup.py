from setuptools import setup, find_packages
from codecs import open
import os
import glob

__author__ = "The Atomistic.net Team"
__email__ = "nartrith@atomistic.net"
__maintainer__ = "Haoyue Guo, Alexander Urban, Nongnuch Artrith"
__maintainer_email__ = "hg2568@columbia.edu, a.urban@columbia.edu, nartrith@atomistic.net"

here = os.path.abspath(os.path.dirname(__file__))
package_name = 'xas_tools'
package_description = 'Tools for X-ray absorption spectroscopy (XAS)'

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as fp:
    long_description = fp.read()

# Get version number from the VERSION file
with open(os.path.join(here, package_name, 'VERSION')) as fp:
    version = fp.read().strip()

setup(
    name=package_name,
    version=version,
    description=package_description,
    long_description=long_description,
    author=__author__,
    author_email=__email__,
    maintainer=__maintainer__,
    maintainer_email=__maintainer_email__,
    license='MPL',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Programming Language :: Python :: 3'
    ],
    keywords=['materials science', 'spectroscopy'],
    packages=find_packages(exclude=['tests']),
    scripts=glob.glob(os.path.join("scripts", "*.py")),
    install_requires=['numpy', 'pymatgen', 'pyyaml']
)
