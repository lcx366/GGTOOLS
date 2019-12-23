import setuptools
from setuptools import setup 

setup(
    name='ggtools',
    version='1.1.1',
    long_description_content_type='text/markdown',
    description='A package to handle the GRACE and GRACE-FO GSM data and GLDAS grid data',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    license='MIT',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/GGTOOLS',
    classifiers=[
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'scipy',
        'numpy',
        'matplotlib',
        'cartopy',
        'pyshtools',
        'GPy',
        'xarray',
        'h5py',
        'requests',
        'astropy',
        'pySphericalPolygon',
        'datetime'
        ],
)
