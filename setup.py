from setuptools import setup, find_packages

setup(
    # Module name (lowercase)
    name='MS_QSP',

    description='',

    license='MIT license',

    # author='',

    # author_email='',

    maintainer='Lara Herriott',

    maintainer_email='herriott@maths.ox.ac.uk',

    # Packages to include
    packages=find_packages(include=('QSP_models', 'QSP_models.*')),

    # List of dependencies
    install_requires=[]
)