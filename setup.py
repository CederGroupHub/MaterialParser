from setuptools import setup, find_packages

setup(name='MaterialParser',
      version='0.1.3',
      description='Synthesis Project',
      url='https://github.com/CederGroupHub/MaterialParser',
      author='CederGroup(http://ceder.berkeley.edu)',
      packages=find_packages(),
      install_requires=[
          'sympy',
          'pubchempy',
          'chemdataextractor==1.3.0',
          'regex',
      ],
      zip_safe=False)
