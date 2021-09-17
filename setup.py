from setuptools import setup, find_packages

setup(
    name='Materials Parser',
    packages=find_packages(),
    version='7.0.0',
    author='Synthesis Project Team',
    author_email='0lgaGkononova@yandex.ru',
    description='Material Synthesis Project',
    zip_safe=False,
    install_requires=[
        "regex",
        "pubchempy",
        "sympy"
    ],
    include_package_data=True
)
