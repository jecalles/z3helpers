import setuptools
import os

with open('README.md', 'r') as handle:
    long_description = handle.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name='z3helpers-jecalles',
    version='0.0.1',
    author='Jonathan Calles',
    author_email='callesjonathan@gmail.com',
    description='ergonomics around z3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # url='https://github.com/jecalles/z3helpers',
    packages=setuptools.find_packages(),
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>3.8',
    install_requires=requirements,
)
