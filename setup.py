#!/usr/bin/env python

import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='capC-MAP',
                 version='1.1.1',
                 description='Analysis of Capture-C data',
                 long_description=long_description,
                 long_description_content_type="text/x-rst",
                 author='Chris Brackley',
                 author_email='C.Brackley@ed.ac.uk',
                 url='https://git.ecdf.ed.ac.uk/cbrackle/capC-MAP',
                 packages=['capC','capC/resources'],
                 include_package_data=True,
                 install_requires=[
                     'biopython>=1.70',
                     'numpy',
                 ],
                 package_data = {
                     '': ['*.txt'],
                 },
                 classifiers=[
                     "Programming Language :: Python :: 2",
                     "Programming Language :: Python :: 3",
                     "Operating System :: OS Independent",
                 ],
             )
