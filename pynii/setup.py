'''
Created on 15 May 2014

@author: glf12
'''
from distutils.core import setup
import os

setup(name='pynii',
      version='0.9',
      description='A simple python nifti interface derive from nibabel',
      author='Gianlorenzo',
      url='https://github.com/gianlo',
      py_modules=[os.path.join('src','pynii')])
