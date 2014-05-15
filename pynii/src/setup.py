'''
Created on 15 May 2014

@author: glf12
'''
from distutils.core import setup
import os

setup(name='pynii',
      version='0.9',
      description='A simple python nifti interface derived from nibabel',
      author='Gianlorenzo Fagiolo',
      url='https://github.com/gianlo',
      license='License :: OSI Approved :: MIT License',
      requires= ['numpy'],
      py_modules=['pynii'])
