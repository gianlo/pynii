'''
Created on 15 May 2014

@author: glf12
'''
from distutils.core import setup

setup(name='pynii',
      version='0.9',
      description='A simple python nifti interface derived from nibabel',
      author='Gianlorenzo Fagiolo',
      url='https://github.com/gianlo/pynii',
      license='License :: OSI Approved :: MIT License',
      requires= ['numpy', 'dicom'],
      py_modules=['pynii'])
