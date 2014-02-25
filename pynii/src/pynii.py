'''
This is a Nifti library adapted/extracted from nibabel
It is therefore distributed with a MIT licence:

The MIT License

Copyright (c) 2013 Gianlorenzo Fagiolo 
Copyright (c) 2009-2011 Matthew Brett <matthew.brett@gmail.com>
Copyright (c) 2010-2011 Stephan Gerhard <git@unidesign.ch>
Copyright (c) 2006-2010 Michael Hanke <michael.hanke@gmail.com>
Copyright (c) 2011 Christian Haselgrove <christian.haselgrove@umassmed.edu>
Copyright (c) 2010-2011 Jarrod Millman <jarrod.millman@gmail.com>
Copyright (c) 2011-2012 Yaroslav Halchenko <debian@onerussian.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


TODO:
Currently does not support nifti extensions.
There isn't much sanity checking and Exception handling 

Created on 9 May 2013

@author: gfagiolo
'''

__version__ = "0.9 ($Revision: 7 $)"

import gzip
import logging
import sys

import numpy as np

#===============================================================================
# Nifti class
#===============================================================================

class Nifti1Data(object):
    '''
    classdocs
    '''
    HEADER = None
    DATA = None
    FILENAME = None
    AFFINE = None
    ENDIANNESS = None
    
    def __init__(self, header=None, data=None, filename=None, endianness=None):
        '''
        Constructor
        '''
        if not filename is None:
            self.setFilename(filename)
        
        if not endianness is None:
            self.setEndianness(endianness)
        else:
            self.setEndianness(native_code)
            
        if not header is None:
            self.setHeader(header)
        else:
            #set an empty header
            self.__setEmptyHeader()
            #set identity affine
            self.setAffine(np.eye(4))
            
        if not data is None:
            self.setData(data)
        else:
            #set an empty data array of type 4 (int16)
            self.setData(np.zeros((), search_dtype(4)))
    
    #===========================================================================
    # GETTER/SETTERS
    #===========================================================================
    
    def setHeader(self, hdr):
        self.HEADER = hdr

    def getHeader(self):
        return self.HEADER

    def setData(self, data):
        self.DATA = data
        self.__updateHeaderDataChanged()

    def getData(self):
        return self.DATA

    def setFilename(self, filename):
        self.FILENAME = filename

    def getFilename(self):
        return self.FILENAME
    
    def setAffine(self, mat):
        self.AFFINE = mat
        self.__updateHeaderAffineChanged()
    
    def getAffine(self):
        if self.AFFINE is None:
            self.__computeAffine()
        return self.AFFINE
    
    def setEndianness(self, endianness):
        self.ENDIANNESS = endianness
        
    def getEndianness(self):
        return self.ENDIANNESS
    
    def setScaleSlope(self, scl_slope):
        self.getHeader()['scl_slope'] = scl_slope

    def getScaleSlope(self):
        return self.getHeader()['scl_slope']

    def setScaleIntercept(self, scl_inter):
        self.getHeader()['scl_inter'] = scl_inter

    def getScaleIntercept(self):
        return self.getHeader()['scl_inter']
    
    def setTemporalDimension(self, dt):
        self.getHeader()['pixdim'][4] = dt
    
    def getTemporalDimension(self):
        return self.getHeader()['pixdim'][4]

    def getDim(self):
        return self.DATA.shape
    
    def getPixdim(self):
        ndims = self.getHeader()['dim'][0] 
        return self.getHeader()['pixdim'][1 : ndims + 1]
    
    def setXYZTUnits(self, xyzt):
        self.getHeader()['xyzt_units'] = xyzt 

    def getXYZTUnits(self):
        return self.getHeader()['xyzt_units'] 
    
    def setXYZunit(self, unit):
        timeunit = self.getTunit()
        try:
            if isinstance(unit, str):
                xyz = search_unit(unit)
            elif isinstance(unit, int):
                xyz =unit
            self.setXYZTUnits(SPACE_TIME_TO_XYZT(xyz, timeunit))
        except Exception:
            logging.error("Couldn't change the spatial units, check that the provided codes are correct.")
    
    def getXYZunit(self, as_string=False):
        xyz = XYZT_TO_SPACE(self.getXYZTUnits())
        if as_string:
            return search_unitname_code(xyz)
        else:
            return xyz
        
    def setTunit(self, unit):
        'sets the time dimention unit (could be either string or a nifti code)'
        xyzunit = self.getXYZunit()
        try:
            if isinstance(unit, str):
                t = search_unit(unit)
            elif isinstance(unit, int):
                t =unit
            self.setXYZTUnits(SPACE_TIME_TO_XYZT(xyzunit, t))
        except Exception:
            logging.error("Couldn't change the spatial units, check that the provided codes are correct.")
    
    def getTunit(self, as_string=False):
        'returns the time dimension unit'
        tunit = XYZT_TO_TIME(self.getXYZTUnits())
        if as_string:
            return search_unitname_code(tunit)
        else:
            return tunit
    
    def getDescription(self, as_string=True):
        if as_string:
            return str(self.getHeader()['descrip'])
        else:
            return self.getHeader()['descrip']

    def setDescription(self, val):
        'description field is at most 80 characters'
        if len(val) > 80:
            self.getHeader()['descrip'] = val[:80]
        else:
            #add remaining \000 characters
            self.getHeader()['descrip'] = val + (80-len(val))*chr(0)
            
    #NIBABEL compatibility... i.e. duck typing...
    def get_data(self):
        'returns image data as a numpy array'
        return self.getData()
    
    def get_affine(self):
        'returns the reference frame affine transformation as a numpy array'
        return self.getAffine()
    
    def get_header(self):
        'returns the nifti header'
        return self.getHeader()
    
    #===========================================================================
    # HELPER METHODS
    #===========================================================================
    
    def __setEmptyHeader(self):
        self.setHeader(np.zeros((), header_dtype))
        hdr = self.getHeader()
        hdr.setflags(write=True)
#            ('sizeof_hdr', 'i4'), # 0; must be 348
        hdr['sizeof_hdr'] = 348
#    ('data_type', 'S10'), # 4; unused
#        self.getHeader()['data_type'] = 4
#    ('dim_info', 'u1'),   # 39; MRI slice ordering code
#        self.getHeader()[''] = 
#    ('intent_code', 'i2'),# 68; NIFTI intent code
        hdr['intent_code'] = 0 
#    ('datatype', 'i2'),   # 70; it's the datatype
#       default datatype int16
        hdr['datatype'] = 4
#    ('bitpix', 'i2'),     # 72; number of bits per voxel
        hdr['bitpix'] = 16 
#    ('slice_start', 'i2'),# 74; first slice index
#    ('pixdim', 'f4', (8,)),  # 76; grid spacings (units below)
#    ('vox_offset', 'f4'), # 108; offset to data in image file
        hdr['vox_offset'] = 352.0
#    ('scl_slope', 'f4'),  # 112; data scaling slope
        hdr['scl_slope'] = 1.0
#    ('scl_inter', 'f4'),  # 116; data scaling intercept
        hdr['scl_inter'] = 1.0
#    ('slice_end', 'i2'),  # 120; last slice index
#    ('slice_code', 'u1'), # 122; slice timing order
#    ('xyzt_units', 'u1'), # 123; inits of pixdim[1..4]
#       set spatial unit to mm
        hdr['xyzt_units'] = 2
#    ('descrip', 'S80'),  # 148; any text
        
        hdr['descrip'] = DEFAULT_DESCRIPTION + (80-len(DEFAULT_DESCRIPTION))*chr(0)
#    ('aux_file', 'S24'), # 228; auxiliary filename
#    ('qform_code', 'i2'), # 252; xform code
#    ('sform_code', 'i2'), # 254; xform code #They use code 2
#(1, 'scanner', "NIFTI_XFORM_SCANNER_ANAT"),
        hdr['sform_code'] = 1
#    ('intent_name', 'S16'), # 328; name or meaning of data
#        self.getHeader()[''] = 
#    ('magic', 'S4')      # 344; must be 'ni1\0' or 'n+1\0'
#       set magic to single file nifti
        hdr['magic'] = 'n+1\0'

    def __updateHeaderDataChanged(self):
        #update dim
        hdr = self.getHeader()
        ndim = self.getData().shape#[::-1]
        dim_offset = 1
        dims = np.ones((8,), np.dtype(np.int16))
        dims[0] = len(ndim)
        dims[dim_offset:len(ndim) + dim_offset] = ndim
        
        if np.any(dims != hdr['dim']):
            hdr['dim'] = dims
        datatypecode = search_datatype_code(np.dtype(self.getData().dtype), self.getEndianness())
        if hdr['datatype'] != datatypecode:
            hdr['datatype'] = datatypecode
    
    def __resetOrientation(self):
        hdr = self.getHeader()
        hdr['sform_code'] = 0
        hdr['qform_code'] = 0
        hdr['pixdim'][0] = 1.
        hdr['quatern_b'] = 0.
        hdr['quatern_c'] = 0. 
        hdr['quatern_d'] = 0.
        hdr['qoffset_x'] = 0.
        hdr['qoffset_y'] = 0.
        hdr['qoffset_z'] = 0.        
    
    def __updateHeaderAffineChanged(self):
        #update relevant parts of the header
        self.__resetOrientation()
        hdr = self.getHeader()
        aff = self.getAffine()
        RS = aff[:3, :3]
        pixdim = np.sqrt(np.diag(np.dot(RS.T, RS)))
        hdr['pixdim'][1:4] = pixdim
        hdr['sform_code'] = 1
        hdr['srow_x'] = aff[0, :]
        hdr['srow_y'] = aff[1, :]
        hdr['srow_z'] = aff[2, :]
    
    def __computeAffine(self):
        if self.getHeader()['qform_code'] > 0:
            #use qform to find affine
            self.AFFINE = get_qform(self.getHeader())
        elif self.getHeader()['sform_code'] > 0:
            self.AFFINE = np.eye(4)
            self.AFFINE[0,:] = self.getHeader()['srow_x']
            self.AFFINE[1,:] = self.getHeader()['srow_y']
            self.AFFINE[2,:] = self.getHeader()['srow_z']
        else:
            pixdim = self.getPixdim()
            self.AFFINE = np.eye(4)
            for i in range(3):
                self.AFFINE[i, i] = pixdim[i]
    
    def write(self, filename):
        Nifti1Data.save(self, filename)
    
    def read(self, filename):
        img = Nifti1Data.load(filename)
        self.setHeader(img.getHeader())
        self.setData(img.getData())
        self.setFilename(img.getFilename())
    
    #===========================================================================
    # STATIC METHODS
    #===========================================================================
    
    @staticmethod
    def load(filename):
        'Load a nifti image'
        def read_header(fileobj):
            hdr_bindata = of.read(header_dtype.itemsize)
            hdr = np.ndarray((), dtype=header_dtype, buffer=hdr_bindata)
            endianness = guessed_endian(hdr)
            if endianness == swapped_code:
                #read header again with correct endianness
                dt = header_dtype.newbyteorder(swapped_code)
                hdr = np.ndarray((), dtype=dt, buffer=hdr_bindata)
            hdr.setflags(write=True)            
            return hdr, endianness
        
        of = None            
        try:
            if filename.endswith('.nii.gz'):
                of = gzip.GzipFile(filename, 'rb')
            elif filename.endswith('.nii'):
                of = open(filename, 'rb')
            else:
                raise ValueError('File should be either in .nii or .nii.gz format')        
            #read header
            hdr, endianness = read_header(of)
            #seek to voxel data
            of.seek(int(hdr['vox_offset']))
            img_bindata = of.read()
            dtcode = hdr['datatype']
            np_dtype = np.dtype(search_dtype(dtcode))
            np_dtype = np_dtype.newbyteorder(endianness)
            dim_offset = 1
            ndims = hdr['dim'][0] + dim_offset
            shape = hdr['dim'][dim_offset:ndims]#[::-1]
            imgdata = np.ndarray(shape, np_dtype, buffer=img_bindata, order=LOAD_DATA_ORDER)
#            imgdata = np.frombuffer(img_bindata,
#                                    dtype=np_dtype).reshape(shape).squeeze()
            #.transpose(range(len(imgdata.shape))[::-1])            
        except IndexError:
            raise ValueError('Nifti file is not in valid data type')
        finally:
            if not of is None:
                of.close()
            
        return Nifti1Data(**{'filename':filename,'header':hdr,
                                    'data':imgdata, 'endianness':endianness})

    @staticmethod
    def save(img, filename):
        'Save Nifti1Data img to file with filename'
        
        def save_to_disk(fileobj, filename, imgobj):
            #saving instance
            #update img filename
            imgobj.setFilename(filename)
            #write header
            fileobj.write(imgobj.getHeader().tostring())
            #write 4 bytes
            fileobj.write(4*chr(0))
            #write data                
#            fileobj.write(imgobj.getData().T.tostring(order=LOAD_DATA_ORDER))
            data = imgobj.getData().transpose()
            fileobj.write(data.tostring(order=SAVE_DATA_ORDER))
            
        of = None
        try:            
            if filename.endswith('.nii.gz'):
                of = gzip.GzipFile(filename, 'wb')
            elif filename.endswith('.nii'):
                of = open(filename, 'rb')
            else:
                raise ValueError('ERROR: choose a filename that ends either by .nii.gz or .nii (%s)'%filename)
            save_to_disk(of, filename, img)
        finally:
            if not of is None:
                of.close()


#===============================================================================
# GLOBALS
#===============================================================================

LOAD_DATA_ORDER = 'F'
SAVE_DATA_ORDER = 'C'
DEFAULT_DESCRIPTION = 'pynii ' + __version__

#===============================================================================
# ENDIANNESS
#===============================================================================

sys_is_le = sys.byteorder == 'little'
native_code = sys_is_le and '<' or '>'
swapped_code = sys_is_le and '>' or '<'

def guessed_endian(hdr):
    'guess image endianness by lookin at sizeof_hdr (i4) which should be 348 (or 1543569408 if endianness is swapped)'
    if hdr['sizeof_hdr'] == 1543569408:
        return swapped_code
    else:
        return native_code 
    
#===============================================================================
# QUATERNION FUNCTIONS
#===============================================================================

MAX_FLOAT = np.maximum_sctype(np.float)
FLOAT_EPS = np.finfo(np.float).eps
# Needed for quaternion calculation
FLOAT32_EPS_3 = -np.finfo(np.float32).eps * 3

def fillpositive(xyz, w2_thresh=None):
    ''' Compute unit quaternion from last 3 values

    Parameters
    ----------
    xyz : iterable
       iterable containing 3 values, corresponding to quaternion x, y, z
    w2_thresh : None or float, optional
       threshold to determine if w squared is really negative.
       If None (default) then w2_thresh set equal to
       ``-np.finfo(xyz.dtype).eps``, if possible, otherwise
       ``-np.finfo(np.float).eps``

    Returns
    -------
    wxyz : array shape (4,)
         Full 4 values of quaternion

    Notes
    -----
    If w, x, y, z are the values in the full quaternion, assumes w is
    positive.

    Gives error if w*w is estimated to be negative

    w = 0 corresponds to a 180 degree rotation

    The unit quaternion specifies that np.dot(wxyz, wxyz) == 1.

    If w is positive (assumed here), w is given by:

    w = np.sqrt(1.0-(x*x+y*y+z*z))

    w2 = 1.0-(x*x+y*y+z*z) can be near zero, which will lead to
    numerical instability in sqrt.  Here we use the system maximum
    float type to reduce numerical instability

    Examples
    --------
    >>> import numpy as np
    >>> wxyz = fillpositive([0,0,0])
    >>> np.all(wxyz == [1, 0, 0, 0])
    True
    >>> wxyz = fillpositive([1,0,0]) # Corner case; w is 0
    >>> np.all(wxyz == [0, 1, 0, 0])
    True
    >>> np.dot(wxyz, wxyz)
    1.0
    '''
    # Check inputs (force error if < 3 values)
    if len(xyz) != 3:
        raise ValueError('xyz should have length 3')
    # If necessary, guess precision of input
    if w2_thresh is None:
        try: # trap errors for non-array, integer array
            w2_thresh = -np.finfo(xyz.dtype).eps * 3
        except (AttributeError, ValueError):
            w2_thresh = -FLOAT_EPS * 3
    # Use maximum precision
    xyz = np.asarray(xyz, dtype=MAX_FLOAT)
    # Calculate w
    w2 = 1.0 - np.dot(xyz, xyz)
    if w2 < 0:
        if w2 < w2_thresh:
            raise ValueError('w2 should be positive, but is %e' % w2)
        w = 0
    else:
        w = np.sqrt(w2)
    return np.r_[w, xyz]

def quat2mat(q):
    ''' Calculate rotation matrix corresponding to quaternion

    Parameters
    ----------
    q : 4 element array-like

    Returns
    -------
    M : (3,3) array
      Rotation matrix corresponding to input quaternion *q*

    Notes
    -----
    Rotation matrix applies to column vectors, and is applied to the
    left of coordinate vectors.  The algorithm here allows non-unit
    quaternions.

    References
    ----------
    Algorithm from
    http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion

    Examples
    --------
    >>> import numpy as np
    >>> M = quat2mat([1, 0, 0, 0]) # Identity quaternion
    >>> np.allclose(M, np.eye(3))
    True
    >>> M = quat2mat([0, 1, 0, 0]) # 180 degree rotn around axis 0
    >>> np.allclose(M, np.diag([1, -1, -1]))
    True
    '''
    w, x, y, z = q
    Nq = w*w + x*x + y*y + z*z
    if Nq < FLOAT_EPS:
        return np.eye(3)
    s = 2.0/Nq
    X = x*s
    Y = y*s
    Z = z*s
    wX = w*X; wY = w*Y; wZ = w*Z
    xX = x*X; xY = x*Y; xZ = x*Z
    yY = y*Y; yZ = y*Z; zZ = z*Z
    return np.array(
           [[ 1.0-(yY+zZ), xY-wZ, xZ+wY ],
            [ xY+wZ, 1.0-(xX+zZ), yZ-wX ],
            [ xZ-wY, yZ+wX, 1.0-(xX+yY) ]])

def get_qform_quaternion(hdr):
    ''' Compute quaternion from b, c, d of quaternion

    Fills a value by assuming this is a unit quaternion
    '''
    bcd = [hdr['quatern_b'], hdr['quatern_c'], hdr['quatern_d']]
    # Adjust threshold to fact that source data was float32
    return fillpositive(bcd, FLOAT32_EPS_3)
 
def get_qform(hdr):
    """ Return 4x4 affine matrix from qform parameters in header
    hdr

    Returns
    -------
    affine : None or (4,4) ndarray
        return affine reconstructed from qform
        quaternion.  
    """
    code = hdr['qform_code']
    if code == 0:
        return None
    quat = get_qform_quaternion(hdr)
    R = quat2mat(quat)
    vox = hdr['pixdim'][1:4].copy()
    if np.any(vox) < 0:
        raise ValueError('pixdims[1,2,3] should be positive')
    qfac = hdr['pixdim'][0]
    if qfac not in (-1, 1):
        raise ValueError('qfac (pixdim[0]) should be 1 or -1')
    vox[-1] *= qfac
    S = np.diag(vox)
    M = np.dot(R, S)
    out = np.eye(4)
    out[0:3, 0:3] = M
    out[0:3, 3] = [hdr['qoffset_x'], hdr['qoffset_y'], hdr['qoffset_z']]
    return out

#===============================================================================
# CONSTANTS/DATAYPES ETC
#===============================================================================


header_dtd = [
    ('sizeof_hdr', 'i4'), # 0; must be 348
    ('data_type', 'S10'), # 4; unused
    ('db_name', 'S18'),   # 14; unused
    ('extents', 'i4'),    # 32; unused
    ('session_error', 'i2'), # 36; unused
    ('regular', 'S1'),    # 38; unused
    ('dim_info', 'u1'),   # 39; MRI slice ordering code
    ('dim', 'i2', (8,)),     # 40; data array dimensions
    ('intent_p1', 'f4'),  # 56; first intent parameter
    ('intent_p2', 'f4'),  # 60; second intent parameter
    ('intent_p3', 'f4'),  # 64; third intent parameter
    ('intent_code', 'i2'),# 68; NIFTI intent code
    ('datatype', 'i2'),   # 70; it's the datatype
    ('bitpix', 'i2'),     # 72; number of bits per voxel
    ('slice_start', 'i2'),# 74; first slice index
    ('pixdim', 'f4', (8,)),  # 76; grid spacings (units below)
    ('vox_offset', 'f4'), # 108; offset to data in image file
    ('scl_slope', 'f4'),  # 112; data scaling slope
    ('scl_inter', 'f4'),  # 116; data scaling intercept
    ('slice_end', 'i2'),  # 120; last slice index
    ('slice_code', 'u1'), # 122; slice timing order
    ('xyzt_units', 'u1'), # 123; inits of pixdim[1..4]
    ('cal_max', 'f4'),    # 124; max display intensity
    ('cal_min', 'f4'),    # 128; min display intensity
    ('slice_duration', 'f4'), # 132; time for 1 slice
    ('toffset', 'f4'),   # 136; time axis shift
    ('glmax', 'i4'),     # 140; unused
    ('glmin', 'i4'),     # 144; unused
    ('descrip', 'S80'),  # 148; any text
    ('aux_file', 'S24'), # 228; auxiliary filename
    ('qform_code', 'i2'), # 252; xform code
    ('sform_code', 'i2'), # 254; xform code
    ('quatern_b', 'f4'), # 256; quaternion b param
    ('quatern_c', 'f4'), # 260; quaternion c param
    ('quatern_d', 'f4'), # 264; quaternion d param
    ('qoffset_x', 'f4'), # 268; quaternion x shift
    ('qoffset_y', 'f4'), # 272; quaternion y shift
    ('qoffset_z', 'f4'), # 276; quaternion z shift
    ('srow_x', 'f4', (4,)), # 280; 1st row affine transform
    ('srow_y', 'f4', (4,)), # 296; 2nd row affine transform
    ('srow_z', 'f4', (4,)), # 312; 3rd row affine transform
    ('intent_name', 'S16'), # 328; name or meaning of data
    ('magic', 'S4')      # 344; must be 'ni1\0' or 'n+1\0'
    ]

# Full header numpy dtype
header_dtype = np.dtype(header_dtd)

# datatypes not in analyze format, with codes
_float128t = np.void
_complex256t = np.void

def search_dtype(code):
    return filter(lambda x:x[0]==code, _dtdefs)[0][2]

def search_datatype_code(npdtype, endianness):
    return filter(lambda x:np.dtype(x[2]).newbyteorder(endianness)==npdtype, _dtdefs)[0][0]

_dtdefs = ( # code, label, dtype definition, niistring
    (0, 'none', np.void, ""),
    (1, 'binary', np.void, ""),
    (2, 'uint8', np.uint8, "NIFTI_TYPE_UINT8"),
    (4, 'int16', np.int16, "NIFTI_TYPE_INT16"),
    (8, 'int32', np.int32, "NIFTI_TYPE_INT32"),
    (16, 'float32', np.float32, "NIFTI_TYPE_FLOAT32"),
    (32, 'complex64', np.complex64, "NIFTI_TYPE_COMPLEX64"),
    (64, 'float64', np.float64, "NIFTI_TYPE_FLOAT64"),
    (128, 'RGB', np.dtype([('R','u1'),
                  ('G', 'u1'),
                  ('B', 'u1')]), "NIFTI_TYPE_RGB24"),
    (255, 'all', np.void, ''),
    (256, 'int8', np.int8, "NIFTI_TYPE_INT8"),
    (512, 'uint16', np.uint16, "NIFTI_TYPE_UINT16"),
    (768, 'uint32', np.uint32, "NIFTI_TYPE_UINT32"),
    (1024,'int64', np.int64, "NIFTI_TYPE_INT64"),
    (1280, 'uint64', np.uint64, "NIFTI_TYPE_UINT64"),
    (1536, 'float128', _float128t, "NIFTI_TYPE_FLOAT128"),
    (1792, 'complex128', np.complex128, "NIFTI_TYPE_COMPLEX128"),
    (2048, 'complex256', _complex256t, "NIFTI_TYPE_COMPLEX256"),
    (2304, 'RGBA', np.dtype([('R','u1'),
                    ('G', 'u1'),
                    ('B', 'u1'),
                    ('A', 'u1')]), "NIFTI_TYPE_RGBA32"),
    )

#_dtdefs = ( # code, label, dtype definition, niistring
#    (0, 'none', np.dtype(np.void), ""),
#    (1, 'binary', np.dtype(np.void), ""),
#    (2, 'uint8', np.dtype(np.uint8), "NIFTI_TYPE_UINT8"),
#    (4, 'int16', np.dtype(np.int16), "NIFTI_TYPE_INT16"),
#    (8, 'int32', np.dtype(np.int32), "NIFTI_TYPE_INT32"),
#    (16, 'float32', np.dtype(np.float32), "NIFTI_TYPE_FLOAT32"),
#    (32, 'complex64', np.dtype(np.complex64), "NIFTI_TYPE_COMPLEX64"),
#    (64, 'float64', np.dtype(np.float64), "NIFTI_TYPE_FLOAT64"),
#    (128, 'RGB', np.dtype([('R','u1'),
#                  ('G', 'u1'),
#                  ('B', 'u1')]), "NIFTI_TYPE_RGB24"),
#    (255, 'all', np.dtype(np.void), ''),
#    (256, 'int8', np.dtype(np.int8), "NIFTI_TYPE_INT8"),
#    (512, 'uint16', np.dtype(np.uint16), "NIFTI_TYPE_UINT16"),
#    (768, 'uint32', np.dtype(np.uint32), "NIFTI_TYPE_UINT32"),
#    (1024,'int64', np.dtype(np.int64), "NIFTI_TYPE_INT64"),
#    (1280, 'uint64', np.dtype(np.uint64), "NIFTI_TYPE_UINT64"),
#    (1536, 'float128', np.dtype(_float128t), "NIFTI_TYPE_FLOAT128"),
#    (1792, 'complex128', np.dtype(np.complex128), "NIFTI_TYPE_COMPLEX128"),
#    (2048, 'complex256', np.dtype(_complex256t), "NIFTI_TYPE_COMPLEX256"),
#    (2304, 'RGBA', np.dtype([('R','u1'),
#                    ('G', 'u1'),
#                    ('B', 'u1'),
#                    ('A', 'u1')]), "NIFTI_TYPE_RGBA32"),
#    )

## Make full code alias bank, including dtype column
#data_type_codes = make_dt_codes(_dtdefs)
#
## Transform (qform, sform) codes
xform_codes = ( # code, label, niistring
                       (0, 'unknown', "NIFTI_XFORM_UNKNOWN"),
                       (1, 'scanner', "NIFTI_XFORM_SCANNER_ANAT"),
                       (2, 'aligned', "NIFTI_XFORM_ALIGNED_ANAT"),
                       (3, 'talairach', "NIFTI_XFORM_TALAIRACH"),
                       (4, 'mni', "NIFTI_XFORM_MNI_152"))


def search_unit(name):
    'name is a string, finds code'
    return filter(lambda x:x[1]==name, unit_codes)[0][0]

def search_unitname_code(unitcode):
    'unitcode is an int, finds a string'
    return filter(lambda x:x[0]==unitcode, unit_codes)[0][1]

# unit codes
unit_codes = ( # code, label
    (0, 'unknown'),
    (1, 'meter'),
    (2, 'mm'),
    (3, 'micron'),
    (8, 'sec'),
    (16, 'msec'),
    (24, 'usec'),
    (32, 'hz'),
    (40, 'ppm'),
    (48, 'rads'))

def XYZT_TO_SPACE(xyzt):
    return xyzt & 0x07
 
def XYZT_TO_TIME(xyzt):
    return xyzt & 0x38

def SPACE_TIME_TO_XYZT(ss,tt):
    return XYZT_TO_SPACE(ss) | XYZT_TO_TIME(tt)

#===============================================================================
# #define XYZT_TO_SPACE(xyzt)       ( (xyzt) & 0x07 )
# #define XYZT_TO_TIME(xyzt)        ( (xyzt) & 0x38 )
# #define SPACE_TIME_TO_XYZT(ss,tt) (  (((char)(ss)) & 0x07)   \
#                                   | (((char)(tt)) & 0x38) )
#===============================================================================

#
#slice_order_codes = Recoder(( # code, label
#    (0, 'unknown'),
#    (1, 'sequential increasing', 'seq inc'),
#    (2, 'sequential decreasing', 'seq dec'),
#    (3, 'alternating increasing', 'alt inc'),
#    (4, 'alternating decreasing', 'alt dec'),
#    (5, 'alternating increasing 2', 'alt inc 2'),
#    (6, 'alternating decreasing 2', 'alt dec 2')),
#                            fields=('code', 'label'))
#
#intent_codes = Recoder((
#    # code, label, parameters description tuple
#    (0, 'none', (), "NIFTI_INTENT_NONE"),
#    (2, 'correlation',('p1 = DOF',), "NIFTI_INTENT_CORREL"),
#    (3, 't test', ('p1 = DOF',), "NIFTI_INTENT_TTEST"),
#    (4, 'f test',
#     ('p1 = numerator DOF', 'p2 = denominator DOF'),
#     "NIFTI_INTENT_FTEST"),
#    (5, 'z score', (), "NIFTI_INTENT_ZSCORE"),
#    (6, 'chi2', ('p1 = DOF',), "NIFTI_INTENT_CHISQ"),
#    # two parameter beta distribution
#    (7, 'beta',
#     ('p1=a', 'p2=b'),
#     "NIFTI_INTENT_BETA"),
#    # Prob(x) = (p1 choose x) * p2^x * (1-p2)^(p1-x), for x=0,1,...,p1
#    (8, 'binomial',
#     ('p1 = number of trials', 'p2 = probability per trial'),
#     "NIFTI_INTENT_BINOM"),
#    # 2 parameter gamma
#    # Density(x) proportional to # x^(p1-1) * exp(-p2*x)
#    (9, 'gamma',
#     ('p1 = shape, p2 = scale', 2),
#     "NIFTI_INTENT_GAMMA"),
#    (10, 'poisson',
#     ('p1 = mean',),
#     "NIFTI_INTENT_POISSON"),
#    (11, 'normal',
#     ('p1 = mean', 'p2 = standard deviation',),
#     "NIFTI_INTENT_NORMAL"),
#    (12, 'non central f test',
#     ('p1 = numerator DOF',
#      'p2 = denominator DOF',
#      'p3 = numerator noncentrality parameter',),
#     "NIFTI_INTENT_FTEST_NONC"),
#    (13, 'non central chi2',
#     ('p1 = DOF', 'p2 = noncentrality parameter',), 
#     "NIFTI_INTENT_CHISQ_NONC"),
#    (14, 'logistic',
#     ('p1 = location', 'p2 = scale',),
#     "NIFTI_INTENT_LOGISTIC"),
#    (15, 'laplace',
#     ('p1 = location', 'p2 = scale'),
#     "NIFTI_INTENT_LAPLACE"),
#    (16, 'uniform',
#     ('p1 = lower end', 'p2 = upper end'),
#     "NIFTI_INTENT_UNIFORM"),
#    (17, 'non central t test',
#     ('p1 = DOF', 'p2 = noncentrality parameter'),
#     "NIFTI_INTENT_TTEST_NONC"),
#    (18, 'weibull',
#     ('p1 = location', 'p2 = scale, p3 = power'),
#     "NIFTI_INTENT_WEIBULL"),
#    # p1 = 1 = 'half normal' distribution
#    # p1 = 2 = Rayleigh distribution
#    # p1 = 3 = Maxwell-Boltzmann distribution.
#    (19, 'chi', ('p1 = DOF',), "NIFTI_INTENT_CHI"),
#    (20, 'inverse gaussian',
#     ('pi = mu', 'p2 = lambda'),
#     "NIFTI_INTENT_INVGAUSS"),
#    (21, 'extreme value 1',
#     ('p1 = location', 'p2 = scale'),
#     "NIFTI_INTENT_EXTVAL"),
#    (22, 'p value', (), "NIFTI_INTENT_PVAL"),
#    (23, 'log p value', (), "NIFTI_INTENT_LOGPVAL"),
#    (24, 'log10 p value', (), "NIFTI_INTENT_LOG10PVAL"),
#    (1001, 'estimate', (), "NIFTI_INTENT_ESTIMATE"),
#    (1002, 'label', (), "NIFTI_INTENT_LABEL"),
#    (1003, 'neuroname', (), "NIFTI_INTENT_NEURONAME"),
#    (1004, 'general matrix',
#     ('p1 = M', 'p2 = N'),
#     "NIFTI_INTENT_GENMATRIX"),
#    (1005, 'symmetric matrix', ('p1 = M',), "NIFTI_INTENT_SYMMATRIX"),
#    (1006, 'displacement vector', (), "NIFTI_INTENT_DISPVECT"),
#    (1007, 'vector', (), "NIFTI_INTENT_VECTOR"),
#    (1008, 'pointset', (), "NIFTI_INTENT_POINTSET"),
#    (1009, 'triangle', (), "NIFTI_INTENT_TRIANGLE"),
#    (1010, 'quaternion', (), "NIFTI_INTENT_QUATERNION"),
#    (1011, 'dimensionless', (), "NIFTI_INTENT_DIMLESS"),
#    (2001, 'time series',
#     (),
#     "NIFTI_INTENT_TIME_SERIES",
#     "NIFTI_INTENT_TIMESERIES"), # this mis-spell occurs in the wild
#    (2002, 'node index', (), "NIFTI_INTENT_NODE_INDEX"),
#    (2003, 'rgb vector', (), "NIFTI_INTENT_RGB_VECTOR"),
#    (2004, 'rgba vector', (), "NIFTI_INTENT_RGBA_VECTOR"),
#    (2005, 'shape', (), "NIFTI_INTENT_SHAPE")),
#    fields=('code', 'label', 'parameters', 'niistring'))
