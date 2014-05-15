'''
Runs some tests comparing pynii results to nibabel

Created on 6 Jun 2013

@author: gfagiolo
'''
from pynii import Nifti1Data, logging, np, sys
ROOT_LOGGER = logging.getLogger()
import os
from glob import glob
import nibabel as nb
from datetime import datetime

TEST_DATA_DIR = os.path.join('..', 'testdata')

def run_test(dd=TEST_DATA_DIR):
    '''
    :param dd: str the test data directory 
    '''    
    for fname in glob(os.path.join(dd, '*.nii.gz')):
        if fname.find('.pynii')>-1:
            continue
        print '\n',fname
        nbii = nb.load(fname)
        nbii.get_data()
        nii = Nifti1Data.load(fname)
        shape0 = nbii.get_data().shape
        shape1 = nii.getData().shape
        oname = fname.replace('.nii', '.pynii.nii')
        Nifti1Data.save(nii, oname)
        nii_saved = Nifti1Data.load(oname)
        shape1s = nii_saved.getData().shape
#        os.remove(oname)
        try:
            print '>>> LOAD TEST:'
            print 'nibabel.dim = %s pynii.dim = %s'%(str(shape0), str(shape1))
            assert shape0 == shape1
            ndiffvox = len(np.flatnonzero(nii.getData()-nbii.get_data()))
            print 'How many different voxels?', ndiffvox
            assert ndiffvox == 0
            affine_different = len(np.flatnonzero(nii.getAffine() - nbii.get_affine()))>0
            print 'Are affine different?', affine_different
            try:
                assert affine_different == False
            except AssertionError:
                print nii.getAffine()
                print nbii.get_affine()
                sys.stdout.flush()
                logging.error('nibabel pynii have different affine matrices')
                ROOT_LOGGER.handlers[0].flush()
            print '<<< SAVE TEST:'
            print 'pynii_saved.dim = %s'%(str(shape1s))
            assert shape1s == shape1
            ndiffvox = len(np.flatnonzero(nii.getData()-nii_saved.getData()))
            print 'How many different voxels?', ndiffvox
            assert ndiffvox == 0
            affine_different = len(np.flatnonzero(nii.getAffine() - nii_saved.getAffine()))>0
            print 'Are affine different?', affine_different
            assert affine_different == False
            print '<<< PASSED >>>\n'
        except AssertionError as e:
            print '!!! NOT PASSED !!!'
            sys.stdout.flush()
            logging.error('Test failed %s'%str(e))
            ROOT_LOGGER.handlers[0].flush()
            
def single_test(fname=os.path.join(TEST_DATA_DIR, 'avg152T1_RL_nifti.nii.gz')):   
    t0 = datetime.now()
    nbii = nb.load(fname)
    nbii.get_data()
    t1 = datetime.now()
    print 'nibabel:', t1-t0, nbii.get_affine()
    t0 = datetime.now()
    nii = Nifti1Data.load(fname)
    t1 = datetime.now()
    print 'pynii:', t1-t0, nii.getAffine()
    oname = fname.replace('.nii','.nibabel.nii')
    nb.save(nbii, oname)
    oname = fname.replace('.nii','.pynii.nii')
    nii.write(oname)
    #plot results
    import pylab as pl
    if len(nii.getDim()) == 4:
        pl.subplot(311)    
        pl.imshow(nii.getData()[:,:,15,1].squeeze())
    #    pl.imshow(nii.getData().T[:,:,15,1].squeeze())
        pl.subplot(312)
        pl.imshow(nbii.get_data()[:,:,15,1].squeeze())
        pl.subplot(313)
        pl.imshow(nii.getData()[:,:,15,1].squeeze()-nbii.get_data()[:,:,15,1].squeeze())
    #    pl.imshow(nii.getData().T[:,:,15,1].squeeze()-nbii.get_data()[:,:,15,1].squeeze())
        pl.show()
        pl.close()
    elif len(nii.getDim()) == 3:
        pl.subplot(311)    
        pl.imshow(nii.getData()[:,:,15].squeeze())
    #    pl.imshow(nii.getData().T[:,:,15,1].squeeze())
        pl.subplot(312)
        pl.imshow(nbii.get_data()[:,:,15].squeeze())
        pl.subplot(313)
        pl.imshow(nii.getData()[:,:,15].squeeze()-nbii.get_data()[:,:,15].squeeze())
    #    pl.imshow(nii.getData().T[:,:,15,1].squeeze()-nbii.get_data()[:,:,15,1].squeeze())
        pl.show()
        pl.close()
            
if __name__ == "__main__":    
    run_test()