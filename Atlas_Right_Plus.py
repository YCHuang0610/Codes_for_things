import sys
import numpy as np
import nibabel as nib

def nii_right_plus(nii, plus):
    """
    Add a constant value to the right half of a Nifti image.

    Parameters:
    nii (nibabel.nifti1.Nifti1Image): Input Nifti image.
    plus (float): Value to be added to the right half of the image.

    Returns:
    nibabel.nifti1.Nifti1Image: Nifti image with the right half modified.

    """
    mask = nii.get_fdata() > 0
    nii_data = nii.get_fdata()
    nii_data[int(np.ceil(nii_data.shape[0] / 2)):, :, :] += plus
    nii_data = nii_data * mask
    return nib.Nifti1Image(nii_data, nii.affine, nii.header)

# read Nifti image
img = sys.argv[1]
plus = sys.argv[2]
out = sys.argv[3]

atlas = nib.load(img)
atlas = nii_right_plus(atlas, plus)

# Save Nifti image
nib.save(atlas, out)