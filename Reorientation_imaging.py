######################################
# by YCH 20240119
# Reorient the position of a monkey from one orientation to another.
# usage: python Reorient_position.py InputFileName.nii.gz OutputFileName.nii.gz init_axcodes(eg: RIA) final_axcodes(eg: RAS)
######################################

import nibabel as nib
import numpy as np
import sys

def compute_orientation(init_axcodes, final_axcodes):
    """Compute the orientation of the image.
    Parameters
    ----------
    init_axcodes : str
        The initial orientation of the image.
    final_axcodes : str
        The final orientation of the image.
    Returns
    -------
    orientation : int
        The orientation of the image.
    """
    init_ornt = nib.orientations.axcodes2ornt(init_axcodes)
    final_ornt = nib.orientations.axcodes2ornt(final_axcodes)
    orientation = nib.orientations.ornt_transform(init_ornt, final_ornt)
    return orientation

def reorient_image(image, orientation):
    """Reorient an image.
    Parameters
    ----------
    image : array, shape (I, J, K)
        The image to reorient.
    orientation : int
        The orientation of the image.
    Returns
    -------
    image : array, shape (I, J, K)
        The reoriented image.
    """
    image = nib.orientations.apply_orientation(image, orientation)
    return image

def reorient_file(fileIn, fileOut, orientation):
    """Reorient a file.
    Parameters
    ----------
    fileIn : str
        The input file.
    fileOut : str
        The output file.
    orientation : int
        The orientation of the image.
    """
    imagineIn = nib.load(fileIn)
    header = imagineIn.header
    affine = imagineIn.affine
    imgData = imagineIn.get_fdata()
    imgReorient = reorient_image(imgData, orientation)
    nib.Nifti1Image(imgReorient, affine=affine, header=header).to_filename(fileOut)

def main():
    """Reorient a file.
    Parameters
    ----------
    fileIn : str
        The input file.
    fileOut : str
        The output file.
    """
    fileIn = sys.argv[1]
    fileOut = sys.argv[2]
    init_axcodes = sys.argv[3]
    final_axcodes = sys.argv[4]
    orientation = compute_orientation(init_axcodes, final_axcodes)
    print(orientation)
    reorient_file(fileIn, fileOut, orientation)

if __name__ == "__main__":
    main()