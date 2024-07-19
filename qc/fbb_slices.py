import os
import sys
import shutil
import warnings
import subprocess
import numpy as np
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from nibabel.orientations import io_orientation, axcodes2ornt
from matplotlib.colors import LinearSegmentedColormap

# Importing the necessary classes from the rablabqc package
rablab_pkg_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(rablab_pkg_path)

from plotter import QCImageGenerator
from processing import ImageProcessor

#  ________________________________________________________________ FILE PATHS _______________________________________________________________ #
tpm_file = os.path.join(rablab_pkg_path,'TPM.nii')

rtpm_path = os.path.join(rablab_pkg_path, 'reslice', 'rT1.nii') # Resliced T1 file to 1mm isotropic resolution
reslice_matlab_script = os.path.join(rablab_pkg_path,'reslice', 'reslice.m')

mask_reslice_matlab_script = os.path.join(rablab_pkg_path, 'reslice', 'mask_reslice.m')

#tmp_folder = os.path.join('/shared/petcore/Projects/LEADS/data_f7p1/summary/piyush_qc/tmp/')
tmp_folder = os.path.join('/tmp/')
#  _____________________________________________________________ CUSTOM COLORMAPS _____________________________________________________________ #
# Notes: 
# 1. The custom colormaps are created using the LinearSegmentedColormap class from matplotlib.colors.
# 2. The alpha = 0 is to set the transparency for the zeros in the image.

def create_colormap(color):
    """
    Parameters
    ----------
    color : tuple
        The RGB values for the color.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        The custom colormap.
    """
    cmap = LinearSegmentedColormap.from_list('custom', [color, color], N=256)
    cmap.set_under(alpha=0)
    return cmap

cmap_red = create_colormap((1, 0, 0))
cmap_yellow = create_colormap((1, 1, 0))
cmap_blues = create_colormap((0, 0, 1))
cmap_pink = create_colormap((1, 0, 1))
cmap_cyan = create_colormap((0, 1, 1))
cmap_orange = create_colormap((1, 0.5, 0))
cmap_green = create_colormap((0, 1, 0))

# create a darker green colormap
cmap_green_dark = create_colormap((0, 0.5, 0))

# create a darker green then cmap_green_dark
cmap_green_darker = create_colormap((0, 0.25, 0))

cmap_turbo = plt.cm.turbo
cmap_turbo.set_under(alpha=0.1)  # Set transparency for zeros


cmap_yellow2 = LinearSegmentedColormap.from_list('custom_color', [(249/255, 195/255, 31/255), (249/255, 195/255, 31/255)], N=256)
cmap_yellow2.set_under(alpha=0)  # Set transparency for zeros

cmap_blue2 = LinearSegmentedColormap.from_list('custom_color', [(46/255, 69/255, 184/255), (46/255, 69/255, 184/255)], N=256)
cmap_blue2.set_under(alpha=0)  # Set transparency for zeros

# _______________________________________________________________ CUSTOM COLORBAR _____________________________________________________________ #
## Important Note: The old versipn of add_colorbar function has been commented out as it is not being used in the current version of the code.

def add_colorbar(fig, ax, cmap='turbo', vmin=0, vmax=2, ticks=[0, 2], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76):
    """
    Add a colorbar to a given axis.

    Parameters:
        fig (matplotlib.figure.Figure): The figure.
        ax (matplotlib.axes.Axes): The axis to add the colorbar to.
        cmap (str): The colormap to use.
        vmin (float): Minimum value of the colorbar.
        vmax (float): Maximum value of the colorbar.
        ticks (list): List of tick values for the colorbar.
        orientation (str): Orientation of the colorbar ('horizontal' or 'vertical').
    """
    sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])  # Dummy array

    # Get position of the axis
    pos = ax.get_position()

    if orientation == 'horizontal':
        cbar_x = cbar_x  # Slightly to the right of the axis
        cbar_y = pos.y0 - cbar_height - 0.00  # Adjust this value to move colorbar closer or farther from the axis
    else:
        cbar_x = cbar_x
        cbar_y = pos.y0
        cbar_width, cbar_height = cbar_height, cbar_width  # Swap width and height for vertical colorbar

    cbar_ax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation=orientation)

    # Set tick values
    cbar.set_ticks(ticks)

    # Set tick labels color to white
    cbar.ax.tick_params(colors='white', labelsize=10)

    return cbar
#  ______________________________________________________________ Setting Scales  _____________________________________________________________ #
mri_vmin = 0
mri_vmax = 160
pet_vmax = 2.5

# ______________________________________________________________ TEMPLATE SLICES  ____________________________________________________________ #

def load_tpm(path, orientation="LAS"):
        """
        Load nifti image with specified orientation
        """
        from nibabel.orientations import io_orientation, axcodes2ornt

        img = nib.load(path)
        img_ornt = io_orientation(img.affine)
        new_ornt = axcodes2ornt(orientation)
        img = img.as_reoriented(img_ornt)
        return img.get_fdata()

tpm_image = load_tpm(tpm_file)[:,:,:,0]

template_axial_slices = [14, 25, 36, 47, 58, 69, 80, 91]
template_sagittal_slices = [44, 55, 77]
template_coronal_slices = [48, 70, 83]

#  _____________________________________________________________ DEFINE FUNCTIONS _____________________________________________________________ #

class FBBQCplots:
    def __init__(self, suvr_img, axial_slices, sagittal_slices, coronal_slices,
                 nu_img = None, aparc_img = None, c1_img = None, 
                 reference_region_1 = None, reference_region_2 = None, reference_region_3 = None,
                 affine_nu_img = None, affine_suvr_img = None,
                 warped_nu_img=None, warped_suvr_img = None,
                 crop_neck = True):
        
        self.suvr_img = suvr_img
        
        self.nu_img = nu_img
        self.aparc_img = aparc_img
        self.c1_img = c1_img
        self.affine_nu_img = affine_nu_img
        self.warped_nu_img = warped_nu_img
        
        self.warped_suvr_img = warped_suvr_img
        self.affine_suvr_img = affine_suvr_img

        self.reference_region_1 = reference_region_1
        self.reference_region_2 = reference_region_2
        self.reference_region_3 = reference_region_3

        self.axial_slices = axial_slices
        self.sagittal_slices = sagittal_slices
        self.coronal_slices = coronal_slices

        self.crop_neck = crop_neck

        if self.crop_neck and self.aparc_img is None:
            self.crop_neck = False
            warnings.warn("The aparc_img is not provided. The neck will not be cropped.", UserWarning)

    # ___________________________________________________ LOAD/READ THE IMAGES ___________________________________________________ #
    def load_nii(self,path, orientation="LAS"):
        """
        Load nifti image with specified orientation
        """

        img = nib.load(path)
        img_ornt = io_orientation(img.affine)
        new_ornt = axcodes2ornt(orientation)
        img = img.as_reoriented(img_ornt)
        return img.get_fdata()

    def generate_matlab_script(self, path, output_script_path):
        """
        This function generates a MATLAB script to reslice the image.
        """
        #with open('/home/mac/pmaiti/Desktop/leads_qc/reslice_test/reslice.m', 'r') as template_file:
        with open(reslice_matlab_script, 'r') as template_file:
            script_content = template_file.read()

        # Replace placeholders with actual paths
        script_content = script_content.replace('<RTPM_PATH>', rtpm_path)
        script_content = script_content.replace('<DATA_PATH>', path)
        
        # Write the modified script to the output path
        with open(output_script_path, 'w') as script_file:
            script_file.write(script_content)


    def generate_mask_reslice_mtlb(self, path, output_script_path):
        """
        This function generates a MATLAB script to reslice the mask image.
        """
        #with open('/home/mac/pmaiti/Desktop/leads_qc/reslice_test/mask_reslice.m', 'r') as template_file:
        with open(mask_reslice_matlab_script, 'r') as template_file:
            script_content = template_file.read()

        script_content = script_content.replace('<RTPM_PATH>', rtpm_path)
        script_content = script_content.replace('<DATA_PATH>', path)
        
        # Write the modified script to the output path
        with open(output_script_path, 'w') as script_file:
            script_file.write(script_content)


    def load_nii_resliced(self, path, orientation="LAS", mask=False):
        """
        Load nifti image with specified orientation
        """
        
        id = path.split('/')[-1].split('.')[0]

        resliced_image_path = os.path.join(tmp_folder, id, 'qc' + id + '.nii')
        
        if not os.path.exists(resliced_image_path):
            
            # Check if the temporary folder with the id exists
            tmp_id_folder = os.path.join(tmp_folder, id)
            if not os.path.exists(tmp_id_folder):
                os.makedirs(tmp_id_folder)
            
            # Copy the image to the temporary id folder
            tmp_file = os.path.join(tmp_id_folder, id + '.nii')
            shutil.copy2(path, tmp_file)
            print(tmp_file)
            
            if mask:
                output_script_path = os.path.join(tmp_id_folder, 'mask_reslice.m')
                self.generate_mask_reslice_mtlb(tmp_file, output_script_path)
            else:
                output_script_path = os.path.join(tmp_id_folder, 'reslice.m')
                self.generate_matlab_script(tmp_file, output_script_path)
                
            # Command to run the MATLAB script
            command = f"matlab -nodisplay -nosplash -r \"run('{output_script_path}');exit;\""
            
            # Run the command
            matprocess = subprocess.run(command, shell=True, capture_output=True, text=True)
            print("Output:\n", matprocess.stdout)

        else:
            print("Resliced image exists for ", id, ". Loading the resliced image...")
            
        img = nib.load(resliced_image_path)
        img_ornt = io_orientation(img.affine)
        new_ornt = axcodes2ornt(orientation)
        img = img.as_reoriented(img_ornt)
        return img.get_fdata()
    
    def load_c1_nii_resliced(self,path, orientation="LAS"):
        """
        Load nifti image with specified orientation
        """
        id = path.split('/')[-1].split('.')[0]
    
        reslice_mask_path = os.path.join(tmp_folder, id, 'qc'+"mask_"+id+".nii")
        
        if not os.path.exists(reslice_mask_path):
            
            # Check if the temporary folder with the id exists
            tmp_id_folder = os.path.join(tmp_folder, id)
            if not os.path.exists(tmp_id_folder):
                os.makedirs(tmp_id_folder)
            
            # Copy the image to the temporary id folder
            tmp_file = os.path.join(tmp_id_folder, id + '.nii')
            shutil.copy2(path, tmp_file)
            print(tmp_file)

            # Loading the nifti image to create a mask of the image
            img = nib.load(tmp_file)
            img_data = img.get_fdata()
            img_affine = img.affine
            img_header = img.header

            # Creating a mask of the image by thresholding
            # Set values to 1 between the threshold of 0.3 and 1
            mask = np.zeros_like(img_data)
            mask[(img_data >= 0.3) & (img_data <= img_data.max())] = 1
            mask = mask.astype(np.uint8)

            # Save the mask as a nifti image
            mask_img = nib.Nifti1Image(mask, img_affine, img_header)
            mask_path = os.path.join(tmp_id_folder,"mask_"+id+".nii")
            nib.save(mask_img, mask_path)
            
            output_script_path = os.path.join(tmp_id_folder, 'mask_reslice.m')
            
            self.generate_mask_reslice_mtlb(mask_path, output_script_path)
            
            # Command to run the MATLAB script
            command = f"matlab -nodisplay -nosplash -r \"run('{output_script_path}');exit;\""
            
            # Run the command
            matprocess = subprocess.run(command, shell=True, capture_output=True, text=True)
            print("Output:\n", matprocess.stdout)

        else:
            print("Resliced image already exists for ", id, ". Loading the resliced image...")
            
        # Resliced mask path
        img = nib.load(reslice_mask_path)
        img_ornt = io_orientation(img.affine)
        new_ornt = axcodes2ornt(orientation)
        img = img.as_reoriented(img_ornt)
        return img.get_fdata()
    

    def load_images(self):
        """
        Reads the image file and returns the image data as a numpy array.

        Parameters
        ----------
        img_path : str
            The path to the image file.

        Returns
        -------
        numpy.ndarray
            The image data as a numpy array.
        """
        self.basename = self.suvr_img.split('/')[-1].split('_')[0]+'_'+self.suvr_img.split('/')[-1].split('_')[1]+'_'+self.suvr_img.split('/')[-1].split('_')[2]
        self.suvr_img_filename = os.path.basename(self.suvr_img)

        self.suvr_img = self.load_nii_resliced(self.suvr_img)

        if self.nu_img is not None:
            self.nu_img_filename = os.path.basename(self.nu_img)
            self.nu_img = self.load_nii_resliced(self.nu_img)

        if self.aparc_img is not None:
            self.aparc_img_filename = os.path.basename(self.aparc_img)
            self.aparc_img = self.load_nii_resliced(self.aparc_img, mask=True)

        if self.c1_img is not None:
            self.c1_img_filename = os.path.basename(self.c1_img)
            self.c1_img = self.load_c1_nii_resliced(self.c1_img)

        if self.reference_region_1 is not None:
            self.reference_region_1_filename = os.path.basename(self.reference_region_1)
            self.reference_region_1 = self.load_nii_resliced(self.reference_region_1, mask=True)

        if self.reference_region_2 is not None:
            self.reference_region_2_filename = os.path.basename(self.reference_region_2)
            self.reference_region_2 = self.load_nii_resliced(self.reference_region_2, mask=True)

        if self.reference_region_3 is not None:
            self.reference_region_3_filename = os.path.basename(self.reference_region_3)
            self.reference_region_3 = self.load_nii_resliced(self.reference_region_3, mask=True)

        if self.affine_nu_img is not None:
            self.affine_nu_img_filename = os.path.basename(self.affine_nu_img)
            self.affine_nu_img = self.load_nii(self.affine_nu_img)

        if self.affine_suvr_img is not None:
            self.affine_suvr_img_filename = os.path.basename(self.affine_suvr_img)
            self.affine_suvr_img = self.load_nii(self.affine_suvr_img)

        if self.warped_nu_img is not None:
            self.warped_nu_img_filename = os.path.basename(self.warped_nu_img)
            self.warped_nu_img = self.load_nii(self.warped_nu_img)

        if self.warped_suvr_img is not None:
            self.warped_suvr_img_filename = os.path.basename(self.warped_suvr_img)
            self.warped_suvr_img = self.load_nii(self.warped_suvr_img)


    # ______________________________________________________________________________________________________________________________________________ #
    # The following functions generate the slices for the provided images. The slices are generated using the QCImageGenerator class from the plotter.py file.

    def suvr_img_slices(self):
        """
        This function generates the suvr_img.

        Returns
        -------
        numpy.ndarray
            The suvr_img slices.

        """
        return QCImageGenerator(
            underlay_img=self.suvr_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            width_padding=50, height_padding=100).generate_qc_images()
        
    def nu_img_slices(self):
        """
        Usage
        -----
        nu_img_slices()
        """
        return QCImageGenerator(
            underlay_img=self.nu_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            crop_neck=self.aparc_img if self.crop_neck else None).generate_qc_images()
    
    def nu_img_lines(self):
        """
        This function generates the lines representing the axial, sagittal, and coronal slices on the nu_img.
        Returns
        -------
        numpy.ndarray : The arrays with the lines representing the slices.
        """
        return QCImageGenerator(
            underlay_img=self.nu_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            crop_neck=self.aparc_img if self.crop_neck else None).generate_lines()
    
    def mri_based_suvr_img_slices(self):
        """
        This funnction generates the SUVR image slices with the MRI as the underlay image.
        
        Returns
        -------
        numpy.ndarray : The suvr image slices
        """
        _, mri_based_suvr_img_slices = QCImageGenerator(
            underlay_img=self.nu_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            overlay_img=self.suvr_img,
            crop_neck=self.aparc_img if self.crop_neck else None).generate_qc_images()
        
        return mri_based_suvr_img_slices
    
    def reference_region_slices(self, reference_region):
        """
        This function generates the slices for the reference region to be overlayed on the nu_img. 
        The function uses nu_img as the underlay image as a base to determine the bounding box for the slices so that it is accurately cropped.

        Returns
        -------
        numpy.ndarray : The Reference Region Slices.
        """
        _, reference_region_slices = QCImageGenerator(
            underlay_img=self.nu_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            overlay_img=reference_region,
            crop_neck=self.aparc_img if self.crop_neck else None).generate_qc_images()
        return reference_region_slices
        
    def c1_image_slices(self):
        """
        1. This function generates the slices for the c1(spm segmentation of the Tissue Probability Maps) i.e. the Grey matter (containing high densities of unmyelinated neurons).
        2. The c1 has been thresholded between 0.3 and 1 to be overlaid on the nu_img.
        3. The function uses nu_img as the underlay image as a base to determine the bounding box for the slices so that it is accurately cropped.
        
        Returns
        -------
        numpy.ndarray : The c1_img slices.
        """
        _,c1_img_slices = QCImageGenerator(
            underlay_img=self.nu_img,
            overlay_img=self.c1_img,
            select_axial_slices=self.axial_slices,
            select_sagittal_slices=self.sagittal_slices,
            select_coronal_slices=self.coronal_slices,
            #mask_lower_threshold=0.8, mask_upper_threshold= self.c1_img.max(),
            crop_neck=self.aparc_img if self.crop_neck else None).generate_qc_images()
        
        return c1_img_slices
    
    def affine_suvr_img_slices(self):
        """
        """
        
        _, affine_suvr_img_slices = QCImageGenerator(
            underlay_img=self.affine_nu_img,
            select_axial_slices= template_axial_slices,
            select_sagittal_slices= template_sagittal_slices,
            select_coronal_slices= template_coronal_slices,
            overlay_img=self.affine_suvr_img,
            width_padding=3).generate_qc_images()
        
        _, tpm_img_slices = QCImageGenerator(
            underlay_img=self.affine_nu_img,
            select_axial_slices= template_axial_slices,
            select_sagittal_slices= template_sagittal_slices,
            select_coronal_slices= template_coronal_slices,
            overlay_img=tpm_image,
            mask_lower_threshold=0.3, mask_upper_threshold=1,
            width_padding=3).generate_qc_images()
        
        return affine_suvr_img_slices, tpm_img_slices
    
    def warped_suvr_img_slices(self):
        """
        """


        _, warped_suvr_image_slices = QCImageGenerator(
            underlay_img=self.warped_nu_img,
            overlay_img=self.warped_suvr_img,
            select_axial_slices= template_axial_slices,
            select_sagittal_slices= template_sagittal_slices,
            select_coronal_slices= template_coronal_slices,
            width_padding=3).generate_qc_images()
        
        _, tpm_img_slices = QCImageGenerator(
            underlay_img=self.warped_nu_img,
            overlay_img=tpm_image,
            select_axial_slices= template_axial_slices,
            select_sagittal_slices= template_sagittal_slices,
            select_coronal_slices= template_coronal_slices,            
            mask_lower_threshold=0.3, mask_upper_threshold=1,
            width_padding=3).generate_qc_images()
        
        return warped_suvr_image_slices, tpm_img_slices

    #  _____________________________________________________________ Defining the plotting functions _____________________________________________________________ #
    # The following functions have been defined to generate the plots for the MRI slices. 
    # The functions take the axes as input and plot the images on the axes.
    
    def plot_suvr_slices(self, axes):
        """
        This function plots the suvr_img slices.
        """
        sns.heatmap(self.suvr_img_slices(), cmap=cmap_turbo, vmin=0.1, vmax =pet_vmax , cbar=False, ax=axes)
        axes.text(10, 30, 'L', fontsize=10, color='white')
        axes.text(150, 30, 'R', fontsize=10, color='white')
        axes.set_title(f" {self.suvr_img_filename}", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')
        
    def plot_mri_slices(self, axes):
        """
        This functions plots only the nu_img slices
        """
        # Plotting the nu_img slices
        sns.heatmap(self.nu_img_slices(), cmap='gray', vmax=mri_vmax, cbar=False, ax=axes) 
        # Plotting the lines representing the slices
        sns.heatmap(self.nu_img_lines(), cmap=cmap_yellow2, vmin=0.5, cbar=False, ax=axes)
        axes.text(10, 30, 'L', fontsize=10, color='white')
        axes.text(150, 30, 'R', fontsize=10, color='white')
        axes.set_title(f"{self.nu_img_filename}", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')
    
    def plot_mri_based_suvr_img_slices(self, axes):
        """
        This functions plots the mri_based_suvr_img_slices
        """
        
        sns.heatmap(self.mri_based_suvr_img_slices(), cmap=cmap_turbo, vmin=0.1, vmax =pet_vmax , cbar=False, ax=axes)
        axes.set_title(f" {self.suvr_img_filename}", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')

    def plot_mri_based_suvr_img_slices_overlaid_on_mri(self, axes):
        """
        This functions plots the mri_based_suvr_img_slices
        """
        
        sns.heatmap(self.nu_img_slices(), cmap='gray', vmax = mri_vmax, cbar=False, ax=axes)
        sns.heatmap(self.mri_based_suvr_img_slices(), cmap=cmap_turbo, vmin=0.1, vmax =pet_vmax , cbar=False, alpha = 0.6, mask=self.mri_based_suvr_img_slices()==0, ax=axes)
        axes.set_title(f" Underlay: {self.nu_img_filename} \n Overlay: {self.suvr_img_filename}", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')

    def plot_c1_img_slices(self, axes):

        sns.heatmap(self.mri_based_suvr_img_slices(), cmap="gray", vmin=0.1, vmax=pet_vmax, cbar=False, ax=axes)
        sns.heatmap(self.c1_image_slices(), cmap=cmap_red, vmax=1, mask=(self.c1_image_slices())==0, cbar=False, ax=axes)
        axes.set_title(f" Underlay: {self.suvr_img_filename} \n Overlay: {self.c1_img_filename} (voxels > 0.3)", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')

    ### !! Important !! -- Original function to plot the reference regions
    """
    def plot_reference_region_slices(self, axes):
        cbl_ref_slices = self.reference_region_slices(self.reference_region_1)
        eroded_subcortwm_ref_slices = self.reference_region_slices(self.reference_region_2)
        brainstem_ref_slices = self.reference_region_slices(self.reference_region_3)

        sns.heatmap(self.nu_img_slices(), cmap='gray', vmax=mri_vmax, cbar=False, ax=axes)
        
        sns.heatmap(brainstem_ref_slices, cmap=cmap_blue2, vmin=0.1, alpha=0.3, cbar=False, mask=brainstem_ref_slices==0, ax=axes)
        sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(brainstem_ref_slices, lower_threshold = 0.1, upper_threshold=brainstem_ref_slices.max())), cbar = False, cmap = cmap_blue2, vmin = 0.1, ax = axes)

        sns.heatmap(cbl_ref_slices, cmap=cmap_orange, vmin=0.1, alpha=0.3, cbar=False, mask=cbl_ref_slices==0, ax=axes)
        sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(cbl_ref_slices, lower_threshold = 0.1, upper_threshold=cbl_ref_slices.max())), cbar = False, cmap = cmap_orange, vmin = 0.1, ax = axes)

        sns.heatmap(eroded_subcortwm_ref_slices, cmap=cmap_green, vmin=0.1, alpha=0.3, cbar=False, mask=eroded_subcortwm_ref_slices==0, ax=axes)
        sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(eroded_subcortwm_ref_slices, lower_threshold = 0.1, upper_threshold=eroded_subcortwm_ref_slices.max())), cbar = False, cmap = cmap_green, vmin = 0.1, ax = axes)

        axes.set_title(f" Underlay: {self.nu_img_filename} \n Overlay: {self.reference_region_1_filename}, {self.reference_region_2_filename}, {self.reference_region_3_filename}", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')
    """
    def plot_reference_region_slices(self, axes):
        
        sns.heatmap(self.nu_img_slices(), cmap='gray', vmax=mri_vmax, cbar=False, ax=axes)
        
        if self.reference_region_1 is not None and self.reference_region_2 is None and self.reference_region_3 is None:
            cbl_ref_slices = self.reference_region_slices(self.reference_region_1)
            sns.heatmap(cbl_ref_slices, cmap=cmap_orange, vmin=0.1, alpha=0.3, cbar=False, mask=cbl_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(cbl_ref_slices, lower_threshold = 0.1, upper_threshold=cbl_ref_slices.max())), cbar = False, cmap = cmap_orange, vmin = 0.1, ax = axes)
            axes.set_title(f" Underlay: {self.nu_img_filename} \n Overlay: {self.reference_region_1_filename}", fontsize=10, color='white', loc='left')
            axes.axis('off')
            axes.set_aspect('equal')

        elif self.reference_region_1 is not None and self.reference_region_2 is not None and self.reference_region_3 is None:
            cbl_ref_slices = self.reference_region_slices(self.reference_region_1)
            eroded_subcortwm_ref_slices = self.reference_region_slices(self.reference_region_2)
            sns.heatmap(cbl_ref_slices, cmap=cmap_orange, vmin=0.1, alpha=0.3, cbar=False, mask=cbl_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(cbl_ref_slices, lower_threshold = 0.1, upper_threshold=cbl_ref_slices.max())), cbar = False, cmap = cmap_orange, vmin = 0.1, ax = axes)

            sns.heatmap(eroded_subcortwm_ref_slices, cmap=cmap_green, vmin=0.1, alpha=0.3, cbar=False, mask=eroded_subcortwm_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(eroded_subcortwm_ref_slices, lower_threshold = 0.1, upper_threshold=eroded_subcortwm_ref_slices.max())), cbar = False, cmap = cmap_green, vmin = 0.1, ax = axes)
            axes.set_title(f" Underlay: {self.nu_img_filename} \n Overlay: {self.reference_region_1_filename}, {self.reference_region_2_filename}", fontsize=10, color='white', loc='left')
            axes.axis('off')
            axes.set_aspect('equal')

        elif self.reference_region_1 is not None and self.reference_region_2 is not None and self.reference_region_3 is not None:
            cbl_ref_slices = self.reference_region_slices(self.reference_region_1)
            eroded_subcortwm_ref_slices = self.reference_region_slices(self.reference_region_2)
            brainstem_ref_slices = self.reference_region_slices(self.reference_region_3)

            sns.heatmap(brainstem_ref_slices, cmap=cmap_blue2, vmin=0.1, alpha=0.3, cbar=False, mask=brainstem_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(brainstem_ref_slices, lower_threshold = 0.1, upper_threshold=brainstem_ref_slices.max())), cbar = False, cmap = cmap_blue2, vmin = 0.1, ax = axes)

            sns.heatmap(cbl_ref_slices, cmap=cmap_orange, vmin=0.1, alpha=0.3, cbar=False, mask=cbl_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(cbl_ref_slices, lower_threshold = 0.1, upper_threshold=cbl_ref_slices.max())), cbar = False, cmap = cmap_orange, vmin = 0.1, ax = axes)

            sns.heatmap(eroded_subcortwm_ref_slices, cmap=cmap_green, vmin=0.1, alpha=0.3, cbar=False, mask=eroded_subcortwm_ref_slices==0, ax=axes)
            sns.heatmap(ImageProcessor.contour_image(ImageProcessor.mask_image(eroded_subcortwm_ref_slices, lower_threshold = 0.1, upper_threshold=eroded_subcortwm_ref_slices.max())), cbar = False, cmap = cmap_green, vmin = 0.1, ax = axes)

            axes.set_title(f" Underlay: {self.nu_img_filename} \n Overlay: {self.reference_region_1_filename}, {self.reference_region_2_filename}, {self.reference_region_3_filename}", fontsize=10, color='white', loc='left')
            axes.axis('off')
            axes.set_aspect('equal')

    def plot_affine_suvr_img_slices(self, axes):
        affine_suvr_slices, afftpm_slices = self.affine_suvr_img_slices()
        ctr_arr = ImageProcessor.contour_image(afftpm_slices)
        sns.heatmap(affine_suvr_slices, cmap="gray", vmin = 0.1, vmax=pet_vmax, cbar=False, ax=axes)
        sns.heatmap(afftpm_slices, cmap=cmap_pink, vmin=0.1, alpha=0.1, mask=afftpm_slices==0, cbar=False, ax=axes)
        sns.heatmap(ctr_arr, cmap=cmap_pink, vmin=0.1, mask=ctr_arr==0, cbar=False,  ax=axes)

        axes.set_title(f" Underlay: {self.affine_suvr_img_filename} \n Overlay: TPM.nii (c1, voxels > 0.3)", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')


    def plot_warped_suvr_img_slices(self, axes):
        warped_suvr_slices, wsuvr_tpm_slices = self.warped_suvr_img_slices()
        ctr_arr = ImageProcessor.contour_image(wsuvr_tpm_slices)
        sns.heatmap(warped_suvr_slices, cmap="gray", vmin = 0.1, vmax=pet_vmax, cbar=False, ax=axes)
        sns.heatmap(wsuvr_tpm_slices, cmap=cmap_pink, vmin=0.1, mask=wsuvr_tpm_slices==0, alpha=0.1, cbar=False, ax=axes)
        sns.heatmap(ctr_arr, cmap=cmap_pink, vmin=0.1, mask=ctr_arr==0, cbar=False, ax=axes)

        axes.set_title(f" Underlay: {self.warped_suvr_img_filename} \n Overlay: TPM.nii (c1, voxels > 0.3)", fontsize=10, color='white', loc='left')
        axes.axis('off')
        axes.set_aspect('equal')

    #  _____________________________________________________________ Plotting the slices _____________________________________________________________ #
    def plot_slices(self,output_path):
        """
        Plotting the slices that are in arrays using seaborn' heatmap function.
        1. The suvr_img_slices will be plotted as default image.
        2. If the nu_img is provided, the first row will be the nu_img_slices, the second row will be mri_based_suvr_img_slices, and the third row will be the mri_based_suvr_img_slices with the nu_img as underlay.
        3. If the c1_img is provided, the fourth row will be the c1_img_slices.
        4. If the reference_region_1 is provided, the fifth row will be the reference_region_1_slices overlayed.
        5. If the reference_region_2 is provided, the fifth row will be the reference_region_1_slices + reference_region_2_slices overlayed.
        6. If the affine_suvr_img and warped_suvr_img are provided, the sixth row will be the affine_suvr_img_slices + TPM and the seventh row will be the warped_suvr_img_slices + TPM
        """

        # Load the images   
        self.load_images()
        plt.figure(facecolor='black')
        
        # If only the suvr_img is provided
        if self.suvr_img is not None and self.nu_img is None and self.c1_img is None and self.reference_region_1 is None and self.reference_region_2 is None and self.affine_suvr_img is None and self.warped_suvr_img is None:
            fig, axes = plt.subplots(1, 1, figsize=(17, 2.5))
            self.plot_suvr_slices(axes)
            add_colorbar(fig, axes[2], cmap='turbo', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)
            

        # If the suvr_img and nu_img are provided
        elif self.suvr_img is not None and self.nu_img is not None and self.c1_img is None and self.reference_region_1 is None and self.reference_region_2 is None and self.affine_suvr_img is None and self.warped_suvr_img is None:
            fig, axes = plt.subplots(3, 1, figsize=(17, 6))

            self.plot_mri_slices(axes[0])
            self.plot_mri_based_suvr_img_slices(axes[1])
            self.plot_mri_based_suvr_img_slices_overlaid_on_mri(axes[2])
            add_colorbar(fig, axes[2], cmap='turbo', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)
            

        # If the suvr_img, nu_img, and c1_img are provided
        elif self.suvr_img is not None and self.nu_img is not None and self.c1_img is not None and self.reference_region_1 is None and self.reference_region_2 is None and self.affine_suvr_img is None and self.warped_suvr_img is None:
            fig, axes = plt.subplots(4, 1, figsize=(17, 8))

            self.plot_mri_slices(axes[0])
            self.plot_mri_based_suvr_img_slices(axes[1])
            self.plot_mri_based_suvr_img_slices_overlaid_on_mri(axes[2])
            add_colorbar(fig, axes[2], cmap='turbo', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)
            self.plot_c1_img_slices(axes[3])
            add_colorbar(fig, axes[3], cmap='gray', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)

        # if all the images are provided

        elif self.suvr_img is not None and self.nu_img is not None and self.c1_img is not None and self.reference_region_1 is not None and self.reference_region_2 is not None and self.affine_suvr_img is not None and self.warped_suvr_img is not None:

            fig, axes = plt.subplots(7, 1, figsize=(17, 15))

            # Adjust the spacing between subplots
            fig.subplots_adjust(hspace=0.4)

            self.plot_mri_slices(axes[0])
            self.plot_mri_based_suvr_img_slices(axes[1])
            self.plot_mri_based_suvr_img_slices_overlaid_on_mri(axes[2])
            add_colorbar(fig, axes[2], cmap='turbo', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)
            self.plot_c1_img_slices(axes[3])
            add_colorbar(fig, axes[3], cmap='gray', vmin=0, vmax=2.5, ticks=[0, 2.5], orientation='horizontal', cbar_width=0.13, cbar_height=0.01, cbar_x=0.76)
            self.plot_reference_region_slices(axes[4])
            self.plot_affine_suvr_img_slices(axes[5])
            self.plot_warped_suvr_img_slices(axes[6])

        fig.patch.set_facecolor('black')
        #plt.subplots_adjust(top=1, wspace=0, hspace=0.3)
        #plt.tight_layout()
        # Removing the 'r' string from the self.basename
        file_name = self.basename.replace('r', '')
        # Adding '_qc' to the file_name at the end
        file_name = file_name + '_qc'
        plt.savefig(os.path.join(output_path, file_name + '.png'), facecolor='black', bbox_inches='tight', dpi=500)
