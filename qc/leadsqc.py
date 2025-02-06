#!/usr/bin/env python

"""
$ leadsqc.py /home/mac/pmaiti/Desktop/leads_qc/mimic_processed_daniel/LDS1770688/FTP_2024-06-04
"""
import os
import sys
import argparse
import warnings
import glob

import nibabel as nib
from nibabel.orientations import io_orientation
from nilearn.image import resample_to_img

# Importing the necessary classes from the rablabqc package
rablab_pkg_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(rablab_pkg_path)

from mri_slices import MRIQCplots
from fbb_slices import FBBQCplots
from ftp_slices import FTPQCplots
from fdg_slices import FDGQCplots
from slice_selector import SliceSelector

warnings.filterwarnings("ignore", category=RuntimeWarning, 
                        message="NaNs or infinite values are present in the data passed to resample")

def build_parser():
    p = argparse.ArgumentParser(
        description="Description : Python Script to generate Quality Control Images\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "path",
        type=str,
        help="Path to the modality directory containing the MRI and PET files",
    )
    return p

rtpm_path = os.path.join(
    rablab_pkg_path, "reslice", "rT1.nii")  # Resliced T1 file to 1mm isotropic resolution


def load_nii_resliced(path, orientation="LAS", mask=False):
    """
    Load nifti image with specified orientation
    """

    id = path.split("/")[-1].split(".")[0]
    out_folder = os.path.dirname(path)
    resliced_image_path = os.path.join(out_folder, "qc" + id + ".nii")

    if not os.path.exists(resliced_image_path):

        template_img = nib.load(rtpm_path)
        source_img = nib.load(os.path.join(out_folder, id + ".nii"))

        if mask:
            resliced_img = resample_to_img(source_img, template_img, interpolation='nearest', copy_header=True, force_resample=True)
        else:
            resliced_img = resample_to_img(source_img, template_img, interpolation='continuous', copy_header=True, force_resample=True)
        nib.save(resliced_img, resliced_image_path)
    else:
        print(
            "Resliced image already exists for ", id, ". Loading the resliced image..."
        )

    img = nib.load(resliced_image_path)
    img_ornt = io_orientation(img.affine)
    img = img.as_reoriented(img_ornt)
    return img.get_fdata()


# _______________________________________________ Function to Process Images _______________________________________________
def process_qc_images(results, modality):

    id = (results.path.rstrip("/")).split("/")[-2]
    folder = results.path.rstrip("/").split("/")[-1]
    folder_path = results.path

    print(
        f"Processing {modality} for ID : {id}, Folder : {folder}, located \n {folder_path}"
    )

    # ______________________ MRI QC ______________________
    if modality == "MRI":
        mri_date = folder.split("_")[-1]

        # Searching for the relevant files
        nu_img = os.path.join(folder_path, f"{id}_MRI-T1_{mri_date}_nu.nii")
        aparc_aseg_img = os.path.join(
            folder_path, f"{id}_MRI-T1_{mri_date}_aparc+aseg.nii"
        )
        c1_img = os.path.join(folder_path, f"c1{id}_MRI-T1_{mri_date}_nu.nii")
        wnu_img = os.path.join(folder_path, f"w{id}_MRI-T1_{mri_date}_nu.nii")
        affinenu_img = os.path.join(folder_path, f"a{id}_MRI-T1_{mri_date}_nu.nii")

        # ______________________ Generating MRI QC Images ______________________
        if all(
            map(os.path.exists, [nu_img, aparc_aseg_img, c1_img, wnu_img, affinenu_img])
        ):
            print("All the files are present, proceeding with slice selection")
            select_axial_slices, select_coronal_slices, select_sagittal_slices = (
                SliceSelector(
                    load_nii_resliced(aparc_aseg_img, mask=True)
                ).select_leads_slices()
            )
            print(
                "Selected slices for MRI QC : ",
                select_axial_slices,
                select_coronal_slices,
                select_sagittal_slices,
            )

            MRIQCplots(
                nu_img=nu_img,
                aparc_img=aparc_aseg_img,
                c1_img=c1_img,
                affine_nu_img=affinenu_img,
                warped_nu_img=wnu_img,
                axial_slices=select_axial_slices,
                coronal_slices=select_coronal_slices,
                sagittal_slices=select_sagittal_slices,
            ).plot_slices(results.path)

            print(" -- MRI QC Image Generated -- ")

            # Delete qc temp files
            files_to_delete = glob.glob(os.path.join(folder_path, "qc*.nii"))
            for file in files_to_delete:
                os.remove(file)
                print(f"Deleted: {file}")
        else:
            print("Some files are missing")
            return

    # ______________________ FBB QC ______________________
    if modality == "FBB":
        fbb_date = folder.split("_")[-1]

        print("Gathering relevant FBB files")
        fbb_suvr_img = os.path.join(folder_path, f"r{id}_FBB_{fbb_date}_suvr-wcbl.nii")
        fbb_affine_suvr_img = os.path.join(
            folder_path, f"ar{id}_FBB_{fbb_date}_suvr-wcbl.nii"
        )
        fbb_warped_suvr_img = os.path.join(
            folder_path, f"wr{id}_FBB_{fbb_date}_suvr-wcbl.nii"
        )

        print("Gathering the relevant MRI files")
        fbb_related_mri = os.readlink(os.path.join(folder_path, "mri"))
        fbb_related_mri_date = fbb_related_mri.rstrip("/").split("_")[-1]

        fbb_nu_img = os.path.join(
            fbb_related_mri, f"{id}_MRI-T1_{fbb_related_mri_date}_nu.nii"
        )
        fbb_aparc_aseg_img = os.path.join(
            fbb_related_mri, f"{id}_MRI-T1_{fbb_related_mri_date}_aparc+aseg.nii"
        )
        fbb_c1_img = os.path.join(
            fbb_related_mri, f"c1{id}_MRI-T1_{fbb_related_mri_date}_nu.nii"
        )
        fbb_wnu_img = os.path.join(
            fbb_related_mri, f"w{id}_MRI-T1_{fbb_related_mri_date}_nu.nii"
        )
        fbb_affinenu_img = os.path.join(
            fbb_related_mri, f"a{id}_MRI-T1_{fbb_related_mri_date}_nu.nii"
        )

        print("Gathering the Reference Region Masks from the MRI")

        wcbl_reference_mask = os.path.join(
            fbb_related_mri, f"{id}_MRI-T1_{fbb_related_mri_date}_mask-wcbl.nii"
        )
        eroded_subcortwm_reference_mask = os.path.join(
            fbb_related_mri,
            f"{id}_MRI-T1_{fbb_related_mri_date}_mask-eroded-subcortwm.nii",
        )
        brainstem_reference_mask = os.path.join(
            fbb_related_mri, f"{id}_MRI-T1_{fbb_related_mri_date}_mask-brainstem.nii"
        )

        if all(
            map(
                os.path.exists,
                [
                    fbb_suvr_img,
                    fbb_affine_suvr_img,
                    fbb_warped_suvr_img,
                    fbb_nu_img,
                    fbb_aparc_aseg_img,
                    fbb_c1_img,
                    fbb_wnu_img,
                    fbb_affinenu_img,
                    wcbl_reference_mask,
                    brainstem_reference_mask,
                    eroded_subcortwm_reference_mask,
                ],
            )
        ):

            print("All the files are present, proceeding with slice selection")
            fbb_axial_slices, fbb_coronal_slices, fbb_sagittal_slices = SliceSelector(
                load_nii_resliced(fbb_aparc_aseg_img, mask=True)
            ).select_leads_slices()

            FBBQCplots(
                suvr_img=fbb_suvr_img,
                affine_suvr_img=fbb_affine_suvr_img,
                warped_suvr_img=fbb_warped_suvr_img,
                nu_img=fbb_nu_img,
                aparc_img=fbb_aparc_aseg_img,
                c1_img=fbb_c1_img,
                affine_nu_img=fbb_affinenu_img,
                warped_nu_img=fbb_wnu_img,
                reference_region_1=wcbl_reference_mask,
                reference_region_2=eroded_subcortwm_reference_mask,
                reference_region_3=brainstem_reference_mask,
                axial_slices=fbb_axial_slices,
                coronal_slices=fbb_coronal_slices,
                sagittal_slices=fbb_sagittal_slices,
            ).plot_slices(results.path)

            print(" -- FBB QC Image Generated -- ")

            # Delete temp qc files
            folders = [folder_path, fbb_related_mri]
            for folder in folders:
                for file in glob.glob(os.path.join(folder, "qc*.nii")):
                    os.remove(file)
                    print(f"Deleted: {file}")
                
        else:
            print("Some files are missing")
            return

    # ______________________ FTP QC ______________________
    if modality == "FTP":
        ftp_date = folder.split("_")[-1]

        print("Gathering relevant FTP files")
        ftp_suvr_img = os.path.join(
            folder_path, f"r{id}_FTP_{ftp_date}_suvr-infcblgm.nii"
        )
        affine_suvr_img = os.path.join(
            folder_path, f"ar{id}_FTP_{ftp_date}_suvr-infcblgm.nii"
        )
        warped_suvr_img = os.path.join(
            folder_path, f"wr{id}_FTP_{ftp_date}_suvr-infcblgm.nii"
        )

        print("Gathering the relevant MRI files")
        ftp_related_mri = os.readlink(os.path.join(folder_path, "mri"))
        ftp_related_mri_date = ftp_related_mri.rstrip("/").split("_")[-1]

        ftp_nu_img = os.path.join(
            ftp_related_mri, f"{id}_MRI-T1_{ftp_related_mri_date}_nu.nii"
        )
        ftp_aparc_aseg_img = os.path.join(
            ftp_related_mri, f"{id}_MRI-T1_{ftp_related_mri_date}_aparc+aseg.nii"
        )
        ftp_c1_img = os.path.join(
            ftp_related_mri, f"c1{id}_MRI-T1_{ftp_related_mri_date}_nu.nii"
        )
        ftp_wnu_img = os.path.join(
            ftp_related_mri, f"w{id}_MRI-T1_{ftp_related_mri_date}_nu.nii"
        )
        ftp_affinenu_img = os.path.join(
            ftp_related_mri, f"a{id}_MRI-T1_{ftp_related_mri_date}_nu.nii"
        )

        print("Gathering the Reference Region Masks from the MRI")
        infcblgm_reference_mask = os.path.join(
            ftp_related_mri, f"{id}_MRI-T1_{ftp_related_mri_date}_mask-infcblgm.nii"
        )
        ftp_eroded_subcortwm_reference_mask = os.path.join(
            ftp_related_mri,
            f"{id}_MRI-T1_{ftp_related_mri_date}_mask-eroded-subcortwm.nii",
        )

        if all(
            map(
                os.path.exists,
                [
                    ftp_suvr_img,
                    affine_suvr_img,
                    warped_suvr_img,
                    ftp_nu_img,
                    ftp_aparc_aseg_img,
                    ftp_c1_img,
                    ftp_wnu_img,
                    ftp_affinenu_img,
                    infcblgm_reference_mask,
                    ftp_eroded_subcortwm_reference_mask,
                ],
            )
        ):

            print("All the files are present, proceeding with slice selection")
            ftp_axial_slices, ftp_coronal_slices, ftp_sagittal_slices = SliceSelector(
                load_nii_resliced(ftp_aparc_aseg_img, mask=True)
            ).select_leads_slices()

            FTPQCplots(
                suvr_img=ftp_suvr_img,
                affine_suvr_img=affine_suvr_img,
                warped_suvr_img=warped_suvr_img,
                nu_img=ftp_nu_img,
                aparc_img=ftp_aparc_aseg_img,
                c1_img=ftp_c1_img,
                affine_nu_img=ftp_affinenu_img,
                warped_nu_img=ftp_wnu_img,
                reference_region_1=infcblgm_reference_mask,
                reference_region_2=ftp_eroded_subcortwm_reference_mask,
                axial_slices=ftp_axial_slices,
                coronal_slices=ftp_coronal_slices,
                sagittal_slices=ftp_sagittal_slices,
            ).plot_slices(results.path)

            print(" -- FTP QC Image Generated -- ")
            
            # Delete temp qc files
            folders = [folder_path, ftp_related_mri]
            for folder in folders:
                for file in glob.glob(os.path.join(folder, "qc*.nii")):
                    os.remove(file)
                    print(f"Deleted: {file}")
        else:
            print("Some files are missing")
            return

    # ______________________ FDG QC ______________________
    if modality == "FDG":
        fdg_date = folder.split("_")[-1]

        print("Gathering relevant FDG files")
        fdg_suvr_img = os.path.join(folder_path, f"r{id}_FDG_{fdg_date}_suvr-pons.nii")
        fdg_affine_suvr_img = os.path.join(
            folder_path, f"ar{id}_FDG_{fdg_date}_suvr-pons.nii"
        )
        fdg_warped_suvr_img = os.path.join(
            folder_path, f"wr{id}_FDG_{fdg_date}_suvr-pons.nii"
        )

        print("Gathering the relevant MRI files")
        fdg_related_mri = os.readlink(os.path.join(folder_path, "mri"))
        fdg_related_mri_date = fdg_related_mri.rstrip("/").split("_")[-1]

        fdg_nu_img = os.path.join(
            fdg_related_mri, f"{id}_MRI-T1_{fdg_related_mri_date}_nu.nii"
        )
        fdg_aparc_aseg_img = os.path.join(
            fdg_related_mri, f"{id}_MRI-T1_{fdg_related_mri_date}_aparc+aseg.nii"
        )
        fdg_c1_img = os.path.join(
            fdg_related_mri, f"c1{id}_MRI-T1_{fdg_related_mri_date}_nu.nii"
        )
        fdg_wnu_img = os.path.join(
            fdg_related_mri, f"w{id}_MRI-T1_{fdg_related_mri_date}_nu.nii"
        )
        fdg_affinenu_img = os.path.join(
            fdg_related_mri, f"a{id}_MRI-T1_{fdg_related_mri_date}_nu.nii"
        )

        print("Gathering the Reference Region Masks from the MRI")
        pons_reference_mask = os.path.join(
            fdg_related_mri, f"{id}_MRI-T1_{fdg_related_mri_date}_mask-pons.nii"
        )

        if all(
            map(
                os.path.exists,
                [
                    fdg_suvr_img,
                    fdg_affine_suvr_img,
                    fdg_affinenu_img,
                    fdg_warped_suvr_img,
                    fdg_nu_img,
                    fdg_aparc_aseg_img,
                    fdg_c1_img,
                    fdg_wnu_img,
                    fdg_affinenu_img,
                    pons_reference_mask,
                ],
            )
        ):

            print("All the files are present, proceeding with slice selection")
            fdg_axial_slices, fdg_coronal_slices, fdg_sagittal_slices = SliceSelector(
                load_nii_resliced(fdg_aparc_aseg_img, mask=True)
            ).select_leads_slices()

            FDGQCplots(
                suvr_img=fdg_suvr_img,
                affine_suvr_img=fdg_affine_suvr_img,
                warped_suvr_img=fdg_warped_suvr_img,
                nu_img=fdg_nu_img,
                aparc_img=fdg_aparc_aseg_img,
                c1_img=fdg_c1_img,
                affine_nu_img=fdg_affinenu_img,
                warped_nu_img=fdg_wnu_img,
                reference_region_1=pons_reference_mask,
                axial_slices=fdg_axial_slices,
                coronal_slices=fdg_coronal_slices,
                sagittal_slices=fdg_sagittal_slices,
            ).plot_slices(results.path)

            print(" -- FDG QC Image Generated -- ")

            # Delete temp qc files
            folders = [folder_path, fdg_related_mri]
            for folder in folders:
                for file in glob.glob(os.path.join(folder, "qc*.nii")):
                    os.remove(file)
                    print(f"Deleted: {file}")
        else:
            print("Some files are missing")
            return


def main():
    parser = build_parser()
    results = parser.parse_args()
    results.path = os.path.abspath(results.path)

    if not os.path.exists(results.path):
        print("Error: Input path does not exist.")
        return

    id = (results.path).split("/")[-2]
    print("Processing ID : ", id)
    print("-----------")

    # Extracting the modality from the path name
    modality = os.path.basename(results.path.rstrip("/")).split("_")[0]

    print("Modality : ", modality)

    if modality == "MRI-T1":
        process_qc_images(results, "MRI")
    elif modality == "FBB":
        process_qc_images(results, "FBB")
    elif modality == "FTP":
        process_qc_images(results, "FTP")
    elif modality == "FDG":
        process_qc_images(results, "FDG")
    else:
        raise ValueError(
            f"Modality {modality} not recognized! Cannot generate QC image"
        )


if __name__ == "__main__":
    main()
