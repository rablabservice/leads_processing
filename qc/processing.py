import cv2
import numpy as np
import matplotlib.pyplot as plt

class ImageProcessor:
    def __init__(self):
        pass

    @staticmethod
    def brain_padding(underlay_image_array, overlay_image_array = None, width_padding=5, height_padding=5):
        """
        Adjusts padding around the brain by finding the brain boundary and adding or removing the padding around it.
        This can be done both on the PET and the MRI images
        Parameters
        ----------
        underlay_image_array : numpy.ndarray
            2D image array.
        
        overlay_image_array : numpy.ndarray, optional
            2D image array to be padded. If not provided, only the underlay image will be padded.
        pad_size : int, optional
            Size of the padding to be added.

        width_padding : int, optional
            Size of the padding to be added to the width.
            
        height_padding : int, optional
            Size of the padding to be added to the height.

        Returns
        -------
        numpy.ndarray
            Padded 2D image array.
        """
        image_array = underlay_image_array.copy()
        # Converting the image to uint8
        image_array_u8 = np.uint8(image_array)

        # Find the contours of the brain using opencv findContours
        contours, _ = cv2.findContours(image_array_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # Find the largest bounding box which is the brain bounding box. The bounding box with the largest area is the brain
        area = 0
        for contour in contours:
            x, y, w, h = cv2.boundingRect(contour)
            if w * h > area:
                area = w * h
                x_max, y_max, w_max, h_max = x, y, w, h
        
        # Note: Interchanging width and height padding to match the images after the orientation change
        x_max -= height_padding
        y_max -= width_padding
        w_max += 2 * height_padding
        h_max += 2 * width_padding

        # Check if the padding is within the image boundaries
        if x_max < 0:
            x_max = 0
        if y_max < 0:
            y_max = 0
        if x_max + w_max > image_array.shape[1]:
            w_max = image_array.shape[1] - x_max
        if y_max + h_max > image_array.shape[0]:
            h_max = image_array.shape[0] - y_max
            
        # Crop the image to the bounding box
        padded_image = image_array[y_max:y_max + h_max, x_max:x_max + w_max]
        padded_image = np.rot90(padded_image)

        # If an overlay image is provided, crop it to the same bounding box
        if overlay_image_array is not None:
            padded_overlay_image = overlay_image_array[y_max:y_max + h_max, x_max:x_max + w_max]
            padded_overlay_image = np.rot90(padded_overlay_image)
            return padded_image, padded_overlay_image

        return padded_image

    @staticmethod
    def mask_image(image, lower_threshold, upper_threshold):
        mask = np.zeros_like(image)
        # Set the pixel values = lower_threshold to 1 and pixel values = upper_threshold to 1, set everything else to 0
        mask[(image >= lower_threshold) & (image <= upper_threshold)] = 1
        # Convert the mask to absolute integers
        mask = mask.astype(np.uint8)
        
        return mask

    @staticmethod
    def contour_image(mask):
        # Create a copy of the input mask image
        mask_with_contours = mask.copy()

        # Find contours in the mask image
        contours, _ = cv2.findContours(mask_with_contours, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

        # Draw contours on the copied mask image
        cv2.drawContours(mask_with_contours, contours, -1, 255, 1)
        
        # change the value of the contours from 255 to 1 and set everything else to 0
        mask_with_contours[mask_with_contours <= 1] = 0
        mask_with_contours[mask_with_contours == 255 ] = 1
        return mask_with_contours

    @staticmethod
    def msk_ctr(image):
        return ImageProcessor.contour_image(ImageProcessor.mask_image(image, lower_threshold=0, upper_threshold=0))