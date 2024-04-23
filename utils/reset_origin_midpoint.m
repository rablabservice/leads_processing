function reset_origin_midpoint(nii_file)
    % Remove any trailing whitespace from the file path
    nii_file = deblank(nii_file);

    % Load the image volume using SPM's spm_vol function
    st.vol = spm_vol(nii_file);

    % Compute the inverse of the affine transformation matrix
    vs = st.vol.mat \ eye(4);

    % Update the translation part of the matrix to the midpoint of the image dimensions
    vs(1:3,4) = (st.vol.dim + 1) / 2;

    % Update the image's affine transformation space using the modified inverse matrix
    spm_get_space(st.vol.fname, inv(vs));
end
