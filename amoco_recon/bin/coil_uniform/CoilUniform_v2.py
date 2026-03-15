import numpy as np
import SimpleITK as sitk
from skimage import filters
import h5py
import argparse
import os


def n4_bias_correction_4dflow(vol: np.ndarray) -> np.ndarray:
    """
    Perform N4 bias field correction for 4D Flow MRI data.

    Args:
        vol (np.ndarray): 5D array with shape (Nx, Ny, Nz, Nt, Nvencs),
                          where Nx, Ny, Nz are spatial dimensions,
                          Nt is the temporal dimension,
                          and Nvencs represents velocity encoding directions.

    Returns:
        np.ndarray: Bias-corrected 5D array with the same shape as the input.
    """

    if vol.ndim != 5:
        raise ValueError(
            f"Input data must have shape (Nx, Ny, Nz, Nt, Nvencs), but got {vol.shape}"
        )

    # 1. Reorder dimensions if necessary
    # vol_reordered = np.transpose(vol, (4, 3, 2, 1, 0))
    vol_reordered = vol

    # 2. Compute average magnitude image
    mag_avg = np.mean(
        np.sqrt(np.sum(np.square(np.abs(vol_reordered)), axis=0)),
        axis=0
    )

    # 3. Convert to SimpleITK image
    raw_img_sitk = sitk.GetImageFromArray(mag_avg)

    # 4. Generate foreground mask using multi-Otsu thresholding
    fg_mask = sitk.GetImageFromArray(
        (mag_avg > filters.threshold_multiotsu(mag_avg, classes=5)[0]).astype(np.uint8)
    )

    # 5. Perform N4 bias field correction
    bias_corrector = sitk.N4BiasFieldCorrectionImageFilter()
    _ = bias_corrector.Execute(raw_img_sitk, fg_mask)
    log_bias_field = bias_corrector.GetLogBiasFieldAsImage(raw_img_sitk)

    # 6. Apply bias correction to the original data
    bias_field_array = sitk.GetArrayFromImage(sitk.Exp(log_bias_field))
    bias_field_array = bias_field_array[np.newaxis, np.newaxis, :, :, :]

    # bias_field_array = np.transpose(bias_field_array, (4,3,2,1,0))

    corrected_vol = np.abs(vol) / bias_field_array * np.exp(1j * np.angle(vol))

    return corrected_vol


def main():
    parser = argparse.ArgumentParser(
        description="Perform N4 bias field correction for 4D Flow MRI data stored in an .h5 file."
    )
    parser.add_argument("input_filepath", type=str, help="Path to the input .h5 file.")
    args = parser.parse_args()

    input_path = args.input_filepath
    base, ext = os.path.splitext(input_path)
    output_path = base + "_corrected.h5"

    if not os.path.exists(input_path):
        print(f"Error: Input file does not exist - {input_path}")
        return

    print(f"Reading data from '{input_path}'...")

    try:
        with h5py.File(input_path, 'r') as f:
            if 'data_real' in f and 'data_imag' in f:
                data_real = f['data_real'][()]
                data_imag = f['data_imag'][()]
                vol_data = data_real + 1j * data_imag
            else:
                raise KeyError(
                    "The HDF5 file must contain 'data_real' and 'data_imag' datasets."
                )

        print(f"Data successfully loaded. Shape: {vol_data.shape}")

        # Perform bias field correction
        print("Starting bias field correction...")
        corrected_vol_data = n4_bias_correction_4dflow(vol_data)
        print("Bias field correction completed.")

        # Save corrected data (split into real and imaginary parts)
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('data_real', data=np.real(corrected_vol_data))
            f.create_dataset('data_imag', data=np.imag(corrected_vol_data))

        print(f"Corrected data successfully saved to '{output_path}'.")

    except Exception as e:
        print(f"An error occurred during processing: {e}")


if __name__ == '__main__':
    main()
