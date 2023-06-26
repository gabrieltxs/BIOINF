import ftplib
import gzip
import os
from tqdm import tqdm

def download_file_from_ftp(date, local_path):
    # FTP server credentials
    ftp_server = "ftp.ncbi.nlm.nih.gov"

    # Remote file to download
    remote_file_path = f"/genbank/daily-nc/nc{date}.flat.gz"

    # Local file path
    local_file_path = os.path.join(local_path, f"nc{date}.flat.gz")

    # Create an FTP session
    ftp = ftplib.FTP(ftp_server)
    ftp.login()

    # Retrieve file information
    file_size = ftp.size(remote_file_path)
    bytes_downloaded = 0

    # Download the remote file
    with open(local_file_path, "wb") as f:
        def callback(chunk):
            nonlocal bytes_downloaded
            bytes_downloaded += len(chunk)
            f.write(chunk)
            pbar.update(len(chunk))  # Update progress bar with the length of the downloaded chunk

        # Initialize tqdm with the total file size and set the unit to bytes
        with tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024, desc="Downloading") as pbar:
            # Retrieve the remote file and write it to the local file using the callback function
            ftp.retrbinary(f"RETR {remote_file_path}", callback)


    # Close the FTP session
    ftp.quit()

    # Unzip the file
    extract_folder = os.path.join(local_path, "flat_extract")
    os.makedirs(extract_folder, exist_ok=True)
    extract_file_path = os.path.join(extract_folder, f"nc{date}.flat")

    # Get the uncompressed file size
    with gzip.open(local_file_path, 'rb') as gz_file:
        uncompressed_size = gz_file.seek(0, os.SEEK_END)

    # Extract the file and track the progress
    with tqdm(total=uncompressed_size, unit='B', unit_scale=True, unit_divisor=1024, desc="Unzipping") as pbar:
        with gzip.open(local_file_path, 'rb') as gz_file:
            with open(extract_file_path, 'wb') as extracted_file:
                chunk_size = 1024
                while True:
                    chunk = gz_file.read(chunk_size)
                    if not chunk:
                        break
                    extracted_file.write(chunk)
                    pbar.update(len(chunk))  # Update unzip progress bar with the length of the extracted chunk


    print("File downloaded and extracted successfully!")

# Example usage
date = "0601"
pathx = "D:\\OneDrive\\Documentos\\OneDrive\\Documentos\\Pesquisa BIOINF - GENOMIC\\Python\\GitHub"


download_file_from_ftp(date, pathx)
