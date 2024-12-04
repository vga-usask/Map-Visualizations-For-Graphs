import os
import shutil

def copy_dot_files(source_folders, destination_folder):
    # Loop through each source folder
    for source_folder in source_folders:
        source_folder_name = os.path.basename(source_folder)  # Extract source folder name
        
        # Find all subfolders in the destination folder
        for root, dirs, files in os.walk(destination_folder):
            for subfolder in dirs:
                if source_folder_name in subfolder:  # Check if subfolder contains the source folder name
                    # Build paths for source and destination
                    src_g_dot = os.path.join(source_folder, "G.dot")
                    src_h_dot = os.path.join(source_folder, "H.dot")
                    dest_folder = os.path.join(root, subfolder)
                    
                    # Copy files if they exist
                    if os.path.exists(src_g_dot):
                        shutil.copy(src_g_dot, dest_folder)
                        print(f"Copied {src_g_dot} to {dest_folder}")
                    else:
                        print(f"{src_g_dot} does not exist, skipping.")

                    if os.path.exists(src_h_dot):
                        shutil.copy(src_h_dot, dest_folder)
                        print(f"Copied {src_h_dot} to {dest_folder}")
                    else:
                        print(f"{src_h_dot} does not exist, skipping.")

# Define source and destination folders
source_folders = [
    "bookland", "musicland", "recipies", "tradeland", "universities",
    "dlpb_100", "dlpb_300", "dlpb_500", "dlpb_700", "dlpb_900"
]
destination_folder = "evaluation2"

# Call the function
copy_dot_files(source_folders, destination_folder)
