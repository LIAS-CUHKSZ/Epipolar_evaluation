import numpy as np
import os
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

def process_pose_file(file_path):
    # Read poses from file
    poses = np.loadtxt(file_path)
    
    # Convert poses to 3x4 transformation matrices
    transforms = []
    for pose in poses:
        T = pose.reshape(3, 4)
        transforms.append(T)
    
    # Calculate relative poses between consecutive frames
    relative_poses = []
    # Get the sequence number from the file path
    seq_num = Path(file_path).stem
    
    for i in tqdm(range(len(transforms)-1), desc=f"Processing {seq_num}", leave=False):
        # Current frame transformation (3x4)
        T_w_i = transforms[i]
        # Next frame transformation (3x4)
        T_w_ip1 = transforms[i+1]
        
        # Convert to 4x4 matrices for easier computation
        T_w_i_4x4 = np.vstack((T_w_i, [0, 0, 0, 1]))
        T_w_ip1_4x4 = np.vstack((T_w_ip1, [0, 0, 0, 1]))
        
        # Calculate relative transformation (i+1 to i)
        # T_i_ip1 = T_w_i^(-1) * T_w_ip1
        T_i_ip1 = np.linalg.inv(T_w_i_4x4) @ T_w_ip1_4x4
        
        # Extract 3x4 transformation
        relative_poses.append(T_i_ip1[:3, :])
    
    return np.array(relative_poses)

def process_sequence(pose_file):
    print(f"\nProcessing sequence {pose_file.stem}...")
    
    # Calculate relative poses
    relative_poses = process_pose_file(pose_file)
    
    # Create output filename
    output_file = pose_file.parent / f"{pose_file.stem}_rel.txt"
    
    # Save relative poses
    np.savetxt(output_file, relative_poses.reshape(-1, 12), fmt='%.6f')
    return len(relative_poses)

def main():
    # Set paths
    pose_dir = Path("/home/neo/Epipolar_evaluation/dataset/poses")
    
    # Process all txt files in the poses directory
    pose_files = list(pose_dir.glob("*.txt"))
    # Filter out files that already end with _rel
    pose_files = [f for f in pose_files if not f.stem.endswith("_rel")]
    
    # Use all available CPU cores
    n_cores = cpu_count()
    print(f"Using {n_cores} CPU cores")
    
    with Pool(processes=n_cores) as pool:
        results = list(tqdm(
            pool.imap_unordered(process_sequence, pose_files),
            total=len(pose_files),
            desc="Processing sequences"
        ))
    
    print(f"\nCompleted processing {len(pose_files)} sequences")
    print(f"Total relative poses generated: {sum(results)}")

if __name__ == "__main__":
    main()