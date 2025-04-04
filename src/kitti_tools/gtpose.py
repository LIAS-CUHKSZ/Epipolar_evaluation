import numpy as np
import os
from pathlib import Path

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
    for i in range(len(transforms)-1):
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

def main():
    # Set paths
    pose_dir = Path("/home/neo/Epipolar_evaluation/dataset/poses")
    
    # Process all txt files in the poses directory
    for pose_file in pose_dir.glob("*.txt"):
        if pose_file.stem.endswith("_rel"):
            continue
            
        print(f"Processing {pose_file.name}...")
        
        # Calculate relative poses
        relative_poses = process_pose_file(pose_file)
        
        # Create output filename
        output_file = pose_file.parent / f"{pose_file.stem}_rel.txt"
        
        # Save relative poses
        np.savetxt(output_file, relative_poses.reshape(-1, 12), fmt='%.6f')
        print(f"Saved relative poses to {output_file.name}")

if __name__ == "__main__":
    main()