import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import random
import multiprocessing
from functools import partial
import cv2

def load_matches(match_file):
    matches = []
    with open(match_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            img_indices = parts[0].split()
            pts1_str = parts[1].strip().split()
            pts2_str = parts[2].strip().split()
            
            pts1 = np.array([[float(pts1_str[i]), float(pts1_str[i+1])] 
                            for i in range(0, len(pts1_str), 2)])
            pts2 = np.array([[float(pts2_str[i]), float(pts2_str[i+1])] 
                            for i in range(0, len(pts2_str), 2)])
            
            matches.append({
                'img1': int(img_indices[0]),
                'img2': int(img_indices[1]),
                'pts1': pts1,
                'pts2': pts2
            })
    return matches

def load_relative_poses(pose_file):
    poses = np.loadtxt(pose_file)
    return poses.reshape(-1, 3, 4)

def compute_fundamental_matrix(pose):
    """Compute fundamental matrix from relative pose [R|t] and calibration"""
    # Camera intrinsics
    K = np.eye(3)
    K[0, 0] = K[1, 1] = 718.856  # fx, fy
    K[0, 2] = 607.193  # cx
    K[1, 2] = 185.216  # cy

    R = pose[:, :3]
    t = pose[:, 3]
    
    # Skew-symmetric matrix of translation
    t_cross = np.array([
        [0, -t[2], t[1]],
        [t[2], 0, -t[0]],
        [-t[1], t[0], 0]
    ])
    
    # Essential matrix = t_cross * R
    E = t_cross @ R
    
    # Fundamental matrix = K'^(-T) * E * K^(-1)
    # Since both cameras have the same K, K' = K
    K_inv = np.linalg.inv(K)
    F = K_inv.T @ E @ K_inv
    
    return F

def point_to_line_distance(point, line):
    """
    Compute signed distance from point to line
    line: coefficients [a, b, c] of ax + by + c = 0
    point: [x, y]
    """
    a, b, c = line
    x, y = point[0]
    dist = (a*x + b*y + c) / np.sqrt(a*a + b*b)
    return dist

def compute_epipolar_distances(matches, poses):
    all_distances = []
    
    for match, pose in zip(matches, poses):
        F = compute_fundamental_matrix(pose)
        
        distances = []
        pts1 = match['pts1'].reshape(-1, 1, 2)
        pts2 = match['pts2'].reshape(-1, 1, 2)
        
        # Compute epipolar lines using cv2
        lines = cv2.computeCorrespondEpilines(pts1, 1, F)
        lines = lines.reshape(-1, 3)
        
        # Calculate distances
        for pt, line in zip(pts2, lines):
            dist = point_to_line_distance(pt, line)
            distances.append(dist)
            
        all_distances.append({
            'img_pair': (match['img1'], match['img2']),
            'distances': np.array(distances)
        })
    
    return all_distances

def plot_histogram(distances, img_pair, output_dir):
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=50, range=(-5, 5))
    plt.title(f'Epipolar Distance Distribution for Image Pair {img_pair}')
    plt.xlabel('Distance (pixels)')
    plt.ylabel('Count')
    plt.grid(True)
    plt.savefig(output_dir / f'hist_{img_pair[0]}_{img_pair[1]}.png')
    plt.close()

def draw_epipolar_lines(img1, img2, pts1, pts2, lines, pair_index):
    """Draw epipolar lines and points on images"""
    if img1.ndim == 2:
        img1 = cv2.cvtColor(img1, cv2.COLOR_GRAY2BGR)
    if img2.ndim == 2:
        img2 = cv2.cvtColor(img2, cv2.COLOR_GRAY2BGR)
        
    h, w = img1.shape[:2]
    
    # Draw points and lines
    for pt1, pt2, line in zip(pts1, pts2, lines):
        x0, y0 = map(int, [0, -line[2]/line[1]])
        x1, y1 = map(int, [w, -(line[2] + line[0]*w)/line[1]])
        
        # Draw epipolar line
        color = tuple(np.random.randint(0, 255, 3).tolist())
        cv2.line(img2, (x0, y0), (x1, y1), color, 1)
        
        # Draw points
        cv2.circle(img1, tuple(map(int, pt1.ravel())), 5, color, -1)
        cv2.circle(img2, tuple(map(int, pt2.ravel())), 5, color, -1)
    
    return img1, img2

def process_sequence(seq_dir, pose_dir):
    print(f"\nProcessing sequence {seq_dir.name}")
    
    # Load matches
    match_file = seq_dir / 'eval' / 'matches.txt'
    if not match_file.exists():
        print(f"No matches file found for sequence {seq_dir.name}")
        return
        
    # Load poses
    pose_file = pose_dir / f"{seq_dir.name}_rel.txt"
    if not pose_file.exists():
        print(f"No pose file found for sequence {seq_dir.name}")
        return
        
    matches = load_matches(match_file)
    poses = load_relative_poses(pose_file)
    
    # Compute distances
    distances = compute_epipolar_distances(matches, poses)
    
    # Create output directory
    output_dir = seq_dir / 'eval' / 'epipolar'
    output_dir.mkdir(exist_ok=True)
    
    # Save all distances
    all_dists = np.concatenate([d['distances'] for d in distances])
    np.save(output_dir / 'all_distances.npy', all_dists)
    
    # Plot histogram and draw epipolar lines for randomly selected 5 pairs
    random_pairs = random.sample(distances, min(5, len(distances)))
    for pair_data in random_pairs:
        plot_histogram(pair_data['distances'], pair_data['img_pair'], output_dir)
        
        # Load images and draw epipolar lines
        img_idx1, img_idx2 = pair_data['img_pair']
        img1 = cv2.imread(str(seq_dir / 'image_0' / f'{img_idx1:06d}.png'))
        img2 = cv2.imread(str(seq_dir / 'image_0' / f'{img_idx2:06d}.png'))
        
        if img1 is not None and img2 is not None:
            match_idx = matches.index(next(m for m in matches if m['img1'] == img_idx1 and m['img2'] == img_idx2))
            pts1 = matches[match_idx]['pts1'].reshape(-1, 1, 2)
            pts2 = matches[match_idx]['pts2'].reshape(-1, 1, 2)
            F = compute_fundamental_matrix(poses[match_idx])
            lines = cv2.computeCorrespondEpilines(pts1, 1, F)
            lines = lines.reshape(-1, 3)
            
            # Draw and save visualization
            img1_draw, img2_draw = draw_epipolar_lines(img1, img2, pts1, pts2, lines, pair_data['img_pair'])
            cv2.imwrite(str(output_dir / f'epipolar_{img_idx1}_{img_idx2}_img1.png'), img1_draw)
            cv2.imwrite(str(output_dir / f'epipolar_{img_idx1}_{img_idx2}_img2.png'), img2_draw)
    
    # Plot overall histogram
    plt.figure(figsize=(10, 6))
    plt.hist(all_dists, bins=60, range=(-7, 7))
    plt.title(f'Overall Epipolar Distance Distribution for Sequence {seq_dir.name}')
    plt.xlabel('Distance (pixels)')
    plt.ylabel('Count')
    plt.grid(True)
    plt.savefig(output_dir / 'hist_overall.png')
    plt.close()
    
    print(f"Processed {len(matches)} image pairs")
    print(f"Mean distance: {np.mean(np.abs(all_dists)):.3f} pixels")
    print(f"Median distance: {np.median(np.abs(all_dists)):.3f} pixels")

def main():
    dataset_dir = Path("/home/neo/Epipolar_evaluation/dataset")
    pose_dir = Path("/home/neo/Epipolar_evaluation/dataset/poses")
    
    # Process each sequence
    sequence_dirs = sorted([d for d in dataset_dir.glob("[0-9][0-9]") if d.is_dir()])
    
    # Create a process pool and process sequences in parallel
    with multiprocessing.Pool() as pool:
        pool.map(partial(process_sequence, pose_dir=pose_dir), sequence_dirs)
    # without parallelism
    # for seq_dir in sequence_dirs:
    #     process_sequence(seq_dir, pose_dir)
    

if __name__ == "__main__":
    main()