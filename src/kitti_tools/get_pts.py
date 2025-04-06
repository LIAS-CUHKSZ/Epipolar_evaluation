import cv2
import numpy as np
from pathlib import Path
import os
from multiprocessing import Pool, cpu_count
from itertools import islice
from tqdm import tqdm
import matplotlib.pyplot as plt
import random

def draw_epipolar_lines(img1, img2, pts1, pts2, F, output_path):
    """Draw epipolar lines on the second image corresponding to points in the first image."""
    # Convert points to homogeneous coordinates
    pts1_h = np.hstack([pts1, np.ones((pts1.shape[0], 1))])
    
    # Compute epipolar lines in the second image
    lines = np.dot(F, pts1_h.T).T
    
    # Normalize the lines for better visualization
    lines /= np.sqrt(lines[:, 0]**2 + lines[:, 1]**2).reshape(-1, 1)
    
    # Draw the lines on the second image
    img2_color = cv2.cvtColor(img2, cv2.COLOR_GRAY2BGR)
    colors = [tuple(np.random.randint(0, 255, 3).tolist()) for _ in range(len(pts1))]  # Random colors
    
    for line, pt1, color in zip(lines, pts1, colors):
        x0, y0 = map(int, [0, -line[2] / line[1]])
        x1, y1 = map(int, [img2.shape[1], -(line[2] + line[0] * img2.shape[1]) / line[1]])
        cv2.line(img2_color, (x0, y0), (x1, y1), color, 1)
        cv2.circle(img2_color, tuple(map(int, pt1)), 5, color, -1)
    
    # Save the result
    cv2.imwrite(output_path, img2_color)

def calculate_epipolar_distances(pts1, pts2, F):
    """Calculate distances between points and their corresponding epipolar lines."""
    pts2_h = np.hstack([pts2, np.ones((pts2.shape[0], 1))])  # Homogeneous coordinates for pts2
    lines = np.dot(F, np.hstack([pts1, np.ones((pts1.shape[0], 1))]).T).T  # Epipolar lines in img2
    
    # Calculate distances
    distances = (lines[:, 0] * pts2_h[:, 0] + lines[:, 1] * pts2_h[:, 1] + lines[:, 2]) / \
                np.sqrt(lines[:, 0]**2 + lines[:, 1]**2)
    return distances

def save_histogram(distances, output_path):
    """Save a histogram of epipolar distances."""
    plt.figure()
    plt.hist(distances, bins=50, range=(-2, 2), color='blue', alpha=0.7)
    plt.title("Epipolar Distance Histogram")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def extract_and_match_features(image_pair, eval_dir):
    img1_path, img2_path = image_pair
    # Read images
    img1 = cv2.imread(str(img1_path), cv2.IMREAD_GRAYSCALE)
    img2 = cv2.imread(str(img2_path), cv2.IMREAD_GRAYSCALE)
    
    # Initialize SIFT detector
    sift = cv2.SIFT_create()
    
    # Find keypoints and descriptors
    kp1, des1 = sift.detectAndCompute(img1, None)
    kp2, des2 = sift.detectAndCompute(img2, None)
    
    # FLANN matcher
    FLANN_INDEX_KDTREE = 1
    index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
    search_params = dict(checks=50)
    flann = cv2.FlannBasedMatcher(index_params, search_params)
    
    matches = flann.knnMatch(des1, des2, k=2)
    
    # Lowe's ratio test
    good_matches = []
    pts1 = []
    pts2 = []
    
    for m, n in matches:
        if m.distance < 0.7 * n.distance:
            good_matches.append(m)
            pts1.append(kp1[m.queryIdx].pt)
            pts2.append(kp2[m.trainIdx].pt)
    
    pts1 = np.float32(pts1)
    pts2 = np.float32(pts2)
    
    # RANSAC
    if len(pts1) >= 8:
        F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, 2.0, 0.999)
        
        # Select only inlier points
        pts1 = pts1[mask.ravel() == 1]
        pts2 = pts2[mask.ravel() == 1]
        
        # Calculate epipolar distances
        distances = calculate_epipolar_distances(pts1, pts2, F)
        
        return img1_path.stem, img2_path.stem, pts1, pts2, distances, F
    
    return img1_path.stem, img2_path.stem, None, None, None, None

def process_image_pair(args):
    """Wrapper function to process an image pair with eval_dir."""
    image_pair, eval_dir = args
    return extract_and_match_features(image_pair, eval_dir)

def process_sequence(sequence_path):
    img_dir = sequence_path / 'image_0'
    # create a eval dir under sequence path
    eval_dir = sequence_path / 'evaltest'
    os.makedirs(eval_dir, exist_ok=True)
    
    output_file = eval_dir / 'matches.txt'
    stats_file = eval_dir / 'stats.txt'
    histogram_file = eval_dir / 'distance_histogram.png'
    
    # Get sorted list of image files
    image_files = sorted(img_dir.glob('*.png'))
    
    # Create pairs of consecutive images
    image_pairs = list(zip(image_files[:-1], image_files[1:]))
    
    # Process pairs in parallel
    n_cores = cpu_count()
    print(f"Using {n_cores} CPU cores")
    
    results = []
    num_points = []  # Track number of points for statistics
    all_distances = []  # Collect all distances for histogram
    
    with Pool(processes=16) as pool:
        # Add tqdm progress bar
        args = [(pair, eval_dir) for pair in image_pairs]
        for result in tqdm(pool.imap_unordered(process_image_pair, args), 
                          total=len(image_pairs), 
                          desc=f"Processing sequence {sequence_path.name}"):
            img1_idx, img2_idx, pts1, pts2, distances, F = result
            if pts1 is not None and len(pts1) > 0:
                results.append((img1_idx, img2_idx, pts1, pts2, F))
                num_points.append(len(pts1))
                all_distances.extend(distances)
    
    # Save histogram of distances
    if all_distances:
        save_histogram(all_distances, histogram_file)
    
    # Randomly select 5 pairs to draw epipolar lines
    random_pairs = random.sample(results, min(5, len(results)))
    for img1_idx, img2_idx, pts1, pts2, F in random_pairs:
        output_path = eval_dir / f"{img1_idx}_to_{img2_idx}_epilines.png"
        draw_epipolar_lines(cv2.imread(str(img_dir / f"{img1_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            cv2.imread(str(img_dir / f"{img2_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            pts1, pts2, F, str(output_path))
    
    # Calculate statistics
    if num_points:
        avg_points = sum(num_points) / len(num_points)
        min_points = min(num_points)
        max_points = max(num_points)
        
        # Save statistics
        with open(stats_file, 'w') as f:
            f.write(f"Sequence Statistics:\n")
            f.write(f"Average points per pair: {avg_points:.2f}\n")
            f.write(f"Minimum points in a pair: {min_points}\n")
            f.write(f"Maximum points in a pair: {max_points}\n")
            f.write(f"Total image pairs processed: {len(num_points)}\n")
    
    # Sort results by image index
    results.sort(key=lambda x: int(x[0]))
    
    # Write results to file
    with open(output_file, 'w') as f:
        for img1_idx, img2_idx, pts1, pts2, _ in results:
            f.write(f"{img1_idx} {img2_idx}")
            f.write(',')
            
            # Write points from first image
            for pt in pts1:
                f.write(f" {pt[0]:.2f} {pt[1]:.2f}")
            f.write(',')
            # Write points from second image
            for pt in pts2:
                f.write(f" {pt[0]:.2f} {pt[1]:.2f}")
                
            f.write("\n")

def main():
    dataset_dir = Path("/home/neo/Epipolar_evaluation/dataset")
    
    # Process each sequence directory
    sequence_dirs = list(dataset_dir.glob("[0-9][0-9]"))
    for i, seq_dir in enumerate(sequence_dirs):
        if seq_dir.is_dir():
            print(f"\nProcessing sequence {i+1}/{len(sequence_dirs)}: {seq_dir.name}...")
            process_sequence(seq_dir)

if __name__ == "__main__":
    main()