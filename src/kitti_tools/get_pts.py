import cv2
import numpy as np
from pathlib import Path
import os
from multiprocessing import Pool, cpu_count
from itertools import islice

def extract_and_match_features(image_pair):
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
        F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, 2.0)
        
        # Select only inlier points
        pts1 = pts1[mask.ravel() == 1]
        pts2 = pts2[mask.ravel() == 1]
        
        return img1_path.stem, img2_path.stem, pts1, pts2
    
    return img1_path.stem, img2_path.stem, None, None

def process_sequence(sequence_path):
    img_dir = sequence_path / 'image_0'
    output_file = sequence_path / 'matches.txt'
    
    # Get sorted list of image files
    image_files = sorted(img_dir.glob('*.png'))
    
    # Create pairs of consecutive images
    image_pairs = list(zip(image_files[:-1], image_files[1:]))
    
    # Process pairs in parallel
    n_cores = cpu_count()
    print(f"Using {n_cores} CPU cores")
    
    results = []
    with Pool(processes=n_cores) as pool:
        for result in pool.imap_unordered(extract_and_match_features, image_pairs):
            img1_idx, img2_idx, pts1, pts2 = result
            if pts1 is not None and len(pts1) > 0:
                results.append((img1_idx, img2_idx, pts1, pts2))
                print(f"Processed {img1_idx}-{img2_idx}: Found {len(pts1)} matches")
    
    # Sort results by image index
    results.sort(key=lambda x: int(x[0]))
    
    # Write results to file
    with open(output_file, 'w') as f:
        for img1_idx, img2_idx, pts1, pts2 in results:
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
    for seq_dir in dataset_dir.glob("[0-9][0-9]"):
        if seq_dir.is_dir():
            print(f"\nProcessing sequence {seq_dir.name}...")
            process_sequence(seq_dir)

if __name__ == "__main__":
    main()