import cv2
import numpy as np
from pathlib import Path
import os
from multiprocessing import Pool, cpu_count
from itertools import islice
from tqdm import tqdm
import matplotlib.pyplot as plt
import random
import argparse
from scipy.stats import norm

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
        # Normalize the line parameters
        line /= np.sqrt(line[0]**2 + line[1]**2)
        
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
    """Save a histogram of epipolar distances with a Gaussian fit."""
    plt.figure()
    # Plot histogram and get bin data
    counts, bins, _ = plt.hist(distances, bins=50, range=(-2, 2), color='blue', alpha=0.7, density=True)
    
    # Discretize data for Gaussian fitting
    bin_centers = (bins[:-1] + bins[1:]) / 2
    weights = counts * np.diff(bins)  # Use bin heights as weights
    mu, std = norm.fit(bin_centers, floc=np.average(bin_centers, weights=weights), fscale=np.sqrt(np.average((bin_centers - np.average(bin_centers, weights=weights))**2, weights=weights)))
    
    # Generate Gaussian curve
    x = np.linspace(bins[0], bins[-1], 100)
    p = norm.pdf(x, mu, std)
    
    # Plot Gaussian curve
    plt.plot(x, p, 'r--', linewidth=2, label=f'Gaussian Fit\n$\mu={mu:.2f}, \sigma={std:.2f}$')
    
    # Add labels and title
    plt.title("Epipolar Distance Histogram")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)
    
    # Save the figure
    plt.savefig(output_path)
    plt.close()

def extract_and_match_features(image_pair, eval_dir, detector, threshold, confidence):
    img1_path, img2_path = image_pair
    # Read images
    img1 = cv2.imread(str(img1_path), cv2.IMREAD_GRAYSCALE)
    img2 = cv2.imread(str(img2_path), cv2.IMREAD_GRAYSCALE)
    
    if detector == "lkof":
        # Optical flow-based feature matching
        feature_params = dict(maxCorners=2500, qualityLevel=0.01, minDistance=7, blockSize=7)
        pts1 = cv2.goodFeaturesToTrack(img1, mask=None, **feature_params)
        if pts1 is None:
            return img1_path.stem, img2_path.stem, None, None, None, None
        pts1 = np.float32(pts1).reshape(-1, 2)
        
        lk_params = dict(winSize=(15, 15), maxLevel=2, criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 0.03))
        pts2, st, err = cv2.calcOpticalFlowPyrLK(img1, img2, pts1, None, **lk_params)
        pts1_back, st_back, err_back = cv2.calcOpticalFlowPyrLK(img2, img1, pts2, None, **lk_params)
        
        d = np.linalg.norm(pts1 - pts1_back, axis=1)
        valid = d < threshold
        pts1 = pts1[valid]
        pts2 = pts2[valid]
        
        # RANSAC to refine optical flow matches
        if len(pts1) >= 8:
            F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, threshold, confidence)
            pts1 = pts1[mask.ravel() == 1]
            pts2 = pts2[mask.ravel() == 1]
            distances = calculate_epipolar_distances(pts1, pts2, F)
            return img1_path.stem, img2_path.stem, pts1, pts2, distances, F
    else:
        # Feature detector-based matching
        if detector == "surf":
            feature_detector = cv2.xfeatures2d.SURF_create(hessianThreshold=400)
        elif detector == "orb":
            feature_detector = cv2.ORB_create(nfeatures=2500)  # Set max features for ORB
        elif detector == "sift":
            feature_detector = cv2.SIFT_create(nfeatures=2500)  # Set max features for SIFT
        else:
            raise ValueError(f"Unsupported feature detector: {detector}")
        
        kp1, des1 = feature_detector.detectAndCompute(img1, None)
        kp2, des2 = feature_detector.detectAndCompute(img2, None)
        
        if detector in ["surf", "sift"]:
            FLANN_INDEX_KDTREE = 1
            index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
            search_params = dict(checks=50)
            flann = cv2.FlannBasedMatcher(index_params, search_params)
            matches = flann.knnMatch(des1, des2, k=2)
            good_matches = [m for m, n in matches if m.distance < 0.7 * n.distance]
        else:
            bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
            good_matches = bf.match(des1, des2)
            good_matches = sorted(good_matches, key=lambda x: x.distance)
        
        pts1 = np.float32([kp1[m.queryIdx].pt for m in good_matches])
        pts2 = np.float32([kp2[m.trainIdx].pt for m in good_matches])
    
    # RANSAC
    if len(pts1) >= 8:
        F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, threshold, confidence)
        pts1 = pts1[mask.ravel() == 1]
        pts2 = pts2[mask.ravel() == 1]
        distances = calculate_epipolar_distances(pts1, pts2, F)
        return img1_path.stem, img2_path.stem, pts1, pts2, distances, F
    
    return img1_path.stem, img2_path.stem, None, None, None, None

def process_image_pair(args):
    """Wrapper function to process an image pair with eval_dir."""
    image_pair, eval_dir, detector, threshold, confidence = args
    return extract_and_match_features(image_pair, eval_dir, detector, threshold, confidence)

def process_sequence(sequence_path, detector, threshold, confidence):
    img_dir = sequence_path / 'image_0'
    eval_dir = sequence_path / f'{threshold}_{int(confidence * 100)}_{detector}'
    os.makedirs(eval_dir, exist_ok=True)
    
    output_file = eval_dir / 'matches.txt'
    stats_file = eval_dir / 'stats.txt'
    histogram_file = eval_dir / 'distance_histogram.png'
    
    image_files = sorted(img_dir.glob('*.png'))
    image_pairs = list(zip(image_files[:-1], image_files[1:]))
    
    n_cores = cpu_count()
    print(f"Using {n_cores} CPU cores")
    
    results = []
    num_points = []
    all_distances = []
    
    with Pool(processes=n_cores) as pool:
        args = [(pair, eval_dir, detector, threshold, confidence) for pair in image_pairs]
        for result in tqdm(pool.imap_unordered(process_image_pair, args), 
                          total=len(image_pairs), 
                          desc=f"Processing sequence {sequence_path.name}"):
            img1_idx, img2_idx, pts1, pts2, distances, F = result
            if pts1 is not None and len(pts1) > 0:
                results.append((img1_idx, img2_idx, pts1, pts2, F))
                num_points.append(len(pts1))
                all_distances.extend(distances)
    
    if all_distances:
        save_histogram(all_distances, histogram_file)
    
    random_pairs = random.sample(results, min(5, len(results)))
    for img1_idx, img2_idx, pts1, pts2, F in random_pairs:
        output_path = eval_dir / f"{img1_idx}_to_{img2_idx}_epilines.png"
        residual_hist_path = eval_dir / f"{img1_idx}_to_{img2_idx}_residual_hist.png"
        draw_epipolar_lines(cv2.imread(str(img_dir / f"{img1_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            cv2.imread(str(img_dir / f"{img2_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            pts1, pts2, F, str(output_path))
        save_histogram(calculate_epipolar_distances(pts1, pts2, F), residual_hist_path)
    
    if num_points:
        avg_points = sum(num_points) / len(num_points)
        min_points = min(num_points)
        max_points = max(num_points)
        
        with open(stats_file, 'w') as f:
            f.write(f"Sequence Statistics:\n")
            f.write(f"Average points per pair: {avg_points:.2f}\n")
            f.write(f"Minimum points in a pair: {min_points}\n")
            f.write(f"Maximum points in a pair: {max_points}\n")
            f.write(f"Total image pairs processed: {len(num_points)}\n")
    
    results.sort(key=lambda x: int(x[0]))
    
    with open(output_file, 'w') as f:
        for img1_idx, img2_idx, pts1, pts2, _ in results:
            f.write(f"{img1_idx} {img2_idx}")
            f.write(',')
            for pt in pts1:
                f.write(f" {pt[0]:.2f} {pt[1]:.2f}")
            f.write(',')
            for pt in pts2:
                f.write(f" {pt[0]:.2f} {pt[1]:.2f}")
            f.write("\n")

def main():
    print("\033[92mif u want to use SURF, please install python3.6 and opencv & opencv-contrib-python==3.4.1.15\033[0m")
    parser = argparse.ArgumentParser(description="Epipolar Geometry Evaluation")
    parser.add_argument("--dataset_dir", type=str, default="../../dataset", help="Path to the dataset directory")
    parser.add_argument("--detector", type=str, choices=["surf", "orb", "sift","lkof"], default="sift", help="Feature detector to use")
    parser.add_argument("--threshold", type=float, default=1.0, help="RANSAC or bi-direction optiflow error threshold")
    parser.add_argument("--confidence", type=float, default=0.999, help="RANSAC confidence level")
    args = parser.parse_args()
    
    print(f"detector: {args.detector}")
    print(f"threshold: {args.threshold}")
    print(f"confidence: {args.confidence}")
    print(f"dataset directory: {args.dataset_dir}")
    
    dataset_dir = Path(args.dataset_dir)
    sequence_dirs = list(dataset_dir.glob("[0-9][0-9]"))
    for i, seq_dir in enumerate(sequence_dirs):
        if seq_dir.is_dir():
            print(f"\nProcessing sequence {i+1}/{len(sequence_dirs)}: {seq_dir.name}...")
            process_sequence(seq_dir, args.detector, args.threshold, args.confidence)

if __name__ == "__main__":
    main()