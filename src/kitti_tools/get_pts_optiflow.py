import cv2
import numpy as np
from pathlib import Path
import os
from tqdm import tqdm
import random
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

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
        cv2.circle(img2_color, tuple(map(int, pt1)), 3, color, -1)
    
    # Save the result
    cv2.imwrite(output_path, img2_color)

def extract_and_match_features_optiflow(image_pair, eval_dir):
    img1_path, img2_path = image_pair
    # Read images
    img1 = cv2.imread(str(img1_path), cv2.IMREAD_GRAYSCALE)
    img2 = cv2.imread(str(img2_path), cv2.IMREAD_GRAYSCALE)

    # Detect good features to track in the first image
    feature_params = dict(maxCorners=1000, qualityLevel=0.01, minDistance=7, blockSize=7)
    pts1 = cv2.goodFeaturesToTrack(img1, mask=None, **feature_params)

    if pts1 is None:
        return img1_path.stem, img2_path.stem, None, None, None, None

    pts1 = np.float32(pts1).reshape(-1, 2)

    # Calculate optical flow to find corresponding points in the second image
    lk_params = dict(winSize=(15, 15), maxLevel=2, criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 0.03))
    pts2, st, err = cv2.calcOpticalFlowPyrLK(img1, img2, pts1, None, **lk_params)

    # Perform inverse optical flow to filter outliers
    pts1_back, st_back, err_back = cv2.calcOpticalFlowPyrLK(img2, img1, pts2, None, **lk_params)
    d = np.linalg.norm(pts1 - pts1_back, axis=1)
    valid = d < 1.0  # Threshold for bidirectional consistency
    pts1 = pts1[valid]
    pts2 = pts2[valid]

    # RANSAC to further refine matches
    if len(pts1) >= 8:
        F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, 2.0, 0.999)
        pts1 = pts1[mask.ravel() == 1]
        pts2 = pts2[mask.ravel() == 1]

        # Calculate epipolar distances
        distances = calculate_epipolar_distances(pts1, pts2, F)
        
        # discard points with distance larger than 2.5 pixels
        valid = np.abs(distances) < 2.5
        pts1 = pts1[valid]
        pts2 = pts2[valid]
        distances = distances[valid]

        return img1_path.stem, img2_path.stem, pts1, pts2, distances, F

    return img1_path.stem, img2_path.stem, None, None, None, None

def process_image_pair(args):
    """Wrapper function to process an image pair with eval_dir."""
    image_pair, eval_dir = args
    return extract_and_match_features_optiflow(image_pair, eval_dir)

def process_sequence(sequence_path):
    img_dir = sequence_path / 'image_0'
    eval_dir = sequence_path / 'eval_optiflow_sac2_99'
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
        args = [(pair, eval_dir) for pair in image_pairs]
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
        draw_epipolar_lines(cv2.imread(str(img_dir / f"{img1_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            cv2.imread(str(img_dir / f"{img2_idx}.png"), cv2.IMREAD_GRAYSCALE),
                            pts1, pts2, F, str(output_path))

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
    dataset_dir = Path("/home/neo/Epipolar_evaluation/dataset")
    sequence_dirs = list(dataset_dir.glob("[0-9][0-9]"))
    for i, seq_dir in enumerate(sequence_dirs):
        if seq_dir.is_dir():
            print(f"\nProcessing sequence {i+1}/{len(sequence_dirs)}: {seq_dir.name}...")
            process_sequence(seq_dir)

if __name__ == "__main__":
    main()