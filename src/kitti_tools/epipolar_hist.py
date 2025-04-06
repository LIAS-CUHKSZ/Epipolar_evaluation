import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import random
from tqdm import tqdm
import cv2

def calculate_epipolar_distances(pts1, pts2, F):
    """Calculate distances between points and their corresponding epipolar lines."""
    pts2_h = np.hstack([pts2, np.ones((pts2.shape[0], 1))])  # Homogeneous coordinates for pts2
    lines = np.dot(F, np.hstack([pts1, np.ones((pts1.shape[0], 1))]).T).T  # Epipolar lines in img2
    
    # Calculate distances
    distances = (lines[:, 0] * pts2_h[:, 0] + lines[:, 1] * pts2_h[:, 1] + lines[:, 2]) / \
                np.sqrt(lines[:, 0]**2 + lines[:, 1]**2)
    return distances

def save_histogram(distances, output_path, title="Epipolar Distance Histogram", bins=50):
    """Save a histogram of epipolar distances."""
    plt.figure()
    plt.hist(distances, bins=bins, range=(-4, 4), color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def load_eth3d_correspondences(images_txt_path):
    """Load ETH3D correspondence points from images.txt."""
    images = {}
    with open(images_txt_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            parts = line.split()
            if len(parts) == 10:  # Image info line
                image_id = int(parts[0])
                file_path = parts[9]
                images[image_id] = {'file_path': file_path, 'observations': {}}
            elif len(parts) > 0:  # Observations line
                observations = images[image_id]['observations']
                for i in range(0, len(parts), 3):
                    x, y, pt3d_id = map(float, parts[i:i+3])
                    if int(pt3d_id) != -1:
                        observations[int(pt3d_id)] = (x, y)
    return images

def process_eth3d_sequence(sequence_path):
    images_txt_path = sequence_path / 'dslr_calibration_undistorted/images.txt'
    eval_dir = sequence_path / 'epipolar_eval'
    eval_dir.mkdir(exist_ok=True)

    if not images_txt_path.exists():
        print(f"Skipping sequence {sequence_path.name}: images.txt not found.")
        return

    images = load_eth3d_correspondences(images_txt_path)
    image_ids = sorted(images.keys())
    all_distances = []

    # Process all pairs for accumulated histogram
    for i in tqdm(range(len(image_ids) - 1), desc=f"Processing {sequence_path.name}"):
        img1_id = image_ids[i]
        img2_id = image_ids[i + 1]

        pts1 = []
        pts2 = []
        for pt3d_id, (x1, y1) in images[img1_id]['observations'].items():
            if pt3d_id in images[img2_id]['observations']:
                x2, y2 = images[img2_id]['observations'][pt3d_id]
                pts1.append((x1, y1))
                pts2.append((x2, y2))

        pts1 = np.float32(pts1)
        pts2 = np.float32(pts2)

        if len(pts1) >= 8:
            F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, 4, 0.999)
            pts1 = pts1[mask.ravel() == 1]
            pts2 = pts2[mask.ravel() == 1]

            distances = calculate_epipolar_distances(pts1, pts2, F)
            all_distances.extend(distances)

    # Save accumulated histogram
    save_histogram(all_distances, eval_dir / 'overall_histogram.png', title="Overall Epipolar Distance Histogram", bins=50)

    # Randomly select 5 pairs to draw individual histograms
    random_pairs = random.sample(list(zip(image_ids[:-1], image_ids[1:])), min(5, len(image_ids) - 1))
    for img1_id, img2_id in random_pairs:
        pts1 = []
        pts2 = []
        for pt3d_id, (x1, y1) in images[img1_id]['observations'].items():
            if pt3d_id in images[img2_id]['observations']:
                x2, y2 = images[img2_id]['observations'][pt3d_id]
                pts1.append((x1, y1))
                pts2.append((x2, y2))

        pts1 = np.float32(pts1)
        pts2 = np.float32(pts2)

        if len(pts1) >= 8:
            F, mask = cv2.findFundamentalMat(pts1, pts2, cv2.RANSAC, 4.0, 0.999)
            pts1 = pts1[mask.ravel() == 1]
            pts2 = pts2[mask.ravel() == 1]

            distances = calculate_epipolar_distances(pts1, pts2, F)
            save_histogram(distances, eval_dir / f"{img1_id}_to_{img2_id}_histogram.png", 
                           title=f"Epipolar Distance Histogram ({img1_id} to {img2_id})", bins=50)

def main():
    dataset_dir = Path("/home/neo/Epipolar_evaluation/dataset")
    sequence_dirs = list(dataset_dir.glob("*"))

    for seq_dir in sequence_dirs:
        if seq_dir.is_dir():
            print(f"Processing sequence: {seq_dir.name}")
            process_eth3d_sequence(seq_dir)

if __name__ == "__main__":
    main()
