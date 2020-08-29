import numpy as np
import cv2


# class Panorama():
#     def __init__(self, images, num_of_interest_points=1000, nn_threshold=0.7):
#         self.num_of_interest_points = num_of_interest_points
#         self.nn_threshold = nn_threshold
#         self.images = images

def cylindrical_projection(image, focal_length):
    # Project image to cylindrical coordinate
    rows, cols, _ = image.shape
    image_projected = np.zeros(image.shape, dtype=np.uint8)

    for r in range(rows):
        for c in range(cols):
            # Use center as (0, 0)
            y = round(r - rows/2)
            x = round(c - cols/2)
            x_proj = focal_length * np.arctan(x/focal_length)
            y_proj = (focal_length * y) / np.sqrt(np.square(x) + np.square(focal_length))

            # Use top left as (0, 0)
            x_proj = round(x_proj + cols/2)
            y_proj = round(y_proj + rows/2)

            image_projected[y_proj, x_proj, :] = image[r, c, :]

    return image_projected


def extract_features(image):
    pass


def multiscale_hearris(image):
    pass


def adaptive_non_maximal_suppression():
    pass


def msop_descriptor():
    pass


def feature_matching():
    pass


def disp_matched_features():
    pass


def image_matching(method='RANSAC'):
    pass


def blend():
    pass


def stitch(images, num_of_interest_points=1000, nn_threshold=0.7):
    num_of_pics = len(images)

    for i in range(num_of_pics):
        images[i] = cylindrical_projection(images[i])


def read_images(src_dir, num_of_pics):
    images = []
    for i in range(1, num_of_pics+1):
        img = cv2.imread(f"{src_dir}/{i}.jpg")
        images.append(img)
    return images


if __name__ == "__main__":
    images = read_images(src_dir="./imageset/", num_of_pics=18)
