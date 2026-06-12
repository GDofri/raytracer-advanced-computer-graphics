#!/usr/bin/env python3

import sys
from PIL import Image
import math
import argparse

def load_images(image_paths):
    return [Image.open(path) for path in image_paths]

def create_pyramid(images):
    # Determine the number of images
    n = len(images)

    # Determine the number of rows needed for the pyramid
    rows = 0
    count = 0
    while count < n:
        rows += 1
        count += rows + 1

    # Adjust rows if we exceed the number of images
    if count > n:
        rows -= 1

    # Determine the width of the output image
    max_width = max(img.width for img in images)
    total_width = max_width * (rows+1)

    # Determine the height of the output image
    max_height = max(img.height for img in images)
    total_height = max_height * rows

    # Create the output image
    output_image = Image.new("RGB", (total_width, total_height), (255, 255, 255))

    # Paste images into the output image in a pyramid structure
    img_index = 0
    for row in range(2, rows + 2):
        y = (row - 2) * max_height
        for col in range(row):
            if img_index >= n:
                break
            x = (total_width // 2) - (row * max_width // 2) + (col * max_width)
            output_image.paste(images[img_index], (x, y))
            img_index += 1

    return output_image

def save_image(image, output_path):
    print("Saving pyramid to: " + output_path)
    image.save(output_path)

def main():
    parser = argparse.ArgumentParser(description="Create a pyramid of images.")
    parser.add_argument('image_paths', nargs='+', help="Paths to the input images.")
    parser.add_argument('--output', required=True, help="Path to the output image.")

    args = parser.parse_args()
    images = load_images(args.image_paths)
    output_image = create_pyramid(images)
    save_image(output_image, args.output)

if __name__ == "__main__":
    main()