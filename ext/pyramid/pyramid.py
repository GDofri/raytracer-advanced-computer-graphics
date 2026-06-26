#!/usr/bin/env python3

import re
import sys
from PIL import Image, ImageDraw, ImageFont
import math
import argparse

def load_images(image_paths):
    return [Image.open(path) for path in image_paths]

def parse_st(path):
    """Extract (s, t) from a filename like ..._sN_tM.png, or return None."""
    m = re.search(r'_s(\d+)_t(\d+)', path)
    return (int(m.group(1)), int(m.group(2))) if m else None

def _load_font(size):
    for name in ["DejaVuSans.ttf", "DejaVuSansMono.ttf", "LiberationSans-Regular.ttf",
                 "FreeSans.ttf", "Arial.ttf", "Helvetica.ttf"]:
        try:
            return ImageFont.truetype(name, size)
        except (IOError, OSError):
            pass
    return ImageFont.load_default()

def label_image(img, text):
    """Return a copy of img with text stamped in the bottom-right corner."""
    img = img.copy()
    draw = ImageDraw.Draw(img)
    font = _load_font(max(12, img.width // 10))
    pad = 4
    try:
        bbox = draw.textbbox((0, 0), text, font=font)
    except ValueError:
        bbox = (0, 0, len(text) * 6, 11)
    # Anchor the text so its actual bottom-right corner sits pad px from the image edge
    tx = img.width  - bbox[2] - 2 * pad
    ty = img.height - bbox[3] - 2 * pad
    rect = [tx + bbox[0] - pad, ty + bbox[1] - pad,
            tx + bbox[2] + pad, ty + bbox[3] + pad]
    draw.rectangle(rect, fill=(0, 0, 0))
    draw.text((tx, ty), text, fill=(255, 255, 255), font=font)
    return img

def create_pyramid(images, labels=None):
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
            img = images[img_index]
            if labels and img_index < len(labels) and labels[img_index]:
                img = label_image(img, labels[img_index])
            output_image.paste(img, (x, y))
            img_index += 1

    return output_image

def save_image(image, output_path):
    print("Saving pyramid to: " + output_path)
    image.save(output_path)

def main():
    parser = argparse.ArgumentParser(description="Create a pyramid of images.")
    parser.add_argument('image_paths', nargs='+', help="Paths to the input images.")
    parser.add_argument('--output', required=True, help="Path to the output image.")
    parser.add_argument('--label', action='store_true',
                        help="Stamp each image with its s/t strategy indices (parsed from filename).")

    args = parser.parse_args()
    images = load_images(args.image_paths)
    labels = None
    if args.label:
        labels = []
        for p in args.image_paths:
            st = parse_st(p)
            labels.append(f"s={st[0]} t={st[1]}" if st else None)
    output_image = create_pyramid(images, labels=labels)
    save_image(output_image, args.output)

if __name__ == "__main__":
    main()