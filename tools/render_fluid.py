#!/usr/bin/python3

# Converts a file

import colorsys
import cv2
import io
import itertools
import math
import numpy as np
import subprocess
from PIL import Image, ImageFile #Changed
import struct
import sys
import os
import time
import matplotlib
import argparse
import multiprocessing

from numpy.lib.stride_tricks import as_strided

from functools import partial


import tempfile


def to_color(frame):
    r = frame.max() - frame.min()
    if (r < 0.01):
        r = 0.01

    frame = (frame - frame.min()) / r # guarantees all values [0. 1)
    #frame = np.where(np.isnan(frame), 0, frame)
    frame = (frame / -3) + 1 # map [0, 1) to (1, 2/3]
    color_size = list(frame.shape)
    color_size.append(3)
    color_size = tuple(color_size) # adds a third dimension to the frames (for rgb)
    hsv_vals = np.ones(color_size)
    hsv_vals[:,:,0] = frame
    rgb_vals = matplotlib.colors.hsv_to_rgb(hsv_vals)
    rgb_vals = (rgb_vals*255).astype(np.uint8)
    return rgb_vals

if __name__ == "__main__":
    WIDTH = 200
    HEIGHT = 100

    ap = argparse.ArgumentParser(description='render binary fluidsim file to video')

    ap.add_argument('input_file', type=str)
    ap.add_argument('output_file', type=str)

    args = ap.parse_args()

    sim_file = open(args.input_file, 'rb')
    next(sim_file)
    next(sim_file)
    next(sim_file)

    # fourcc = cv2.VideoWriter_fourcc(*"XVID")
    fourcc = cv2.VideoWriter_fourcc('p', 'n', 'g', ' ')

    video = cv2.VideoWriter(args.output_file, fourcc, 20.0, (WIDTH, HEIGHT))

    i = 0
    for chunk in iter(partial(sim_file.read, 3 * WIDTH * HEIGHT * 4), b''):
        print(f"writing frame {i}\r", end='')

        # Read array of floats
        FRAME_LEN = WIDTH*HEIGHT*4
        offs = 0
        curl = np.array(struct.unpack('f'*(WIDTH*HEIGHT), chunk[FRAME_LEN*(offs):FRAME_LEN*(offs + 1)])).reshape((HEIGHT,WIDTH))
        print(curl)

        # Turn it into a picture and write it to the frame
        frame = to_color(curl)

        video.write(frame)

        i = i + 1

    video.release()
