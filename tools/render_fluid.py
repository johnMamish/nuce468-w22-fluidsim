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
import mmap
import matplotlib
import argparse
import multiprocessing

from numpy.lib.stride_tricks import as_strided

from functools import partial


import tempfile


def to_color(frame):
    rMin = 0.0
    rMax = 0.5

    # Handle NaN
    for x, row in enumerate(frame):
        for y, el in enumerate(row):
            if math.isnan(el):
                frame[x][y] = 0.0

    # frame[frame == math.inf] = rMax
    # frame[frame == -1*math.inf] = rMin
    # frame[math.isnan(frame)] = rMin

    # Restrict to a set range
    # frame[frame < rMin] = rMin
    # frame[frame > rMax] = rMax

    # frame /= rMax

    r = frame.max() - frame.min()
    if (r < 0.1):
        r = 0.1

    frame = (frame - frame.min()) / r # guarantees all values [0. 1)
    # if r != 0:
    #     frame /= r# guarantees all values [0. 1)
    #     print(r)
    #frame = np.where(np.isnan(frame), 0, frame)
    # frame = (frame / -3) + 1 # map [0, 1) to (1, 2/3]
    color_size = list(frame.shape)
    color_size.append(3)
    color_size = tuple(color_size) # adds a third dimension to the frames (for rgb)
    hsv_vals = np.ones(color_size)
    hsv_vals[:,:,0] = frame
    rgb_vals = matplotlib.colors.hsv_to_rgb(hsv_vals)
    rgb_vals = (rgb_vals*255).astype(np.uint8)
    return rgb_vals

def boolean_to_bw(frame):
    color_size = list(frame.shape)
    color_size.append(3)
    color_size = tuple(color_size) # adds a third dimension to the frames (for rgb)
    rgb_vals = np.zeros(color_size)
    inv = np.invert(frame)
    rgb_vals[:,:,0] = inv.astype(np.uint8)
    rgb_vals[:,:,1] = inv.astype(np.uint8)
    rgb_vals[:,:,2] = inv.astype(np.uint8)
    return rgb_vals

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='render binary fluidsim file to video')

    ap.add_argument('input_file', type=str)
    ap.add_argument('output_file', type=str)

    args = ap.parse_args()

    # parse "header" from file
    # memory-map the file, size 0 means whole file
    sim_file = open(args.input_file, 'r+b')
    sim_mmap = mmap.mmap(sim_file.fileno(), 0, access=mmap.ACCESS_READ)
    MAGIC_STR = b'FSIM'
    headers_end = sim_mmap.find(MAGIC_STR)
    if (headers_end == -1):
        print("malformed file, string b\'FSIM\' not found to terminate headers")
        sys.exit(-1)

    # Very sloppy parsing
    header_strings = [h.strip() for h in sim_mmap[0:headers_end].split(b';')]
    headers = []
    for hb in header_strings:
        if b':' not in hb:
            continue

        hs = hb.decode('utf-8')
        name = hs.split(':')[0]
        datatype = hs.split(':')[1].strip()[0]

        # super hacky; not trying to write a compiler here
        dimstr = hs.split(':')[1].strip()[1:]
        dimstart = dimstr.find('[')
        dimend = dimstr.find(']')
        dimstr = dimstr[dimstart:dimend]
        dims = [int(d) for d in dimstr[1:].split(',')]

        headers.append((name, datatype, dims))

    print("Found headers: ")
    for h in headers:
        print(f"{h[0]}: {h[1]} {h[2]}")
    print()

    WIDTH = headers[0][2][0]
    HEIGHT = headers[0][2][1]

    # fourcc = cv2.VideoWriter_fourcc(*"XVID")
    fourcc = cv2.VideoWriter_fourcc('p', 'n', 'g', ' ')
    video = cv2.VideoWriter(args.output_file, fourcc, 40.0, (WIDTH, HEIGHT))

    # skip headers
    sim_file.read(headers_end + len(MAGIC_STR))

    # size of each frame
    field_decode_strs = []
    field_sizes = []
    frame_size = 0
    for field in headers:
        field_decode_str = field[1]*(field[2][0] * field[2][1])
        field_decode_strs.append(field_decode_str)
        field_sizes.append(struct.calcsize(field_decode_str))

    print(f"each frame is calculated to be {sum(field_sizes)} bytes")

    frame_count = 0
    for chunk in iter(partial(sim_file.read, sum(field_sizes)), b''):
        print(f"writing frame {frame_count}\r", end='')

        # Read array of floats
        FRAME_LEN = WIDTH*HEIGHT*4
        offs = 0

        frame_data = {}
        data_offset = 0
        for i in range(len(headers)):
            field_name = headers[i][0]
            subchunk = chunk[data_offset:data_offset + field_sizes[i]]
            frame_data[field_name] = np.array(struct.unpack(field_decode_strs[i], subchunk)).reshape(HEIGHT, WIDTH)
            data_offset += field_sizes[i]

        # Turn it into a picture and write it to the frame
        frame = (to_color(frame_data['curl']) * boolean_to_bw(frame_data['barriers'])).astype(np.uint8)
        video.write(frame)

        frame_count += 1
    print()
    video.release()
