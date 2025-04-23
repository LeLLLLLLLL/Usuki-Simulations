#!/usr/bin/env python
# coding: utf-8

# In[1]:


#USING SEABORN/WITHOUT PLOTLY
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.ticker as ticker
import cv2
import imageio
import re

frames_folder = "frames"
os.makedirs(frames_folder, exist_ok=True)
w = np.loadtxt("waves.txt")
t = np.loadtxt("tr_b.txt")
energy = t[:, 0]
transmission = t[:, 2]
total_frames = 199 #link variable (this is ne from prev code)
rows = 119 #link variable (this is rows from prev code)
nsl = 69 #link variable (this is nsl + 1 from prev code)

all_values = []
for frame in range(total_frames):
    w2 = np.zeros((rows, nsl))
    for i in range(nsl):
        for j in range(rows):
            w2[j, i] = w[nsl * rows * frame + i * rows + j]
    all_values.append(w2.flatten())

all_values = np.concatenate(all_values)
cmin = 0
cmax = np.percentile(all_values, 95)

filenames_combined = []

for frame in range(total_frames):
    w2 = np.zeros((rows, nsl))
    for i in range(nsl):
        for j in range(rows):
            w2[j, i] = w[nsl * rows * frame + i * rows + j]

    x = np.arange(nsl)
    y = np.arange(rows)
    X, Y = np.meshgrid(x, y)

    fig_3D = plt.figure(figsize=(5, 2.5))
    ax = fig_3D.add_subplot(121, projection='3d')

    surf = ax.plot_surface(X, Y, w2, cmap='jet', vmin=cmin, vmax=cmax, alpha=0.9, rstride=3, cstride=3)

    ax.set_title(f"Frame {frame + 1} - 3D View", fontsize=8)
    ax.set_xlabel("X Axis", fontsize=6, labelpad=8)
    ax.set_ylabel("Y Axis", fontsize=6, labelpad=8)
    ax.set_zlabel("Wave Amplitude", fontsize=6, labelpad=10)
    ax.tick_params(axis='both', labelsize=6)
    ax.w_zaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.1e}'))
    ax.view_init(elev=30, azim=150)
    local_cmax = np.max(w2)
    ax.set_zticks(np.linspace(0, local_cmax, 6))
    ax.tick_params(axis='z', pad=6, labelsize=6)
    ax.set_box_aspect([1, 1.5, 0.5])

    filename_3D = os.path.join(frames_folder, f"temp_3D_{frame:03d}.png")
    plt.tight_layout()
    plt.savefig(filename_3D, dpi=150, bbox_inches='tight')
    plt.close(fig_3D)

    fig_top = plt.figure(figsize=(5, 2.5), facecolor = 'white')
    ax = fig_top.add_subplot(122)
    sns.heatmap(w2, xticklabels=False, yticklabels=False, cmap='jet', cbar_kws={'label': 'Wave Amplitude', 'ticks': np.linspace(cmin, cmax, 5), 'format': '%.1e'}, vmin=cmin, vmax=cmax)
    
    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=6)
    colorbar.set_label('Wave Amplitude', fontsize=6)
    
    ny, nx = w2.shape
    x_ticks = np.linspace(0, nx-nx%10, int(nx/10) + 1)
    y_ticks = np.linspace(0, ny-ny%40, int(ny/40) + 1)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.tick_params(axis='both', labelsize=6)
    
    ax.set_title(f"Frame {frame + 1} - Top View", fontsize=8)
    ax.set_xlabel("X Axis", fontsize=6, labelpad=8)
    ax.set_ylabel("Y Axis", fontsize=6, labelpad=8)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    filename_top = os.path.join(frames_folder, f"temp_top_{frame:03d}.png")
    plt.tight_layout()
    plt.savefig(filename_top, dpi=150, bbox_inches='tight')
    plt.close(fig_top)
    
    fig_trans = plt.figure(figsize=(4, 2.5), facecolor='white')
    ax_trans = fig_trans.add_subplot(111)
    sns.lineplot(x=energy, y=transmission, ax=ax_trans, color='blue', linewidth=1)
    ax_trans.axvline(x=energy[frame], color='red', linestyle='--', linewidth=1)

    ax_trans.set_title("Transmission vs Energy", fontsize=8)
    ax_trans.set_xlabel("Energy", fontsize=6, labelpad=8)
    ax_trans.set_ylabel("Transmission", fontsize=6, labelpad=8)
    ax_trans.tick_params(axis='x', labelsize=6)
    ax_trans.tick_params(axis='y', labelsize=6)
 
    transmission_filename = os.path.join(frames_folder, f"temp_transmission_{frame:03d}.png")
    plt.tight_layout()
    plt.savefig(transmission_filename, dpi=150, bbox_inches='tight')
    plt.close(fig_trans)

    img_3D = Image.open(filename_3D)
    img_top = Image.open(filename_top)
    img_transmission = Image.open(transmission_filename)

    combined_height = img_3D.height
    aspect_ratio_3D = img_3D.width / img_3D.height
    aspect_ratio_top = img_top.width / img_top.height
    aspect_ratio_transmission = img_transmission.width / img_transmission.height
    
    new_width_3D = int(combined_height * aspect_ratio_3D)
    new_width_top = int(combined_height * aspect_ratio_top)
    new_width_transmission = int(combined_height * aspect_ratio_transmission)
    
    img_3D = img_3D.resize((new_width_3D, combined_height))
    img_top = img_top.resize((new_width_top, combined_height))
    img_transmission = img_transmission.resize((new_width_transmission, combined_height))

    combined_width = img_3D.width + img_top.width + img_transmission.width - 10
    combined_img = Image.new("RGB", (combined_width, combined_height))
    combined_img.paste(img_3D, (0, 0))
    combined_img.paste(img_top, (img_3D.width - 5, 0))
    combined_img.paste(img_transmission, (img_3D.width + img_top.width - 10, 0))

    filename_combined = os.path.join(frames_folder, f"frame_combined_{frame:03d}.png")
    combined_img.save(filename_combined)
    filenames_combined.append(filename_combined)

    os.remove(filename_3D)
    os.remove(filename_top)
    os.remove(transmission_filename)

#GIF
frames = "frames"
output_gif = "wire_electron_density.gif" #adjust name as needed
frame_paths = [os.path.join(frames, f"frame_combined_{str(i).zfill(3)}.png") for i in range(total_frames)]
with imageio.get_writer(output_gif, mode='I', duration=0.1) as writer:
    for frame_path in frame_paths:
        image = imageio.imread(frame_path)
        writer.append_data(image)

#MP4
img = cv2.imread(filenames_combined[0])
height, width, layers = img.shape
video_filename = "wire_electron_density.mp4" #adjust name as needed
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video_writer = cv2.VideoWriter(video_filename, fourcc, 10, (width, height))
for filename in filenames_combined:
    img = cv2.imread(filename)
    video_writer.write(img)
video_writer.release()

#Clean up
for filename in filenames_combined:
    os.remove(filename)
os.rmdir(frames_folder)


# In[ ]:




