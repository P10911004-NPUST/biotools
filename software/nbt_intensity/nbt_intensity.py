import os
import re
import cv2
import numpy as np
import pandas as pd
import skimage.filters as skf_filters
from PIL import Image, ImageDraw, ImageFont
from urllib.request import urlopen

from filters import *

truetype_url = "https://github.com/P10911004-NPUST/fonts/blob/main/harmonyos-sans/HarmonyOS_Sans_Medium.ttf?raw=true"

# font = {"family": "sans", "weight": "bold", "size": 5}
# matplotlib.rc("font", **font)
# matplotlib.rc("image", cmap="gray")
# matplotlib.use("agg")


def nbt_intensity(img):
    img_name = "NA"
    nbt_area = None
    total_nbt_intensity = None
    nbt_intensity_per_area = None
    raw_img, roi, contours = None, None, None

    try:
        if isinstance(img, str):
            raw_img = cv2.imread(img, cv2.IMREAD_COLOR)[:, :, ::-1]
            img_name = os.path.basename(img)
            img_dir = os.path.dirname(img)

            output_img_name = "OUT_" + img_name
            output_img_dir = os.path.join(
                os.path.dirname(img_dir), "OUT_" + os.path.basename(img_dir)
            )

            if not os.path.exists(output_img_dir):
                os.mkdir(output_img_dir)

        elif isinstance(img, np.ndarray):
            raw_img = img
        else:
            print("Input should be an image directory or a numpy ndarray object")

        h, w, d = raw_img.shape
        _, _, B = cv2.split(raw_img)
        B = cv2.bitwise_not(B)

        try:
            thresh = skf_filters.threshold_multiotsu(B, classes=3)
        except:
            thresh = [0, 0, 255]

        roi = B >= np.max(thresh)

        # The filtering process requires uint8 values
        roi = np.multiply(roi, 255).astype(np.uint8)

        # Filtering process
        for _ in range(5):
            roi = cv2.medianBlur(roi, 7)
            roi = min_filter(roi, (5, 5), iteration=1)
            roi = cv2.medianBlur(roi, 7)
            roi = max_filter(roi, (5, 5), iteration=1)

        # Create contour
        contours, hierarchy = cv2.findContours(
            roi, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )

        # convert back to float values of 0 or 1
        roi = roi / 255.0

        # NBT staining area (pixels)
        nbt_area = np.sum(roi)

        # NBT staining intensity (sum of digital number of the staining area)
        total_nbt_intensity = np.sum((B * 1.0) * roi)

        # Minimize the impact of root size variation to the NBT intensity values
        nbt_intensity_per_area = total_nbt_intensity / nbt_area if nbt_area > 0 else 0

        # Plotting
        ## Draw contour line on a blank image
        contours = cv2.drawContours(B * 0, contours, -1, (255, 0, 0), thickness=3)
        ## Show the coutour as red color
        contours = cv2.merge([contours, contours * 0, contours * 0])
        ## Combined the original image with the contour line
        contours = cv2.addWeighted(raw_img, 1, contours, 1, 0)

        if isinstance(img, str):
            contours = Image.fromarray(contours)
            draw = ImageDraw.Draw(contours)
            font = ImageFont.truetype(urlopen(truetype_url), size=70)
            draw.text(
                (30, 10),
                f"Intensity: {round(total_nbt_intensity / 1_000_000, 4)} M",
                (255, 0, 0),
                font=font,
            )
            contours.save(os.path.join(output_img_dir, output_img_name))
            print(f"Saving to: {img_dir}/{img_name}")
    except:
        print(f"Problematic image: {img_dir}/{img_name}")
        contours = B * 0
        # The image is all white, no staining
        if np.max(B) == 0:
            nbt_area = 0
            total_nbt_intensity = 0
            nbt_intensity_per_area = 0
        else:
            # if the image is problematic
            nbt_area = "NA"
            total_nbt_intensity = "NA"
            nbt_intensity_per_area = "NA"

        if isinstance(img, str):
            contours = Image.fromarray(contours)
            draw = ImageDraw.Draw(contours)
            font = ImageFont.truetype(urlopen(truetype_url), size=70)
            draw.text(
                (30, 10),
                f"Intensity: {round(total_nbt_intensity / 1_000_000, 4)} M",
                (255, 0, 0),
                font=font,
            )
            contours.save(os.path.join(output_img_dir, output_img_name))

    return img_name, nbt_area, total_nbt_intensity, nbt_intensity_per_area


# def nbt_intensity_multiproc(img_list):
#    if len(img_list) == 0 or img_list is None:
#        pass

#    img_dir = os.path.dirname(img_list[0])

#    img_list = [
#        os.path.join(img_dir, i)
#        for i in img_list
#        if i.lower().endswith((".jpg", ".jpeg", "tif", "tiff", "png", ".bmp"))
#    ]

#    output_dir = os.path.join(
#        os.path.dirname(img_dir),
#        "OUT_" + os.path.basename(img_dir),
#    )

#    df0 = {
#        "img_name": [],
#        "nbt_area": [],
#        "nbt_intensity": [],
#        "nbt_intensity_per_area": [],
#    }

#    pool = mp.Pool()
#    res = pool.map(nbt_intensity, img_list)

#    df0["img_name"] = [i[0] for i in res]
#    df0["nbt_area"] = [i[1] for i in res]
#    df0["nbt_intensity"] = [i[2] for i in res]
#    df0["nbt_intensity_per_area"] = [i[3] for i in res]
#    raw_img = [i[4] for i in res]
#    roi = [i[5] for i in res]
#    contours = [i[6] for i in res]

#    df0 = pd.DataFrame.from_dict(df0)
#    df0.to_csv(
#        os.path.join(output_dir, os.path.basename(output_dir) + ".csv"), index=False
#    )
#    print(f"Save to {os.path.join(output_dir, os.path.basename(output_dir) + '.csv')}")


## def nbt_intensity_multithread(img_list):
##    if len(img_list) == 0 or img_list is None:
##        pass
##    # input_folder_dir = folder_path.get()
##    input_folder_dir = os.path.dirname(img_list[0])
##    img_list = [
##        os.path.join(input_folder_dir, i)
##        for i in os.listdir(input_folder_dir)
##        if i.lower().endswith((".jpg", ".jpeg", "tif", "tiff", "png", ".bmp"))
##    ]

##    output_dir = os.path.join(
##        os.path.dirname(input_folder_dir),
##        "OUT_" + os.path.basename(input_folder_dir),
##    )

##    df0 = {
##        "img_name": [],
##        "nbt_area": [],
##        "nbt_intensity": [],
##        "nbt_intensity_per_area": [],
##    }

##    with ThreadPoolExecutor() as executor:
##        res = executor.map(nbt_intensity, img_list)

##    res = [i for i in res]

##    df0["img_name"] = [i[0] for i in res]
##    df0["nbt_area"] = [i[1] for i in res]
##    df0["nbt_intensity"] = [i[2] for i in res]
##    df0["nbt_intensity_per_area"] = [i[3] for i in res]
##    df0 = pd.DataFrame.from_dict(df0)
##    df0.to_csv(
##        os.path.join(output_dir, os.path.basename(output_dir) + ".csv"), index=False
##    )
##    print(f"Save to {os.path.join(output_dir, os.path.basename(output_dir) + '.csv')}")
