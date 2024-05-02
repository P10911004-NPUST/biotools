import cv2
import math

def min_filter(img, kernel=(3, 3), iteration=1):
   shape = cv2.MORPH_RECT
   kernel = cv2.getStructuringElement(shape, kernel)
   res = cv2.erode(img, kernel, iterations=iteration)
   return res

def max_filter(img, kernel_size=(3, 3), iteration=1):
   shape = cv2.MORPH_RECT
   kernel = cv2.getStructuringElement(shape, kernel_size)
   res = cv2.dilate(img, kernel)
   return res

def sd_filter(img, kernel_size=(3, 3)):
   img1 = cv2.blur(img ** 2, kernel_size)
   img2 = cv2.blur(img, kernel_size) ** 2
   res = math.sqrt(img1 - img2)
   return res