# -*- coding: utf-8 -*-
"""Edge detection.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/10ZIvyVgDjGlWd09LJZzwboQs-4RPlCut

## Edge detection in OpenCV and skimage
"""

!wget https://pns2019.github.io/images/Lenna.png

"""### Edge detection in OpenCV"""

import cv2
import numpy as np
from matplotlib import pyplot as plt

# read image
img = cv2.imread("Lenna.png", 0)
# Find edge with Canny edge detection
edges = cv2.Canny(img, 100, 200)

# display results
plt.subplot(121), plt.imshow(img, cmap='gray')
plt.title('Original Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122), plt.imshow(edges, cmap='gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

plt.show()

"""### Edge detection in skimage"""

import numpy as np
from skimage.io import imread
from skimage.feature import canny

import matplotlib.pyplot as plt

# read image
img = imread("Lenna.png", as_grey=True)

# find edge with Canny edge detection
edges = canny(img)

# display results
plt.subplot(121), plt.imshow(img, cmap='gray')
plt.title('Original Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122), plt.imshow(edges, cmap='gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

plt.show()