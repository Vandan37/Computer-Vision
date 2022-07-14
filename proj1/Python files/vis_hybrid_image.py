import numpy as np
import skimage
from scipy.misc import imresize

def vis_hybrid_image(hybrid_img):
    """Observe hybrid image by concatenation in various scales"""

    scales_num = 5
    # Number of downsampled versions to create

    # By how much factor to downsamplee
    padding = 5
    # how many pixels to pad.

    original_height = hybrid_img.shape[0]
    num_colors = hybrid_img.shape[2]
    #counting color channels in the input
    output = hybrid_img[:]
    current_image = hybrid_img[:]

    for i in range(1, scales_num):
        #add padding
        output = np.concatenate((output,
                                 np.ones((original_height, padding, num_colors))),
                                axis=1)
        #Dowsampling the image;

        size = np.array(current_image.shape)
        # resize the image to half of it's original size
        newsize = (size[:2]*0.5).astype(float)
        # new rwsized iamge
        current_image = skimage.transform.resize(current_image, newsize)
        #pad the top and append to the output
        tmp = np.concatenate((np.ones((original_height-current_image.shape[0],
                                       current_image.shape[1],num_colors)), current_image), axis=0)
        output = np.concatenate((output,tmp), axis=1)

    return output