import numpy as np
def normalize(img):
    """ Function to normalize an input array to 0-1 """
    return (img - np.min(img)) / (np.max(img) - np.min(img))
    """
    Below is the original code, it was used in a different way due to local issues
    """
    # img_min = img.min()
    # img_max = img.max()
    # return (img - img_min) / (img_max - img_min)
