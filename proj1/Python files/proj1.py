

import os
#import scipy
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.image as mtpimg
"""importing the math, matlab libraries. Matlab library used to display the results"""
#from PIL import Img

from my_imfilter import my_imgfilter
from vis_hybrid_image import vis_hybrid_image
from normalize import normalize
from gauss2D import gauss2D
"""importing other python files for use in the proj1 main funtion"""
def main():
    # function to create hybrid imgs
    main_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    #main_path = '/Users/vvcj2/Desktop/proj1/Proj1'
    # __file__returns the pathname from which the file was loaded, abspath() is used to covert it to an absolute path/
    img1 = mtpimg.imread(os.path.join(main_path, 'data', 'plane.bmp'))
    img2 = mtpimg.imread(os.path.join(main_path, 'data', 'bird.bmp'))
    img1 = img1.astype(np.float32)/255
    img2 = img2.astype(np.float32)/255
    """
    img1 will be converted to a low frequency img, i.e. it will be blurred.
    img2 will be converted to a high frequency img, i.e. it will be sharpened"""

    # Applying filter to both the imgs, constructing a hybrid img.
    cutoff_frequency_1 = 5
    """cutppf_frequency is the standard deivation for the gaussian blur
    Upon trying different values for the cutoff_frequency, the higher the value,
    the more blurry the low pass img becomes. Choosing a value lower than 6, 
    the high pass img can be seen but not to an extent observable. Hence, choosing 7.
    The best result for me was taking 6.5, but keeping it an integer.  
    """
    gaussian_filter_1 = gauss2D(shape=(cutoff_frequency_1*4+1,cutoff_frequency_1*4+1),
                                sigma = cutoff_frequency_1)
    #We were told to use sigma as 1, but using a greater value to get better results
    #Calling the gauss2D funtion as gaussian_filer_1 from gauss2D.py, refer to the same for the function

    low_frequencies = my_imgfilter(img1, gaussian_filter_1)
    #creating a low frequency img, applying the gauss2D the img1
    ############################################################################
    # Remove the low frequencies from img2. The easiest way to do this is to #
    # subtract a blurred version of img2 from the original version of img2.#
    # This will give you an img centered at zero with negative values.       #
    ############################################################################
    cutoff_frequency_2 = 5
    gaussian_filter_2 = gauss2D(shape=(cutoff_frequency_2*4+1,cutoff_frequency_2*4+1),
                                sigma = 1)
    """
    The cutoff frequencies are taken as different variables, but taken same values
    A higher value for cutoff_frequency_2 means the img will be more sharper, and it will
    supercede the blur img. The final result will have more of a sharpened look, with only the colours 
    of the blur img being visible. 
    """
    """
    Since all the imgs had differnt results when blurred and sharpened, different values are used for each set
    Dog-Cat pair (7, 7)
    Einstein-marilyn (3, 3), etc
    """
    low_frequencies_2 = my_imgfilter(img2, gaussian_filter_2)
    high_frequencies = img2 - low_frequencies_2
    #Applying filter to img2, creating a high frequency img by subtraction
    high_frequencies=normalize(high_frequencies)
    #Calling the normalize function from normalize.py, refer to the same for the function

    # print(np.min(low_frequencies))
    # print(np.max(low_frequencies))
    #Un-comment the preceding statements to see the normalized values for the low frequency imgs

    #Combining the high frequencies and low frequencies

    hybrid_img = low_frequencies + high_frequencies
    hybrid_img= normalize(hybrid_img)
    #Creating the hybrid img, normalizing it

    # Visualize and save outputs
    plot.figure(1)
    plot.imshow(low_frequencies)
    plot.figure(2)
    plot.imshow(high_frequencies+0.5)
    plot.figure(3)
    plot.imshow(hybrid_img)
    #Plotting the matrices using matplot library

    vis = vis_hybrid_image(hybrid_img)
    plot.figure(4)
    plot.imshow(vis)
    #Calling the vis_hybrid_img from vis_hybrid_img.py, refer to the same for the function

    plot.imsave(os.path.join(main_path, 'results', 'low_frequencies.png'), low_frequencies, dpi=95)
    plot.imsave(os.path.join(main_path, 'results', 'high_frequencies.png'), high_frequencies, dpi=95)
    plot.imsave(os.path.join(main_path, 'results', 'hybrid_img.png'), hybrid_img, dpi=95)
    plot.imsave(os.path.join(main_path, 'results', 'hybrid_img_scales.png'), vis, dpi=95)
    #Saving the results in the main path/results, 'results' joined
    plot.show()

if __name__ == '__main__':
    main()
