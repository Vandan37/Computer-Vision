import numpy as np

def gauss2D(shape=(3, 3),sigma=0.5):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    --This is a code snippet obtained from stackoverflow
    It does the same functioning as OpenCV Gaussian blur
    """
    a, b = [(ss-1.)/2. for ss in shape]
    #print(a, b) --> 14.0, 14.0
    """
    the values of shape in the function argument are given, but
    they are provided in the main (). For this implementation, we use the value as 7
    Hence, the value of ss, for the first cycle of the array will be...
    Going back to cutoff freq equation, shape = cof*4 + 1
    Hence, 'shape' will be 14, 14 for a, b, in the first cycle.  
    """
    i, j = np.ogrid[-a:a+1,-b:b+1]
    #print(i, j)
    """
    i, j represent the column and row resp. FOr this, the values of m will cycle from
    -a to a (not a+1), and -b to b. -14 to 14
    So it will be
    """
    c = np.exp(-(j*j + i*i) / (2.*sigma*sigma))
    """
    What this does is np.exp (-(j[0]*j[0]+i[0]*i[0])
    """
    #print(h, h.shape)--> 29, 29
    c[ c < np.finfo(c.dtype).eps*c.max() ] = 0
    """
    This is a sequence to compare the values of eps*c.max with h. 
    If the value is less, it will be kept the same, but if the value is 
    greater, it will be initialized to 0
    eps- epsfloat numpy function
    The difference between 1.0 and the next smallest representable float larger 
    than 1.0. For example, for 64-bit binary floats in the IEEE-754 standard, eps = 2**-52, approximately 2.22e-16.
    """
    sumc = c.sum()
    # print(sumc)

    if sumc != 0:
        c /= sumc
        "Returning the summed up matrix to the main function"
    return c
