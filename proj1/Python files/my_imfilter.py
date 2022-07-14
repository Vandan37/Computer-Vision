import numpy as np

def my_imgfilter(image, imgfilter):

    R = image[:,:,0]
    G = image[:,:,1]
    B = image[:,:,2]
    """
    The ":" is used to get all values. Image[:, :, 0] means get all rows,
     all columns, and the first color channel; index 0, which is the red channel
    This converts all the images to greyscale, and stores them as greyscale
    The value is not normalized or something. THe grey value will be corresponding to the 
    value of Red in the color image
    """
    H_image = R.shape[0]
    #print --> 360, for cut off freq-7

    W_image = R.shape[1]
    #print --> 410 for cof - 7
    "For the first case"
    H_filt = imgfilter.shape[0]
    #print --> 29
    W_filt = imgfilter.shape[1]
    # print --> 29
    H_padding = int((H_filt-1)/2)
    # print --> 14
    W_padding = int((W_filt-1)/2)
    # print --> 14

    npad = ((H_padding, H_padding), (W_padding, W_padding))
    RGB_padded = []
    #Initializing a new array RGB_pad
    RGB_padded.append(np.pad(R, pad_width=npad, mode='reflect'))
    # print((np.pad(R, pad_width=npad, mode='reflect').shape) will give out the size (389, 438)
    #'reflect' is used to pad with the reflection of the vector mirrored on the first
    # and last values of vector along each axis
    RGB_padded.append(np.pad(G, pad_width=npad, mode='reflect'))
    RGB_padded.append(np.pad(B, pad_width=npad, mode='reflect'))
    #Appending RGB_pad with RGB values. THere will be 3 elements, and each with a shape of (389, 438)
    output = np.zeros_like(R)
    # print(np.zeros_like(R))
    # numpy Zeros_like returns an array of zeros with the same shape and type as a given array

    for each in RGB_padded:
        RGB_padded_new = []

        for m in range(H_image):
            # convolution of the whole matrix
            for n in range(W_image):
                # convolution of each small matrix
                total = 0
                total = np.sum(np.multiply(each[m:m+H_filt, n:n+W_filt], imgfilter))
                #multiplying each respective elemen, summing them
                RGB_padded_new.append(total)

        RGB_padded_new = np.asarray(RGB_padded_new)
        #Appending to the new RGB array
        RGB_padded_new = RGB_padded_new.reshape(H_image, W_image)


        """Reshaping to a new shape of H_img and W_img, using a numpy.reshape function
        numpy.reshape(a, newshape, order='C')
        """

        "This is used to combine RGB channel into 3D array"
        output = np.dstack((output, RGB_padded_new))

    output = output[:, :, 1:]
    "Removinh the zeros array"
    return output