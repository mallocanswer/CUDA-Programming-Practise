# CUDA PROGRAMMING PRACTISE
This repository is intended for cuda programming practise.               
I just find thinking in parallel has a lot of fun and wish to do more pracises to refine my skills :)
## Image Blurring
In image processing, a Gaussian blur (also known as Gaussian smoothing) is the result of blurring an image by a Gaussian function. The visual effect of this blurring technique is a smooth blur resembling that of viewing the image through a translucent screen, distinctly different from the bokeh effect produced by an out-of-focus lens or the shadow of an object under usual illumination.       
In this practise, first separate the different color channels so that each color is stored contiguously instead of being interleaved, then compute a blurred pixel value by multiplying all neighbor pixels with their corresponding weights. 
## HDR Tone-mapping
Tone mapping is a technique used in image processing to map one set of colors to another to approximate the appearance of high dynamic range images in a medium that has a more limited dynamic range. Tone mapping addresses the problem of strong contrast reduction from the scene radiance to the displayable range while preserving the image details and color appearance important to appreciate the original scene content.     
In this practise, first adopt reduce function to compute the maximum and minimum value of input image, then use histogram and scan to compute thecumulative distribution of the luminance values.
