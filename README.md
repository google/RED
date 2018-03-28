# RED: version 1.0.0

Demonstration of the image restoration experiments conducted in 
Y. Romano, M. Elad, and P. Milanfar, "[The Little Engine that Could: 
Regularization by Denoising (RED)](https://epubs.siam.org/doi/10.1137/16M1102884)", 
SIAM Journal on Imaging Sciences, 10(4), 1804–1844, 2017 
[[arXiv](https://arxiv.org/abs/1611.02862)].

![Examples](images/examples.png?raw=true)

The code was tested on Windows 7 and Windows 10, with Matlab R2016a.

This is not an official Google product.


## Quick start

1. Clone or download this repository.
2. Download a lightweight version of TNRD denoising from
   https://drive.google.com/file/d/0B9L0NyTobx_3NHdJLUtPQWJmc2c/view?usp=sharing.
3. Unzip the file tnrd_denoising.zip (containing 7 files in total) to the folder 
   located in `YOUR_LOCAL_PATH/RED/tnrd_denoising`.
4. Download test images from
   https://drive.google.com/drive/folders/0BzCe024Ewz8aZ3pNQXNqYXJseTQ?usp=sharing
   and copy to the folder located in `YOUR_LOCAL_PATH/RED/test_images`.
5. Open Matlab and change the current folder to `YOUR_LOCAL_PATH/RED` (use "cd" function).
6. In Matlab's command window run
   ```
   >> main
   ```
   to start a demo.


## Third-party software components needed

The code uses an image denoising algorithm called TNRD, described in 
Yunjin Chen, and Thomas Pock, "Trainable Nonlinear Reaction Diffusion: 
A Flexible Framework for Fast and Effective Image Restoration", 
IEEE TPAMI 2016. The TNRD code is available in the authors website
http://www.escience.cn/people/chenyunjin
and downloaded from
https://www.dropbox.com/s/8j6b880m6ddxtee/TNRD-Codes.zip?dl=0.

Note1: For a fast execution, please enable "parfor" using Matlab command "matlabpool".

Note2: The TNRD software contains mex files. If you encounter problems when running the
denoiser, try to download the TNRD code from the above dropbox link.
Then, copy the files from the directory 
`YOUR_LOCAL_PATH/TNRD-Codes/TestCodes(denoising-deblocking-SR)/GaussianDenoising`
to `YOUR_LOCAL_PATH/RED/tnrd_denoising` directory.
If you are using mex files for the first time, please run the following in Matlab's
command window:

```
>> mex -setup
You will see the following question:
"Would you like mex to locate installed compilers [y]/n?"
>> y
Now you should choose a compiler. Simply choose the first one by clicking
>> 1
Then verify by
>> y
```

The test images, used to degrade the image, are taken from on the NCSR software. 
For more details, please refer to Weisheng Dong, Lei Zhang, Guangming Shi, and Xin Li 
"Nonlocally Centralized Sparse Representation for Image Restoration", IEEE-TIP, 2013. 
The NCSR code is available in http://www4.comp.polyu.edu.hk/~cslzhang/NCSR.htm.
Also, to have a fair comparison, we use a similar degradation process as done in NCSR.

Good luck!
