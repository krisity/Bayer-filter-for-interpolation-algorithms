# Bayer-filter-for-interpolation-algorithms

Main task: 

**Perform the following algorithms using a Bayer filter:<br />**

♦ bilinear interpolation;<br />
♦ median interpolation (Freeman);<br />
♦ interpolation using variable number of gradients (Chang);<br />
♦ gradient interpolation (Laroche-Prescott);<br />
♦ adaptive interpolatiob (Hamilton-Adams).<br />

*The Bayer filter that was used has the following pattern:*

![](Screenshots/GRBG1.JPG)

*Screenshots of the program:*

- Main screen that appears on opening:
![](Screenshots/1.JPG)
- Dialog window that appears when the user clicks on the "Browse" button:
![](Screenshots/2.JPG)
- Selected image appears on the left side of the screen:
![](Screenshots/3.JPG)
- Interpolated image appears on the right after the user clicks on one of the interpolation algorithms (each algorithm is displayed with a button):
![](Screenshots/4.JPG)
- MSE/PSNR calculation appear on the right side after each interpolation on the image:
![](Screenshots/4.1.JPG)
                                                                                 
- The user can click on the "CLEAR" button on the right to remove the interpolated image. After that the user can perform another interpolation on the same image that was previously uploaded:

![](Screenshots/5 - clear button.JPG)

- Dialog window that appears after the user clicks on the "SAVE" button:
![](Screenshots/6.JPG)

- The user can zoom on part of the image:
![](Screenshots/7 Zoom.JPG)

- Zoomed part of the original image (left) and interpolated image (right). Interpolation errors can be seen clearly.
![](Screenshots/7 Zoom 2.JPG)

