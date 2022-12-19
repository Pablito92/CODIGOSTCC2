# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:03:14 2022

@author: pabli

PARA CADA PASSO NO ALGORITMO UMA IMAGEM É GERADA, MOSTRADA E SALVA NO MESMO DIRETÓRIO DE EXECUÇÃO DESTE CÓDIGO.

O ARQUIVO GERADO "eq_gauss_limiar.jpg" TEM UMA REGIÃO REPRESENTADA NA FIGURA 11 NO TEXTO DO TCC.
NO CASO É REPRESENTADO O ANTES E O DEPOIS DA REGIÃO AO SER APLICADA UMA OPERAÇÃO MORFOLÓGICA DE ABERTURA.
"""

from userlib import os, plt, io, figax
from math import pi
from matplotlib import image
import matplotlib.patches as mpatches
from skimage.exposure import equalize_adapthist, histogram
from skimage.feature import blob_dog, blob_doh, blob_log, canny
from skimage.filters import gaussian, threshold_local, threshold_mean, threshold_otsu, threshold_multiotsu
from skimage.color import rgb2gray
from skimage.util import img_as_float
from skimage.measure import label, regionprops
from skimage.morphology import square, disk, star, diamond
from skimage.morphology import remove_small_objects, binary_opening, binary_erosion, binary_dilation
from skimage.segmentation import clear_border

import numpy as np

plt.close('all')
cwd = os.getcwd()
print(cwd +'\\liga.jpg')

image = rgb2gray(image.imread(cwd +'\\liga.jpg'))
plt.figure("image"); plt.imshow(image, cmap='gray')

image_eq = equalize_adapthist(image, kernel_size = 48, clip_limit = 0.01) #, nbins= 256)
plt.figure("image equalized"); plt.imshow(image_eq, cmap='gray')
plt.imsave('eq.jpg', image_eq, cmap='gray')

image_eq = gaussian(image_eq, sigma = 1)
plt.figure("image equalized gaussian"); plt.imshow(image_eq, cmap='gray')
plt.imsave('eq_gauss.jpg', image_eq, cmap='gray')

thresh = threshold_local(image_eq, block_size = 13, method = 'mean')
# thresh = threshold_local(image_eq, 21, 'gaussian')

binary = image_eq > thresh ; binary = 1 - img_as_float(binary)
# plt.figure("image binary threshold local"); plt.imshow(binary, cmap='gray')
plt.imsave('eq_gauss_limiar.jpg', binary, cmap='gray')

binary = binary_opening(binary, disk(5))
plt.figure("image binary opening"); plt.imshow(binary, cmap='gray')
plt.imsave('eq_gauss_limiar_abertura.jpg', binary, cmap='gray')

label_image = label(binary) ; label_image = clear_border(label_image)
#plt.figure("labeled image"); plt.imshow(label_image)

test = label_image > 0
plt.imsave('eq_gauss_limiar_abertura_clearborder.jpg', test, cmap='gray')

regions = regionprops(label_image, image_eq)
for region in regions:
    if region.feret_diameter_max < 15:
        for point in region.coords:
            label_image[point[0]][point[1]] = 0

label_image = label_image > 0
plt.imsave('eq_gauss_limiar_abertura_clearborder_diametro.jpg', label_image, cmap='gray')
plt.imsave('kernelsize48_blocksize13.jpg', label_image, cmap='gray')

label_image = label(label_image)
regions = regionprops(label_image, image_eq)
# plt.figure("labeled image 2"); plt.imshow(label_image)
# pixel_list = []

rdn = 0 ; rdn_mean = 0.00 ; roundness_list = []
sf = 0 ; sf_mean = 0.00 ; sf_list = []
arealist = [] ; minlist = [] ; meanlist = [] ; maxlist = [] ; perimlist = []

count = 0

fig, ax = plt.subplots()
ax.imshow(label_image)
plt.show()
    
for region in regions:
    minl, minc, maxl, maxc = region.bbox    
    largura = maxc - minc
    altura = maxl - minl
    removida = False
    
    sf = region.area / ( pi * ( ( region.feret_diameter_max/2 )**2 ) ) 
        
    if sf < 0.7:
        for point in region.coords:
            label_image[point[0]][point[1]] = 0
        removida = True
    
    if not removida : #se a figura não foi marcada para ser removida
        count = count + 1
        
        rdn = 4 * pi * region.area/(region.perimeter**2)
        rdn_mean = rdn_mean + rdn
        roundness_list.append(rdn)        
        
        sf_mean = sf_mean + sf
        sf_list.append(sf)
        
        arealist.append(region.area)
        minlist.append(region.intensity_min)
        meanlist.append(region.intensity_mean)
        maxlist.append(region.intensity_max)
        
        roundness_list.append(rdn)
        sf_list.append(sf)
        

print()
print('# de regioes segmentadas consideradas nódulos =', count)
print('# de regioes segmentadas totais =', len(regions))
print('Grau de nodularide do material = ', count/len(regions))

fig, ax = plt.subplots()
ax.imshow(label_image)    
plt.show()

test = label_image >0
#SALVA IMAGEM CONTENDO SOMENTE OS NÓDULOS CONSIDERADOS REDONDOS
plt.imsave('eq_gauss_limiar_abertura_clearborder_diametro_final.jpg', test, cmap = 'gray')


