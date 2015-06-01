#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import matplotlib.image as mpimg
from PIL import Image
import configparser


#Create PSF shape on a grid based on primary mirror size, input angle, wavelength, position on detector

#Calculate energy per pixel given exposure time,

#Create Zernike function given Zernike coefficients

#Show PSF with wave front errors applied

#Calculate energy per pixel

#Send PSF through fourier optics

#Calculate dpdf

class PSF(object):
    def __init__(self,config_file):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.dia = self.config.getfloat('Instrument','Diam')
        self.det_nx = self.config.getint('Instrument','det_nx')
        self.det_ny = self.config.getint('Instrument','det_ny')

        self.FWHM = self.config.getfloat('Observation','FWHM')
        self.flux = self.config.getfloat('Observation','flux')

        self.x0 = 0 #where the center of PSF is
        self.y0 = 0

         


    def grid(self):
        xx = np.arange(0,self.det_nx,1,dtype='complex')
        yy = np.arange(0,self.det_ny,1,dtype='complex')
        
        x, y = np.meshgrid(xx,yy)
        
        return x, y

psf = PSF('/Users/parkerf/Research/WFIRST/PSF/config_file.txt')

#define an obstruction and np.where around it
def obstruction(imagefile):
    img = mpimg.imread(imagefile)
    #print (img.shape, img.dtype)
    x,y = img.shape[0],img.shape[1]

    #create a bindary file
    if img.dtype == 'uint8':
        img = img[:,:,0]
        img[np.where(img<np.amax(img))] = 1
        img[np.where(img==np.amax(img))] = 0
    else:
        pass

    #ind = np.where(image==1)
    return img


#need to establish where FWHM is defined
def gaussian():
    #Could make this take any obstruction file but easier like this for now.
    obs = obstruction('/Users/parkerf/Research/WFIRST/WFIRST_IFU_Pupil_Mask_for_R1.jpg')
    ox, oy = obs.shape[0],obs.shape[1]

    NX, NY = psf.grid()
    NX, NY = NX-(psf.det_nx/2.), NY-(psf.det_ny/2.)

    alpha = psf.FWHM/2.35
    gauss = psf.flux*((alpha**2)*(2*np.pi))**(-1.)*np.exp(-(1./2.)*((NX-psf.x0)/alpha)**2.)*np.exp(-(1./2.)*((NY-psf.y0)/alpha)**2.)
    
    #cut down on empty space
    intx,inty = np.where(gauss>0)
    #print (np.amin(intx),np.amax(inty))
    if np.amin(intx) > 200 & np.amin(inty) > 200:
        if np.amax(inty)<(len(NY)-200) & np.amax(intx)<(len(NX)-200):
            iminx = np.amin(intx)-200
            imaxx = np.amax(intx)+200
    
            iminy = np.amin(inty)-200
            imaxy = np.amax(inty)+200

            gauss=gauss[iminx:imaxx,iminy:imaxy]
            
    else:
        pass
    #print (iminx,iminy,imaxx,imaxy)

    #Getting the obstruction to be the right size. Wonky I know!
    A = np.zeros((ox/2.+len(NX)/2.,oy/2.+len(NY)/2.))
    A[:ox,:oy]=obs
    B=np.zeros((len(NX),len(NY)))
    B[(len(NX)-len(A[:,0])):,(len(NY)-len(A[0,:])):]=A

    #plt.imshow(B)

    #Avoid the obscuration
    ind = np.where(B>0)
    gauss[ind] = 0

    
    #FFT it!
    FG = fftpack.fft2(gauss)
    FGa = fftpack.fftshift(FG)
    power = np.abs(FGa)**2.
        
    fig = plt.figure(figsize=(6, 3.2))

    ax = fig.add_subplot(311)
    ax.set_title('PSF')
    plt.imshow(power)
    ax = fig.add_subplot(312)
    plt.imshow(np.log(power))
    ax = fig.add_subplot(321)
    plt.imshow(np.abs(gauss))
    plt.show()

#gaussian()
#PSF = MonochromaticSource('/Users/parkerf/Research/WFIRST/PSF/config_file.txt')
#PSF.gaussian()


class WaveFrontErrors(PSF):
    def __init__(self):
        PSF.__init__(self,'/Users/parkerf/Research/WFIRST/PSF/config_file.txt')

        #self.fx = config.getint('Instrument','fx')   #Distance to the image plane
        self.x0 = self.config.getfloat('Observation','x0')  #Angle of entry into telescope
        self.y0 = self.config.getfloat('Observation','y0')  #Angle of entry into telescope
        self.wave = self.config.getfloat('Observation','wave')
        self.n_zernike = 22 #number of zernike coefficients 
    #Not really sure where this was taken from, but got from IDL code.
    #Supposed to be a correction coefficient in front of the zernike_coefficient
   

    def approx_inputs(self):
        #This would approximate input angles and wavelength to align with values given in zernike_coeffs doc
        #These values are specific to the values given by GSFC
        p = np.zeros(3)
        ##wavelength##
        if self.wave >= 1.75:
            p[0] = 2.
        elif 1.75>self.wave >= 1.2:
            p[0] = 1.448
        elif 1.2 > self.wave >= 0.8:
            p[0] = 0.977
        elif 0.8 > self.wave >= 0.7:
            p[0] = 0.728
        elif self.wave < 0.7:
            p[0] = 0.6
        ## x0 ##
        if self.x0 >= 0.0001979175*3.:
            p[1] = 0.00079167
        elif 0.0001979175 > self.x0 >= 0.0001979175*3:
            p[1] = 0.000395835
        elif -0.0001979175 > self.x0 >= 0.0001979175:
            p[1] = 0
        elif -0.0001979175*3 > self.x0 >= -0.0001979175:
            p[1] = -0.000395835
        elif self.x0 < -0.0001979175*3:
            p[1] = -0.00079167
        ## y0 ##
        if self.y0 >= .0022:
            p[2] = 0.002639
        elif 0.0024 > self.y0 >= 0.002:
            p[2] = 0.002222 
        elif 0.002 > self.y0 >= 0.0015:
            p[2] = 0.001806
        elif 0.0015 > self.y0 >= 0.0012:
            p[2] = 0.001389
        elif 0.0012 > self.y0 >= 0.0008:
            p[2] = 0.0009722
        elif 0.0008 > self.y0 >= 0.0001041675*3:
            p[2] = 0.00041667
        elif 0.0001041675*3 > self.y0 >= 0.0001041675:
            p[2] = 0.000208335
        elif 0.0001041675 > self.wave >= -0.0001041675:
            p[2] = 0
        elif -0.0001041675 > self.y0 >= -0.0001041675*3:
            p[2] = -0.000208335
        elif self.y0 < -0.0001041675*3:
            p[2] = -0.00041667

        return p[0],p[1],p[2]
            

    def zernike_coeffs(self):
        wave = self.approx_inputs()[0]
        x0 = self.approx_inputs()[1]
        y0 = self.approx_inputs()[2]

        print ('wave:'+str(wave)+'; x0:'+str(x0)+'; y0:'+str(y0))
        
        z_file = np.genfromtxt('/Users/parkerf/Research/WFIRST/PSF/PSF_zernike_calculation.csv',
                               delimiter=',',skip_header=11)
        z_wave = z_file[:,1]
        z_x0 = z_file[:,4]
        z_y0 = z_file[:,5]
        ind = np.where((z_wave==wave) & (z_x0==x0) & (z_y0==y0))

        coeffs = np.zeros(self.n_zernike)
        for i in range(0,self.n_zernike-1):
            coeffs[i] = z_file[ind,9+i]

        return coeffs

    def zernike_function(self):
        #copied from IDL code essentially
        #What is pas and rpups?

        X, Y = psf.grid()
        X, Y = X-(psf.det_nx/2.), Y-(psf.det_ny/2.)
        ai = self.zernike_coeffs()

        coeff_correction = np.array([1.,4**(1./2),4**(1./2),3**(1./2),6**(1./2),6**(1./2),8**(1./2),
                        8**(1./2),8**(1./2),8**(1./2),5**(1./2),10**(1./2),10**(1./2),
                        10**(1./2),10**(1./2),12**(1./2),10**(1./2),12**(1./2),12**(1./2),
                        12**(1./2),12**(1./2),7**(1./2)])
                        #14**(1./2),14**(1./2),14**(1./2),
                        #14**(1./2),14**(1./2),14**(1./2),16**(1./2),16**(1./2),16**(1./2),
                        #16**(1./2),16**(1./2),16**(1./2),16**(1./2),16**(1./2),9**(1./2)])

        ai = ai*coeff_correction

        #How to calculate these??
        pas = 1.
        rpups = 1.

        p=np.sqrt(((X-psf.x0/2)*pas)**2+((Y-psf.y0/2)*pas)**2)/rpups
        a=np.arctan((Y-psf.y0/2)/(X-psf.x0/2)) 
        p2 = p*p
        p3 = p2*p
        p4 = p3*p
        p5 = p4*p
        p6 = p5*p
        sina = np.sin(a)
        cosa = np.cos(a)
        sin2a = np.sin(2.*a)
        cos2a = np.cos(2.*a)
        sin3a = np.sin(3.*a)
        cos3a = np.cos(3.*a)
        sin4a = np.sin(4.*a)
        cos4a = np.cos(4.*a)
        cos5a = np.cos(5.*a)
        sin5a = np.sin(5.*a)
        cos6a = np.cos(6.*a)
        sin6a = np.sin(6.*a)
        Ze = ai[0] + ai[1]*(p*cosa) + ai[2]*(p*sina) + ai[3]*(2.*p2-1) + \
          ai[4]*(p2*sin2a) + ai[5]*(p2*cos2a) + ai[6]*((3.*p3-2.*p)*sina) + \
          ai[7]*((3*p3-2*p)*cosa) + ai[8]*((p3)*sin3a) + ai[9]*((p3)*cos3a) + \
          ai[10]*((6*p4-6*p2+1)) + ai[11]*((4*p4-3*p2)*cos2a) + ai[12]*((4*p4-3*p2)*sin2a) + \
          ai[13]*((p4)*cos4a) + ai[14]*((p4)*sin4a) + ai[15]*((10*p5-12*p3+3*p)*cosa) + \
          ai[16]*((10*p5-12*p3+3*p)*sina) + ai[17]*((5*p5-4*p3)*cos3a) + ai[18]*((5*p5-4*p3)*sin3a) + \
          ai[19]*((p5)*cos5a) + ai[20]*((p5)*sin5a) + ai[21]*(20*p6-30*p4+12*p2-1)
          #ai[22]*((15*p6-20*p4+6*p2)*sin2a)])
          #+ai[23]*((15*p6-20*p4+6*p2)*cos2a)
          #+ai[24]*((6*p6-5*p4)*sin4a)+ai[25]*((6*p6-5*p4)*cos4a)
          #+ai[26]*((p6)*sin6a)+ai[27]*((p6)*cos6a)))/1000.

        Zfft = fftpack.fft2(Ze)
        Zfft = fftpack.fftshift(Zfft)
        Z_fft = np.abs(Zfft)**2.

        #plot it
        fig = plt.figure(figsize=(6, 3.2))

        ax = fig.add_subplot(111)
        ax.set_title('Wave Front Errors')

        plt.imshow(Z_fft)

        #plt.imshow(Ze+self.gaussian())

        plt.show()
        #return Ze
    
        
WaveFrontErrors().zernike_function()
   


'''
class FourierOptics(WaveFrontErrors):
    def __init__(self):
        WaveFrontErrors.__init__(self)

class dpdf_output(FourierOptics):
    def __init__(self):
        FourierOptics.__init__(self)
        
def make_grid
'''

