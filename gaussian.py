# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 17:26:01 2016

@author: mears
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from typing import Union

def calc_axes_from_cov(cov,inverse=False,verbose=False):
    # this function calculates the major and minor axes and the rotation angle
    # from the inverse covariance matrix

    if inverse:
        invcov = cov
    else:
        invcov = np.linalg.inv(cov)

    eigenvalues,eigenvectors = np.linalg.eigh(invcov)
    major = 2.355*math.sqrt(1.0/(eigenvalues[0]))
    minor = 2.355*math.sqrt(1.0/(eigenvalues[1]))
    if verbose:
        print('eigenvalues: ',eigenvalues)
        print('eigenvectors: ',eigenvectors)
    phi = 180.0/math.pi*math.atan2(-eigenvectors[0,1],eigenvectors[0,0])
    return major,minor,phi

class Gaussian():
    """Base footprint class"""
    # x = x location of center, relative to point lon0,lat0
    # y = y location of center
    # major = width of major axis
    # minor = wind of minor axis
    # phi = rotation angle.  0 means major axis is on X axis
    
    def __init__(self,
                 _x0:float |int | None=0.0,
                 _y0:float | int | None=0.0,
                 _major:float | int |None=None,
                 _minor:float | int |None=None,
                 _phi:float | int |None=None,
                 invcov:np.ndarray |None=None):
        
        self.x0 = _x0
        self.y0 = _y0
        
        if all(isinstance(item, int | float) for item in [_major,_minor,_phi]):
            # this is the case where the major and minor 3dB sizes and phi
            # are specified
            self.major = _major
            self.minor = _minor
            self.phi = _phi
            self.get_abc()
            self.mu = np.array([self.x0,self.y0])
            self.invcov = np.array([[self.a,self.b],[self.b,self.c]])
            self.cov = np.linalg.inv(self.invcov)
            self.norm = 1.0/(2.0*np.pi*np.sqrt(np.linalg.det(self.cov)))
        elif invcov is not None:
            # this is the case where the inverse covariance matrix is specified
            # [a,b]
            # [b,c]

            self.invcov = invcov
            self.a = invcov[0,0]
            self.b = invcov[0,1]
            self.c = invcov[1,1]
            self.cov = np.linalg.inv(self.invcov)
            self.major,self.minor,self.phi = calc_axes_from_cov(self.invcov,inverse=True)
            self.mu = np.array([self.x0,self.y0])
            self.norm = 1.0/(2.0*np.pi*np.sqrt(np.linalg.det(self.cov)))
        else:
            raise ValueError('Must specify either cov or major,minor,phi')
        
        self.weight = 1.0
        # self.norm is the amplitude that makes the integral equal to 1
        self.gridded_output = None

    def __repr__(self):
        return f'Gaussian({self.x0},{self.y0},{self.major},{self.minor},{self.phi})'

    def get_abc(self):
        sigma_x = self.major/2.355
        sigma_y = self.minor/2.355
        
        cos_azimuth = math.cos((math.pi/180.0)*self.phi)
        sin_azimuth = math.sin((math.pi/180.0)*self.phi)
        sin_2azimuth = math.sin(2.0*(math.pi/180.0)*self.phi)
    
        self.a = 2.0*((cos_azimuth*cos_azimuth)/(2.0*sigma_x*sigma_x) + 
                  (sin_azimuth*sin_azimuth)/(2.0*sigma_y*sigma_y))
        self.b = 2.0*(-1.0*sin_2azimuth/(4.0*sigma_x*sigma_x) + 
                       sin_2azimuth/(4.0*sigma_y*sigma_y))
        self.c = 2.0*((sin_azimuth*sin_azimuth)/(2.0*sigma_x*sigma_x) + 
                  (cos_azimuth*cos_azimuth)/(2.0*sigma_y*sigma_y))
        

        
    def multiply(self,other):

        # this function multiplies two footprints together
        # the result is a new footprint

        result_inv_cov = self.invcov + other.invcov
        result_cov = np.linalg.inv(result_inv_cov)

        result_mu = np.dot(result_cov,(np.dot(self.invcov,self.mu) + 
                                       np.dot(other.invcov,other.mu)))
        
        const = (1.0/(2.0*np.pi))*np.sqrt(
                            np.linalg.det(result_cov)*
                       1.0/(np.linalg.det(self.cov)*
                            np.linalg.det(other.cov)))
        exponent = -0.5*(
            np.dot(np.dot(self.mu,self.invcov),self.mu) +
            np.dot(np.dot(other.mu,other.invcov),other.mu) -
            np.dot(np.dot(result_mu,result_inv_cov),result_mu))
        
        result = Gaussian(result_mu[0],result_mu[1],invcov=result_inv_cov)
        result.weight = const*np.exp(exponent)
        return result
    
    def overlap(self,other):    
        normalization_cov = self.cov + other.cov
        normalization_cov_inv = np.linalg.inv(normalization_cov)
        norm_gaussian = Gaussian(self.x0,self.y0,invcov=normalization_cov_inv)
        return norm_gaussian.output_at(other.x0,other.y0)*norm_gaussian.norm

    def output_at(self, x, y):
        return self.weight*math.exp(0.5*(-1.0*self.a*(x-self.x0)*(x-self.x0) - 
                                     2.0*self.b*(x-self.x0)*(y-self.y0) - 
                                     1.0*self.c*(y-self.y0)*(y-self.y0)))
    
    def gridded(self,x0,y0,dx=1.0,dy=1.0,nx=400,ny=400):

        # this function returns a 2D array of the footprint
        # x0,y0 is the center of the grid
        # dx,dy is the grid spacing
        # nx,ny is the number of points in the grid
        # the grid is returned as a 2D numpy array

        x = np.arange(x0 - (nx/2.0)*dx,x0 + (nx/2.0)*dx,dx)
        y = np.arange(y0 - (ny/2.0)*dy,y0 + (ny/2.0)*dy,dy)
        xgrid,ygrid = np.meshgrid(x,y)
        output = self.weight*np.exp(0.5*(-1.0*self.a*(xgrid-self.x0)*(xgrid-self.x0) - 
                                     2.0*self.b*(xgrid-self.x0)*(ygrid-self.y0) - 
                                     1.0*self.c*(ygrid-self.y0)*(ygrid-self.y0)))
        self.xgrid=xgrid
        self.ygrid=ygrid
        self.dx=dx
        self.dy=dy
        self.gridded_output = output
        
        return output
    
    def gridded_total(self):
        if self.gridded_output is None:
            self.gridded(self.x0,self.y0)
            
        return np.sum(self.gridded_output)*self.dx*self.dy*self.norm
    
    def plot(self,x0,y0,dx,dy,nx,ny,ax_in=None,fig_in=None):

        x = np.arange(x0 - (nx/2.0)*dx,x0 + (nx/2.0)*dx,dx)
        y = np.arange(y0 - (ny/2.0)*dy,y0 + (ny/2.0)*dy,dy)
        xgrid,ygrid = np.meshgrid(x,y)
        values_to_plot = self.gridded(x0,y0,dx,dy,nx,ny)

        if fig_in is None:
            fig = plt.figure(figsize=(10,10))
        else:
            fig = fig_in

        if ax_in is None:
            ax = fig.add_subplot(111)
        else:
            ax = ax_in

        ax.contour(xgrid,ygrid,values_to_plot)
        ax.set_aspect('equal')
        ax.set_xlim(x0 - (nx/2.0)*dx,x0 + (nx/2.0)*dx)
        ax.set_ylim(y0 - (ny/2.0)*dy,y0 + (ny/2.0)*dy)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Gaussian Footprint')

        return fig,ax


if __name__ == "__main__":

    from rss_plotting.plot_2d_array import plot_2d_array

    #g0 = Gaussian(0.0,0.0,100.0,50.0,30.0)
    #g0.calc_axes_from_invcov(verbose=True)

    #fig,ax = g0.plot(0.0,0.0,1.0,1.0,400,400)
    #print(g0.gridded_total())
    
    x0 = 0.0
    y0 = 0.0
    g1 = Gaussian(x0,y0,47.0,39.0,0.0)
    fig,ax = g1.plot(0.0,0.0,1.0,1.0,400,400)
    # print(g1.output_at(0.0,12.5))
    # print(g1.output_at(25.0,0.0))
    # print(g1.gridded_total())

    x0 = 5.0
    y0 = 0.0
    g2 = Gaussian(x0,y0,70.0,70.0,0.0)
    fig,ax = g2.plot(0.0,0.0,1.0,1.0,400,400)

    normalization_cov = g1.cov + g2.cov
    normalization_cov_inv = np.linalg.inv(normalization_cov)
    norm_gaussian = Gaussian(0.0,0.0,invcov=normalization_cov_inv)

    wt_array = np.zeros((101,181))
    for itht,tht in enumerate(np.arange(-90.0,90.1,1.0)):
        for idst,dst in enumerate(np.arange(0.0,100.1,1.0)):
            x0 = dst*math.cos(tht*math.pi/180.0)
            y0 = dst*math.sin(tht*math.pi/180.0)
            wt_array[idst,itht] = norm_gaussian.output_at(x0,y0)*norm_gaussian.norm
        
            g2 = Gaussian(x0,y0,15.0,30.0,45.0)
            g3 = g1.multiply(g2)
            # fig,ax = g1.plot(0.0,0.0,1.0,1.0,400,400)
            # fig,ax = g2.plot(0.0,0.0,1.0,1.0,400,400)
            # fig,ax = g3.plot(0.0,0.0,1.0,1.0,400,400)
            # plt.show()

            #gridded_g1 = g1.gridded(0.0,0.0,1.0,1.0,400,400)*g1.norm
            #gridded_g2 = g2.gridded(0.0,0.0,1.0,1.0,400,400)*g2.norm

    
            #gridded_g3_direct = gridded_g1*gridded_g2
            #print(f'Total overlap at {x0:.1f},{y0:.1f}: {g3.gridded_total():.8f} {np.sum(gridded_g3_direct):.8f} {norm_gaussian.output_at(x0,y0)*norm_gaussian.norm:.8f}')
            #print(f'Total overlap at {x0:.1f},{y0:.1f}: {norm_gaussian.output_at(x0,y0)*norm_gaussian.norm:.8f}')
            #wt_array[iy,ix] = norm_gaussian.output_at(x0,y0)*norm_gaussian.norm

    fig,ax = plot_2d_array(wt_array,
                np.arange(-90.0,90.1,1.0),
                np.arange(0.0,100.1,1.0),
                norm='linear',
                title='Weight',
                xtitle='Relative Azimuth Angle',
                ytitle='Distance (km)',
                cmap='viridis',
                zrange=(0.0,np.nanmax(wt_array)))
    ax.contour(np.arange(-90.0,90.1,1.0), np.arange(0.0,100.1,1.0), wt_array,[np.nanmax(wt_array)/2.0], colors='k',interpolation='none')
    print(norm_gaussian)
    print(norm_gaussian.output_at(0.0,norm_gaussian.minor/2.0))
    print(norm_gaussian.output_at(norm_gaussian.major/2.0,0.0))
    plt.show()
    print('Done')
        
