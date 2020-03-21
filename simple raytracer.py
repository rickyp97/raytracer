##simple raytracer v1##

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as spc
import random as rand

# ----------------------------- classes -------------------------------


class Ray2D(object):
    
    def __init__(self, x = 0, y = 0, theta = 0, lamda = 560, time_delay = 0): 
        if theta > 90 or theta < -90:
            raise ValueError("angle must be between -90 and +90 (deg), wrt x axis")
        self.y = y
        self.x = x
        self.angle = theta
        self.wavelenght = lamda
        self.m = np.tan(theta/180*np.pi)
        self.c = y
        self.positions_x = [self.x]
        self.positions_y = [self.y]
        self.time = time_delay
        return
    
    
class Collimated_beam(object):
    
    def __init__(self, rays_number = 10, y = 0, theta = 0, lamda0 = 560, time_width = 0, bandwidth = 100, shape = "gaussian" ): #bandwidth and time width are FWHM
        self.rays = []
        self.width = time_width
        self.shape = shape
        
        if self.shape == "gaussian":
            counter = 0
            while counter < rays_number:
                counter = counter + 1
                t = rand.gauss(0,time_width/2.35482)
                self.rays.append(Ray2D(y, theta, lamda = rand.gauss(lamda0,bandwidth/2.3548),time_delay = t))
            return
        else:
            counter = 0
            while counter < rays_number:
                counter = counter + 1
                self.rays.append(Ray2D(y, theta, lamda = np.random.randint(lamda0-bandwidth/2,lamda0+bandwidth/2)))
            return
    
    def plot_rays (self, color = "visible"):
        plt.figure("rays")
        for ray in self.rays:
            if color == "visible":
                plt.plot(ray.positions_x,ray.positions_y, color = wl2col(ray.wavelenght))
            else:
                plt.plot(ray.positions_x,ray.positions_y, color = color)
        return
    
    def time_diff(self):
        for ray in self.rays:
            path = calculate_path(ray)
            ray.time = ray.time + path/spc.c
        return   
    
    def propagate(self,elements):
        for element in elements:
            for ray in self.rays:
                element.propagate(ray)
        return
      
    def plot_chirp(self,color = 'b', figure = 'chirp'):
        plt.figure(figure)
        for ray in self.rays:
            plt.plot(ray.wavelenght, ray.time, c = color, marker = '.')
        return
    
    def plot_time_shape(self,color = 'b', figure = 'initial shape', bins = 100):
        plt.figure(figure)
        time = []
        for ray in self.rays:
            time.append(ray.time)  
        plt.hist(time,bins = bins, density = False)
        plt.xlabel('Time(s)')
        plt.ylabel("incident rays count")
        print (np.std(time)*2.3582)
        return
    
    def plot_shapes(self,color = 'b', figure = 'initial shape', bins = 100):
        plt.figure(figure)
        time = []
        wavelength = []
        for ray in self.rays:
            time.append(ray.time)
            wavelength.append(ray.wavelenght)
        plt.subplot(121,xlabel = 'Time(s)',ylabel = "incident rays count")    
        plt.hist(time,bins = bins, density = False)
        plt.subplot(122,xlabel = 'Wavelength(nm)',ylabel = "incident rays count")  
        plt.hist(wavelength,bins = bins, density = False)
        print (np.std(time)*2.3582)
        print (np.std(wavelength)*2.3582)
        return
            
    
    
class Grating(object):

    def __init__(self, x = 20, y = 0, angle = 45, gratingsep = 0.00004, facing = "up", order = 1):  # angle is between grating and x
        if facing != "up" and facing != "down":
            raise ValueError("facing must be 'up' or 'down'")    
        if angle > 90 or angle < -90:
            raise ValueError("angle must be between -90 and +90 (deg), wrt x axis")
        self.n_dir = facing
        self.y = y
        self.x = x
        self.angle = angle
        self.m = np.tan(angle/180*np.pi)       
        self.c = self.y-self.m*self.x
        self.d = gratingsep
        self.order = order
        return
        
    def propagate(self,ray): #devide in 4 cases cause i'm bad and can't find general relationship
        if self.n_dir == "up" and ray.angle < self.angle:       #grating facing up, ray incoming from left
            incident_angle = 90+ray.angle-self.angle
            out_angle = 180/np.pi*np.arcsin(np.sin(incident_angle*np.pi/180)-(self.order*ray.wavelenght*10**-9/self.d))
            out_wrt_x_axis = 90+self.angle-out_angle
            
        elif self.n_dir == "up" and ray.angle > self.angle:       #grating facing up, ray incoming from right
            incident_angle = 90+self.angle-ray.angle
            out_angle = 180/np.pi*np.arcsin(np.sin(incident_angle*np.pi/180)-(self.order*ray.wavelenght*10**-9/self.d))
            out_wrt_x_axis = out_angle-(90-self.angle)
            
        elif self.n_dir == "down" and ray.angle < self.angle:    #grating facing down, ray incoming from right
             incident_angle = 90+ray.angle-self.angle
             out_angle = 180/np.pi*np.arcsin(np.sin(incident_angle*np.pi/180)-(self.order*ray.wavelenght*10**-9/self.d))
             out_wrt_x_axis = 90+self.angle-out_angle
        
        elif self.n_dir == "down" and ray.angle > self.angle:  #grating facing down and ray coming from left
             incident_angle = 90+self.angle-ray.angle
             out_angle = 180/np.pi*np.arcsin(np.sin(incident_angle*np.pi/180)-(self.order*ray.wavelenght*10**-9/self.d))
             out_wrt_x_axis = out_angle-(90-self.angle)
        
        if out_wrt_x_axis > 90:
            out_wrt_x_axis = out_wrt_x_axis - 180
        elif out_wrt_x_axis < -90:
            out_wrt_x_axis = out_wrt_x_axis + 180
        if np.abs(out_wrt_x_axis) == 90 and self.n_dir == "down":
            out_wrt_x_axis = -90
        elif np.abs(out_wrt_x_axis) == 90 and self.n_dir == "up":
            out_wrt_x_axis = 90
            

#        print (incident_angle)
#        print(out_wrt_x_axis)
        
        ray.angle = out_wrt_x_axis
        ray.x,ray.y = self.intercept(ray)
        ray.m = np.tan(ray.angle/180*np.pi)
        ray.c = ray.y-ray.m*ray.x
        ray.positions_x.append(ray.x)
        ray.positions_y.append(ray.y)
        return
    
    def intercept(self,ray):
        if np.absolute(ray.angle) == 90:
            x_int = ray.x
            y_int = self.m*ray.x+self.c
        else:
            x_int = (ray.c-self.c)/(-ray.m+self.m)
            y_int = (-ray.m*self.c+self.m*ray.c)/((-ray.m+self.m))
        return x_int,y_int



class Verticalmirror(object):
    
    def __init__(self, x = 3):
        self.x = x
        return
    
    def propagate(self,ray):
        ray.x = self.x
        ray.y = ray.m*self.x+ray.c 
        ray.positions_x.append(ray.x)
        ray.positions_y.append(ray.y)
        ray.m = -ray.m
        ray.c = ray.y-ray.m*ray.x
        return
        
    
class Verticalscreen:
    
    def __init__(self,x = 40):
        self.x = x
        return
    
    def propagate(self,ray):
        ray.x = self.x
        ray.y = ray.m*self.x+ray.c
        ray.positions_x.append(ray.x)
        ray.positions_y.append(ray.y)
        return
    
    
    
class Thinlens(object):               #doesn't work for horizontal lenses yet
    
    def __init__(self, x = 0.5, y = 1, angle = 0, f = 1):  # angle of optical axis of lens with x axis. betweem -90 and 90
        self.x = x
        self.y = y
        self.theta = angle
        if angle > 0:
            self.angle = angle - 90
        elif angle <= 0:
            self.angle = 90+angle
        self.m = np.tan(self.angle/180*np.pi)       
        self.c = self.y-self.m*self.x
        self.focal = f
        self.dist = 0
        return
    
    
    def propagate(self, ray):
        ray.x,ray.y = self.intercept(ray)
        incident_angle = (ray.angle - self.theta)*np.pi/180
        if ray.y >= self.y:
            self.dist = np.sqrt((ray.x-self.x)**2+(ray.y-self.y)**2)
        elif ray.y < self.y:
            self.dist = -np.sqrt((ray.x-self.x)**2+(ray.y-self.y)**2)    
        out_angle = (-1/self.focal*(self.dist)+(incident_angle)) *180/ np.pi
        ray.angle = out_angle+self.theta
        ray.m = np.tan(ray.angle/180*np.pi)
        ray.c = ray.y-ray.m*ray.x
        ray.positions_x.append(ray.x)
        ray.positions_y.append(ray.y)
        return 
    
    def intercept(self,ray):
        if np.absolute(self.angle) == 90:
            x_int = self.x
            y_int = ray.m*self.x+ray.c
        else:
            x_int = (ray.c-self.c)/(-ray.m+self.m)
            y_int = (-ray.m*self.c+self.m*ray.c)/((-ray.m+self.m))
        return x_int,y_int


#---------------- functions -----------------------------------------

def calculate_path(ray):
    path = 0
    for i in range(len(ray.positions_x)-1):
        l = np.sqrt((ray.positions_x[i+1]-ray.positions_x[i])**2+(ray.positions_y[i+1]-ray.positions_y[i])**2)
        path = path+l
    return path

def find_x_given_y(elements, y, ray_y = 0, ray_theta = 0, ray_lamda = 560):      #find coords of optical axis after elements
    trial_ray = Ray2D(ray_y,ray_theta, lamda = ray_lamda)
    for element in elements:
        element.propagate(trial_ray)
    x_pos = (y-trial_ray.c)/trial_ray.m
    return x_pos

def find_y_given_x(elements, x, ray_y = 0, ray_theta = 0, ray_lamda = 560):      #find coords of optical axis after elements
    trial_ray = Ray2D(ray_y,ray_theta, lamda = ray_lamda)
    for element in elements:
        element.propagate(trial_ray)
    y_pos = x*trial_ray.m+trial_ray.c
    return y_pos

def path_between_elements(all_elements, elem_for_path, ray_y = 0, ray_theta = 0, ray_lamda = 560):  #calculate lenght of path (through elements in elem_for_path) , used to then calculate time
     trial_ray = Ray2D(ray_y,ray_theta, lamda = ray_lamda)
     for element in all_elements:
         element.propagate(trial_ray)
     path = 0
     index_1 = all_elements.index(elem_for_path[0])
     index_2 = all_elements.index(elem_for_path[1])
     for i in range(index_1+1,index_2+1):
         l = np.sqrt((trial_ray.positions_x[i+1]-trial_ray.positions_x[i])**2+(trial_ray.positions_y[i+1]-trial_ray.positions_y[i])**2)
         path = path+l
     return path


    
def wl2col(w):  #convert wavelenght into rgb color for plots
    
    if w >= 380 and w < 440:
        R = -(w - 440) / (440 - 350)
        G = 0
        B = 1
    elif w >= 440 and w < 490:
        R = 0
        G = (w - 440) / (490 - 440)
        B = 1
    elif w >= 490 and w < 510:
        R = 0
        G = 1
        B = -(w - 510) / (510 - 490)
    elif w >= 510 and w < 580:
        R = (w - 510) / (580 - 510)
        G = 1
        B = 0
    elif w >= 580 and w < 645:
        R = 1
        G = -(w - 645) / (645 - 580)
        B = 0
    elif w >= 645 and w <= 780:
        R = 1
        G = 0
        B = 0
    else:                   # this makes everything out of the visible range black
        R = 0
        G = 0
        B = 0

    # intensity correction, colors fade at the edge of the visible range
    if w >= 380 and w < 420:
        SSS = 0.3 + 0.7*(w - 350) / (420 - 350)
    elif w >= 420 and w <= 700:
        SSS = 1
    elif w > 700 and w <= 780:
        SSS = 0.3 + 0.7*(780 - w) / (780 - 700)
    else:
        SSS = 0
    

    return ((SSS*R), (SSS*G), (SSS*B))



###########################code here to save time #########################################

if __name__ == "__main__":

   
#------------ null stretcher? ----------------    
#    beam = Collimated_beam(2000, time_width = 200E-15)
#    
#    g = Grating(1, 0, 80, 5.75E-5, order = -1)  
#    l1 = Thinlens( find_x_given_y([g],0.35), y =0.35, angle = -20, f = -1 )
#    l2 = Thinlens( find_x_given_y([g],1.06), y =1.06, angle = -20 , f = -1)
#    g2 = Grating(find_x_given_y([g],1.42),1.42,80, 5.75E-5, facing = "down", order = 1)
#    s1 =Verticalscreen(1)
#    beam.plot_chirp('g','initialchirp')
#    beam.plot_shapes(figure = "initial shape")
#    beam.propagate([g,l1,l2,g2,s1])
#    beam.plot_rays()
#    beam.time_diff() 
#    beam.plot_chirp('g','chirp3')
#    beam.plot_time_shape(figure = "final shape")
    
# -------------------- martinez stretcher ---------------------------    
    
    beam = Collimated_beam(2000, lamda0 = 1053, time_width = 250E-15, bandwidth = 100)
    
    m1 = Grating(1, 0, 70, 5.75E-5, order = 0)  
    g1 = Grating(0, find_y_given_x([m1],0, ray_lamda = 1053), 69.44165, 5.75E-5, facing = "down", order = 1)
    l1 = Thinlens(1, find_y_given_x([m1,g1], 1, ray_lamda = 1053), angle = 0, f = 1 )
    l2 = Thinlens(3, find_y_given_x([m1,g1,l1], 3, ray_lamda = 1053), angle = 0 , f = 1)
    g2 = Grating(3.5, find_y_given_x([m1,g1,l1,l2],3.5, ray_lamda = 1053),  -69.44165, 5.75E-5, facing = "down", order = 1)
    m2 = Grating(find_x_given_y([m1,g1,l1,l2,g2],0, ray_lamda = 1053), 0, -50, 5.75E-5, order = 0) 
    l3 = Thinlens(3, find_y_given_x([m1,g1,l1,l2,g2,m2,g2], 3, ray_lamda = 1053), angle = 0 , f = -1)
    l4 = Thinlens(1, find_y_given_x([m1,g1,l1,l2,g2,m2,g2,l3], 1, ray_lamda = 1053), angle = 0, f = -1 )
    v = Verticalscreen(-1)
    
    beam.plot_chirp('g','initialchirp')
    beam.plot_shapes(figure = "initial shape")
    beam.propagate([m1,g1,l1,l2,g2,m2,g2,l3,l4,g1,m1,v])
    print(beam.rays[12].angle)
    beam.plot_rays()
    beam.time_diff() 
    beam.plot_chirp('g','chirp3')
    beam.plot_time_shape(figure = "final shape")
    
# -------------------- cerberus parameters (not working atm) --------------------------
#
#    beam = Collimated_beam(1000, lamda0 = 1053, time_width = 250E-15, bandwidth = 10)
#    g1 = Grating(0.003,0,10, 5.743E-7 , order = 1)
#    g2 = Grating(0,find_y_given_x([g1],0,ray_lamda=1053), 10, 5.743E-7,facing = "down", order = 1)
#    m1 = Verticalmirror(x = 0.02)
#    g3 = Grating(-1, 0, -20, 5.743E-7, order = 1)
#    g4 = Grating(0, find_y_given_x([g1,g2,m1,g2,g1,g3],0, ray_lamda = 1053), -20, 5.743E-7, facing = "down", order = 1)
#    m2 = Verticalmirror(-0.5)
#    
#    g5 = Grating(0, 0, -30, 5.743E-7, facing = "down", order = 1)
#    m3 = Grating(find_x_given_y([g5],-0.3,ray_lamda = 1053), -0.3, -7.6775, 5.743E-7, facing = "up", order = 0)
#    l1 = Thinlens(1.5, find_y_given_x([g5,m3], 1.5, ray_lamda = 1053), angle = 0, f = 1.5 )
#    l2 = Thinlens(4.5, find_y_given_x([g5,m3,l1], 4.5, ray_lamda = 1053), angle = 0, f = 1.5 )
#
#    
#    beam.plot_chirp('g','initialchirp')
#    beam.plot_shapes(figure = "initial shape")
#    beam.propagate([g1,g2,m1,g2,g1,g3,g4,m2,g4,g3,g5,m3,l1,l2,g5])
#    print (beam.rays[20].angle)
#    beam.plot_rays()
#    beam.time_diff() 
#    beam.plot_chirp('g','first chirp')
#    beam.plot_time_shape(figure = "first compressor")


