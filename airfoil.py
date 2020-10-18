import numpy as np
import matplotlib.pyplot as plt
import math

def naca4(name,grid_points,closed=True):
    '''
    Parameters
    ----------
    name : String
        name of the airfoil
    grid_points : int
        number of points to consider on the X-axis
    closed : bool, optional
        To define if the trailing edge is closed or not. The default is True.

    Returns
    -------
    4 values: x_upper, y_upper, x_lower, y_lower
    '''
    
    # parameters defining an airfoil
    M = float(name[0])/100      # max camber as percentage of chord 
    P = float(name[1])/10       # position of max camber as a tenth of chord
    T = float(name[2:])/100     # max thickness as percent of chord

    # constants 
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    
    if not closed:
        a4 = -0.1015
    else:
        a4 = -0.1036 # for closed trailing edge
    
    #grid
    x = np.linspace(0,1,grid_points)
    
    # Camber and gradient 
    y_camber = np.zeros(grid_points)
    dyc_dx = np.zeros(grid_points)
    theta = np.zeros(grid_points)
    
    # thickness distribution
    yt = (T*5)*(a0*np.sqrt(x) + a1*(x) + a2*(x)**2 + a3*(x)**3 + a4*(x)**4)
    
    if P==0:
        xu = x
        xl = x
        yu = yt
        yl = -yt
    else:
        for i in range(grid_points):
            if x[i]>=0 and x[i]<P:
                y_camber[i] = (M/P**2)*(2*P*x[i]-x[i]**2)                 # lower than P and greater than 0
                dyc_dx[i] = ((2*M)/(P**2))*(P-x[i])
            elif x[i]>=P and x[i]<=1:                                           
                y_camber[i] = (M/(1-P)**2)*(1-2*P+ 2*P*x[i] - x[i]**2)    # greater than or equal to P and lesser than or equal to 1
                dyc_dx[i] = ((2*M)/(1-P)**2)*(P-x[i])                          
        
            theta[i] = math.atan(dyc_dx[i])
    
        # upper surface 
        xu = x - yt*np.sin(theta)
        yu = y_camber + yt*np.cos(theta)
        
        # lower surface
        xl = x + yt*np.sin(theta)
        yl = y_camber - yt*np.cos(theta)

    return xu,yu,xl,yl 
    
def naca5(name,grid_points,closed=True):
    '''
    Parameters
    ----------
    name : string
        name of the airfoil
    grid_points : int
        number of points to consider on the X-axis
    closed : bool, optional
        To define if the trailing edge is closed or not. The default is True.

    Returns
    -------
    4 values: x_upper, y_upper, x_lower, y_lower
    '''    
    
    # constants 
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    
    if not closed:
        a4 = -0.1015
    else:
        a4 = -0.1036 # for closed trailing edge
    
    #Airfoil Parameters
    cld = int(name[0])*(1.5)/10.0            # First digit is the design lift coefficient in tenths when multiplied by 3/2.
    p = 0.5*int(name[1:3])/100.0             # The next two digits, when divided by 2, give the position of the maximum camber (p) in tenths of chord.
    t = int(name[3:])/100.0                  # The final two digits indicate the maximum thickness (t) in percentage of chord.
    
    P = np.array([0.05,0.1,0.15,0.2,0.25])
    M = np.array([0.0580,0.1260,0.2025,0.2900,0.3910])
    K = np.array([361.4,51.64,15.957,6.643,3.230])
    
    m = np.interp(p,P,M)
    k1 = np.interp(p,P,K)
    
    #grid
    x = np.linspace(0,1,grid_points) 
    
    # thickness distribution
    yt = (t*5)*(a0*np.sqrt(x) + a1*(x) + a2*(x)**2 + a3*(x)**3 + a4*(x)**4)
    
    # Camber and gradient 
    y_camber = np.zeros(grid_points)
    dyc_dx = np.zeros(grid_points)
    theta = np.zeros(grid_points)
    
    if p==0:
        xu = x
        xl = x
        yu = yt
        yl = -yt
    else:
        for i in range(grid_points):
            if x[i]>=0 and x[i]<p:
                y_camber[i] = (k1/6)*((x[i])**3 - 3*m*(x[i])**2 + m**2*(3-m)*(x[i]))             # lower than P and greater than 0
                dyc_dx[i] = (k1/6)*(3*(x[i])**2 - 6*m*(x[i]) + m**2*(3-m))
            elif x[i]>=p and x[i]<=1:                                           
                y_camber[i] = ((k1*m**3)/6)*(1-(x[i]))     # greater than or equal to P and lesser than or equal to 1
                dyc_dx[i] = (-k1*m**3)/6                     
    
        y_camber = (cld/0.3) * y_camber
        dyc_dx = (cld/0.3) * dyc_dx
        theta = [math.atan(dyc_dx[i]) for i in range(dyc_dx.shape[0])]
        
        # upper surface 
        xu = x - yt*np.sin(theta)
        yu = y_camber + yt*np.cos(theta)
        
        # lower surface
        xl = x + yt*np.sin(theta)
        yl = y_camber - yt*np.cos(theta)

    return xu,yu,xl,yl
    
def generate_airfoil(number,gridpoints,closed=False,visualize=True):
     '''
     Generates the coordinates for the Airfoil.
     If visualize is TRUE: will draw the Airofoil on a graph and display it
    
        Returns
        -------
        4 values: x_upper, y_upper, x_lower, y_lower
        
        '''
    
     if len(number)==4:
        xu,yu,xl,yl= naca4(nameNACA,gridpoints,closed)
        if visualize:
            plt.figure(figsize=(10,10))
            plt.plot(xl,yl,'k-')
            plt.plot(xu,yu,'r-')
            plt.title("NACA %d"%int(number))
            plt.gca().set_aspect('equal')
            plt.grid(b=True)
            plt.show()
        return xu,yu,xl,yl
     else:
        xu,yu,xl,yl= naca5(nameNACA,gridpoints,closed)
        if visualize:
            plt.figure(figsize=(10,10))
            plt.plot(xl,yl,'k-')
            plt.plot(xu,yu,'r-')
            plt.title("NACA %d"%int(number))
            plt.gca().set_aspect('equal')
            plt.grid(b=True)
            plt.show()
        return xu,yu,xl,yl

if __name__=="__main__":
    nameNACA = '2408'
    grid_points = 500
    generate_airfoil(nameNACA,grid_points,closed=True)
