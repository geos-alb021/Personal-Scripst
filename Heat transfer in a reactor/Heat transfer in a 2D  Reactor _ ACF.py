###-------------------Heat transfer in a 2D reactor-------------------###

###Authorship: Alberto F. Coto Fonseca

#Description: This script determines the process of heat transfer if a 2D reactor using implicit and explicit approximations. 
#The implicit approximation is used for the calculation of the external nodes of the walls that are exposed to the air and in the reactor's base (isolated side), 
#while the explicit approximatios are implemented for the inner nodes of the reactor.
#The inner temperature of the reactor is 60°C and the environemtal temperature is 20°C. The process ends until the temperature
#of the reactor reaches 50°C. The heat transfer occurs in the left, upper and right wall of the reactor, due to the base been isolated.

#Importing libraries
import numpy as np
import matplotlib.pyplot as plt

#Characteristics of the 2D reactor
X=1.5                      #Width, meters
Y=2						   #Height, meters
alfa=22.18*(1/(1000*1000))   #m2/s
k = 50                     # W/m*k Steel's heat conductivity 
h = 5                      # W/m^2*K Heat transfer coefficient

dx=0.10   #distance between nodes in x axis, meters 
dy=dx     #distance between nodes in y axis, meters                    
	  
x_vec=int(X/dx) #Number of nodes in X axis
y_vec=int(Y/dy) #Number of nodes in Y axis

N=2*(x_vec-2)+2*y_vec #External total nodes (used for determining the total of nodes in the border of the reactor)

#Timesteps conditions
dt=5    #timesteps, seconds

#Initial temperature conditions
Tamb=20     #Environment temperature, °C
Temp = np.zeros((y_vec,x_vec))+60   #Initial temperature vector
Temp2 = np.zeros((2*(x_vec-2)+2*y_vec,1)) #External temperature calculated values vector
K=np.zeros((2*(x_vec-2)+2*y_vec,1))       #K's value vector

#Calculation of the necessary parameters for the implicit and explicit approximations
A=2*h*(dx/k)*Tamb #Corners
B=2*(h*(dx/k)+2) #Walls and Laterals
C = 2*(h*(dx/k)+1) #Corners
F = alfa*(dt/dx**2) #Fourier number for the inner temperature nodes

#Creating the matrix based on the dimensions of the temperature vector
Va = np.zeros((1,(2*(x_vec-2)+2*y_vec)-1))+1 
Vb = np.zeros((1,2*(x_vec-2)+2*y_vec))

#Filling the square matrix for the calculation of the temperatures
for n in range(N): 
  if n==0:  #Reactor's Left Lower corner
    Vb[0,0]=-C
  else:
      if n==y_vec:  #Reactor's Left Upper corner
       Vb[0,y_vec-1]=-C  
      else:
          if n==y_vec+x_vec:   #Reactor's Right Upper corner
           Vb[0,y_vec+x_vec-2]=-C
          else:
              if n==N-x_vec+2:   #Reactor's Right Lower corner
               Vb[0,N-x_vec+1]=-C
              else:
                  if n>=1 and n<y_vec-1:   #Reactor's left wall 
                   Vb[0,1:y_vec-1]=-B
                  else:
                      if n>=y_vec and n<y_vec+x_vec:   #Reactor's upper wall
                       Vb[0,y_vec:y_vec+x_vec-2]=-B
                      else:
                          if n>y_vec+x_vec-1 and n<N-x_vec+2:   #Reactor's righ wall
                           Vb[0,y_vec+x_vec-1:N-x_vec+2]=-B
                          else:
                              n>=N-x_vec+2 and n<=N   #Reactor's lower wall (isolated) 
                              Vb[0,N-x_vec+2:N]=-4
                              #print(n)
                              
#Creating the matrices with individual diagonals
ma = np.zeros((66,66))
MA = np.fill_diagonal(ma,Vb)

mb1 = np.zeros((66,66))
i,j = np.indices(mb1.shape)
mb1[i==j-1]=1

mb2 = np.zeros((66,66))
i2, j2 = np.indices(mb2.shape)
mb2[i2-1==j2]=1

#Creating the complete matrix
M = ma+mb1+mb2
M[0,N-1]=1
M[N-1,0]=1

### Determining the external nodes temperature through implicit approximation ###
Tmean=100
tdx=1

while Tmean>50:  #The process is repeated until the mean temperature reaches 50°C
    for idx in range(N-1):  #Change along the border of the reactor
        if idx==0:
            K[0,0] = -A  #Reactor's Left Lower corner
        else: 
            if idx==y_vec:
                K[y_vec-1,0]=-A  #Reactor's Left Upper corner
            else:
                if idx==y_vec+x_vec:
                    K[y_vec+x_vec-2,0]=-A   #Reactor's Right Upper corner
                else:
                    if idx==N-x_vec+2:
                        K[N-x_vec+1,0]=-A  #Reactor's Right Lower corner
                    else:
                        if idx>=1 and idx<y_vec-2:
                            K[1:y_vec-1,0] = 2*Temp[idx,0]  #Reactor's left wall 
                        else:
                            if idx>=y_vec-1 and idx<=y_vec+x_vec-2:
                                 K[y_vec:y_vec+x_vec-2,0]=2*Temp[0,idx-(y_vec-1)] #Reactor's upper wall
                            else:
                                if idx>=y_vec+x_vec-1 and idx<(N-x_vec):
                                    K[y_vec+x_vec-1:N-x_vec+1,0]=2*Temp[idx-(N//2),x_vec-1]  #Reactor's right wall 
                                else:
                                    if idx>N-x_vec+2 and idx<=N:   #Reactor's lower wall (isolated)  
                                         K[N-x_vec+2:N,0]=2*Temp[y_vec-1,(N+2)-idx]
                                                                                                                      
#Creating the vector that contains the external temperature on time n 
    InvM = np.linalg.inv(M)
    Temp2 = np.dot(InvM,K)*-1

#Creating the matrix wthat contains the external nodes temperatures
    Text = np.zeros((y_vec,x_vec))+Temp #External nodes matrix (20x15)
    Text[1:y_vec-1, 0]=np.flip(Temp2[1:y_vec-1,0])  #Left wall 
    Text[0, 0:x_vec-1]=Temp2[y_vec-1:y_vec+x_vec-2,0]  #Upper wall 
    Text[0:y_vec-1, x_vec-1]=Temp2[y_vec+x_vec-2:N-x_vec+1,0] #Right wall
    Text[y_vec-1, 0:x_vec-1]=Temp[y_vec-1,0:x_vec-1]  #Lower wall

    Temp = Text  #Overwrites the temperature original matrix with the new external temperatures 
    
                                       
### Determining the internal nodes with explicit aproximation ###

    for idX in range(1,y_vec-2):  #Variation in rows
      for idY in range (1,x_vec-1):  #Variation in columns
                                                 
            Temp[idX,idY]= F*(Temp[idX+1,idY]+Temp[idX-1,idY]+Temp[idX,idY+1]+Temp[idX,idY-1])+(1-4*F)*Temp[idX,idY]
                                             
    Tmean=(np.mean(Temp))
    tdx = tdx+1
    print("Time step:", tdx) #Prints the timesteps until the process is concluded.
    
#Plotting the reactor's temperature condition at 50°C
plt.imshow(Temp, cmap='jet')
plt.xlabel("Width (m)")
plt.ylabel("Height (m)")
plt.colorbar(label="Temperature (°C)")
plt.show()   
                       















