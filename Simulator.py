# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 23:05:08 2016

@author: Deveshnie
"""

import numpy as np
from numpy.fft import fft,ifft
import random 

class Simulator:
    def __init__(self, partNo,size,periodic):
        self.dt = 0.05  
        self.partNo = partNo  #Particle Count
        self.size = size #no of blockks in 1 column or row of grid
        self.matrix = [[0 for i in xrange(7)] for i in xrange(self.partNo)] #Particle Property Array
        self.density = np.zeros([self.size, self.size])  #Density matrix
        self.potential = np.zeros([self.size, self.size]) #Potential matrix
        self.softPotential = np.zeros([self.size, self.size]) #Softened Potential matrix
        self.periodic = periodic #1 is periodic, 2 is non periodic
        """
        0 = mass
        1 = position x
        2 = position y
        3 = force x
        4 = force y
        5 = velocity x
        6 = velocity y
        """

    def populate(self):
        if(self.partNo == 1):
            self.matrix [0][0] = 1
            self.matrix [0][1] = self.size/2.0
            self.matrix [0][2] = self.size/2.0
            self.matrix [0][3] = 0
            self.matrix [0][4] = 0
            self.matrix [0][5] = 0
            self.matrix [0][6] = 0

        if(self.partNo == 2):
            self.matrix [0][0] = 1
            self.matrix [0][1] = self.size * 0.67   #x position
            self.matrix [0][2] = self.size * 0.67   #y position
            self.matrix [0][3] = 0
            self.matrix [0][4] = 0
            self.matrix [0][5] = 0
            self.matrix [0][6] = 0
            self.matrix [1][0] = 1
            self.matrix [1][1] = self.size * 0.33
            self.matrix [1][2] = self.size * 0.33
            self.matrix [1][3] = 0
            self.matrix [1][4] = 0
            self.matrix [1][5] = 0
            self.matrix [1][6] = 0

        if(self.partNo >2 ):
            for i in range(0, self.partNo):
                self.matrix [i][0] = 1
                self.matrix [i][1] = random.uniform(0, self.size - 1)   #generate random numbers
                self.matrix [i][2] = random.uniform(0, self.size - 1)
                self.matrix [i][3] = 0
                self.matrix [i][4] = 0
                self.matrix [i][5] = 0
                self.matrix [i][6] = 0

    def populateDensity(self):
        for i in range(len(self.matrix)):   
            x = int(self.matrix[i][1])  #converting number to in for x
            y = int(self.matrix[i][2])

            if (x > self.size - 1):
                x = self.size - 2
            if (y >  self.size - 1):
                y = self.size - 2

            self.density[x][y] += 1 #making density matrix by counting no of particles in each cell
            
    def populatePotential(self):
        for n in range(0,self.size/2):  #making softened potential matrix on left half of grid
            for m in range(0,self.size/2):
                if(n == 0 and m == 0):
                    self.softPotential[0][0] = 1
                    self.softPotential[0][self.size - 1] = 1
                else:
                    self.softPotential[n][m] = 1. / np.sqrt(n**2 + m**2)
                    mm=-1*(m+1)
                    self.softPotential[n][self.size + mm] = 1. /  np.sqrt(n**2 + mm**2)
        count=-1*(self.size/2)-1    #making softened potential matrix on right half of grid
        while(count<-1):
            count+=1
            for n in range(0,self.size/2):   
                self.softPotential[count][n] = 1. / np.sqrt(count**2 + n**2)
                nn=-1*(n+1)
                self.softPotential[count][self.size+nn]=1./np.sqrt(nn**2 + count**2)
                
        fft1 = fft(self.density)
        fft2 = fft(self.softPotential)
        self.potential = np.real(ifft(fft1 * fft2))    #potential matrix by convolving density and softened potential matrix

    def update(self):

        i = 0
        while (i < len(self.matrix) - 1):
            deleted = False
            j = i
            for i in range(j, len(self.matrix)):
                n = len(self.matrix)
                x = self.matrix[i][1]   #using the x component of a particle
                y = self.matrix[i][2]   #using the y component of a particle
                
                #here it for non periodic, the potentials on the left, right, above and below the cell are determined.
                #this is done for 4 corners and 4 sides of the grid. It is also done for the middle of the grid.
                #If there is no block to the left, right above or below a cell, the potential on the block excludes 
                #those blocks(which are not there)
                if (self.periodic == 2):    #Non-periodic
                    #cell on left upper corner of grid
                    if ((x - 1 < 0) and (y - 1 < 0)):
                        
                        #left and above the cell- potentials dont exist
                        l = 0
                        r = self.potential[x + 1][y]    #potential of the block on the right exists
                        u = 0
                        d = self.potential[x][y + 1]    #potential below this block exists
                    
                    #cell on left lower corrner of grid
                    elif ((x - 1 < 0) and (y + 1 > self.size)):
                        #cells to the left and below the cell dont exist
                        l = 0
                        r = self.potential[x + 1][y]    #the potential to the right has x+1 because x increases to the right
                        u = self.potential[x][y - 1]    #the potential above uses y-1 because (0,0 ) 
                                                        #is on left upper corner so y increases downward
                        d = 0

                    #cell on right uppper corner
                    elif ((x + 1 > self.size) and (y - 1 < 0)):
                        #right and up are 0 since they dont exist
                        l = self.potential[x - 1][y]    #left potential is at x-1
                        r = 0
                        u = 0
                        d = self.potential[x][y + 1]    #potential cell below this cell exists at y+1
                        
                    #cell on right lower corner of grid
                    elif ((x + 1 > self.size) and ((y + 1 > self.size))):

                        l = self.potential[x - 1][y]
                        r = 0
                        u = self.potential[x][y - 1]
                        d = 0
                    
                    #uppper boundary of grid
                    elif ((y - 1 < 0) and (0 <= x <= self.size)):

                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = 0
                        d = self.potential[x][y + 1]
                        
                    #lower boundary of grid
                    elif ((y + 1 > self.size) and (0 <= x <= self.size)):

                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = 0
                    
                    #left boundary
                    elif ((0 <= y <= self.size) and ((x - 1 < 0))):

                        l = 0
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]

                    #right boundary
                    elif ((0 <= y <= self.size) and ((x + 1 > self.size))):

                        l = self.potential[x - 1][y]
                        r = 0
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]

                    else:
                        #not on boundary
                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]

                #the same thing is done for peiodic boundary conditions.
                #In this case, if the block left, right, above or below a cell does not exist(corners or boundaries)
                #Instead of making the potentials 0, we use the corresponding cells on the opposite boundaries of the grid
                if (self.periodic == 1):    #Periodic
                    #left upper corner
                    if ((x - 1 < 0) and (y - 1 < 0)):

                        l = self.potential[self.size - 1][y]    #the potential cell on the left side of the grid doesnt exist
                                                                #so the potential of the upper right corner of the grid is taken instead
                        r = self.potential[x + 1][y]            #potential to the right exists so it is at cell x+1
                        u = self.potential[x][self.size - 1]    #potential above the grid cell doesnt exist so the potential 
                                                                #on the left bottom corner is used instead
                        d = self.potential[x][y + 1]            #potential below the cell exists. Since y increases downward, we take
                                                                #potential y+1
                    #left lower corner
                    elif ((x - 1 < 0) and (y + 1 > self.size)):

                        l = self.potential[self.size - 1][y]    #potential to left doesnt exist. Potential on right lower corner used
                        r = self.potential[x + 1][y]            #potential to right exists at x+1
                        u = self.potential[x][y - 1]            #potential above exists at y+1
                        d = self.potential[x][0]                #potential below cell doesnt exist so potential on left upper corner used
                                                                #instead
                    #right upper corner
                    elif ((x + 1 > self.size) and (y - 1 < 0)):

                        l = self.potential[x - 1][y]            #potential to left exists at x-1
                        r = self.potential[0][y]                #the cell on the right doesnt exist so we use the left upper corner
                        u = self.potential[x][self.size - 1]    #potential cell above this cell doesnt exist. Choose the cell on the bottom right 
                        d = self.potential[x][y + 1]            #cell below this one exists at y+1
                    #right lower corner
                    elif ((x + 1 > self.size) and ((y + 1 > self.size))):

                        l = self.potential[x - 1][y]
                        r = self.potential[0][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][0]
                        
                    #potential on upper grid boundary 
                    elif ((y - 1 < 0) and (0 <= x <= self.size)):

                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][self.size - 1]
                        d = self.potential[x][y + 1]

                    #lower boundary
                    elif ((y + 1 > self.size) and (0 <= x <= self.size)):

                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][0]
                    #left boundary
                    elif ((0 <= y <= self.size) and ((x - 1 < 0))):

                        l = self.potential[self.size - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]

                    #right boundary
                    elif ((0 <= y <= self.size) and ((x + 1 > self.size))):

                        l = self.potential[x - 1][y]
                        r = self.potential[0][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]

                    #not on boundaries: in middle of grid
                    else:

                        l = self.potential[x - 1][y]
                        r = self.potential[x + 1][y]
                        u = self.potential[x][y - 1]
                        d = self.potential[x][y + 1]
                
                #Calc Force
                self.matrix[i][3] = -1 * (l - r) / 2.  #force x compoent potential on left - potential on right of cell
                self.matrix[i][4] = -1 * (u - d) / 2.  #force y component potential above-potential below cell

                #Calc Velocity
                self.matrix[i][5] += self.matrix[i][3] * self.dt  #velocity x component is updated using force x component
                self.matrix[i][6] += self.matrix[i][4] * self.dt  #velocity y component is updated using force y component
                #force is used to calculate velocity because mass=1 so acceleration=force in this case
                
                #Calc Position
                self.matrix[i][1] += self.matrix[i][5] * self.dt  #velocity x component is used to update x component of particle
                self.matrix[i][2] += self.matrix[i][6] * self.dt

                if (self.periodic == 1):    #Periodic
                    #these conditions ensure that if the particle exits the grid on one boundary,
                    #that it will reenter the grid on another boundary 
                    #this is done for x and y componenets
                    while (self.matrix[i][1] <= 0): #x component, on or beyond left boundary
                        self.matrix[i][1] = self.size + self.matrix[i][1]       #move x component to right side of grid
                    while (self.matrix[i][1] >= self.size): #x, right boundary
                        self.matrix[i][1] = self.matrix[i][1] - self.size   #move particle to left side
                    while (self.matrix[i][2] <= 0): #y, on/above upper boundary
                        self.matrix[i][2] = self.size + self.matrix[i][2]   #place y near lower boundary
                    while (self.matrix[i][2] >= self.size): #on/below lower boundary
                        self.matrix[i][2] = self.matrix[i][2] - self.size   #place near upper boundary

                if (self.periodic == 2):    #Non periodic
                    #this ensures that if a particle exits the grid on one boundary or sits on the boundary, that it
                    #is deleted 
                    #this is done for x and y componenets of the particle
                    while (self.matrix[i][1] <= 0):
                        del self.matrix[i]
                        n = n - 1
                        deleted = True
                        i-=1
                       
                    while (self.matrix[i][1] >= self.size):
                        del self.matrix[i]
                        n = n - 1
                        deleted = True
                        i-=1
                        
                        
                    while (self.matrix[i][2] <= 0):
                        del self.matrix[i]
                        n = n - 1
                        deleted = True
                        i-=1
                     
                        
                    while (self.matrix[i][2] >= self.size):
                        del self.matrix[i]
                        n = n - 1
                        deleted = True
                        i-=1
                        
                        
                    if deleted:
                        break