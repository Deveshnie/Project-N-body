# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 23:06:35 2016

@author: Deveshnie
"""

import Simulator
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #check that you're not entering a positive no of particlees, not words etc.
    while True:
        try:
            part = int(input ("Enter The Number Of Particles: \nNote: \n1  --> Single Particle Simulation\n2  --> Double Particle Simulation\n1000 --> N-Body Simulation\n"))
            break
        except:
            print "Enter a positive number"
    
    
    while (part<=0):
        print"Enter a number above 0"
        part = int(input ("Enter The Number Of Particles: \nNote: \n1  --> Single Particle Simulation\n2  --> Double Particle Simulation\n1000 --> N-Body Simulation\n"))
    
    #Entering either 1 or 2 for periodic or non periodic
    while True:
        try:
            periodic=int(input("Enter 1 for Periodic \nEnter 2 For Non-Periodic\n")) 
            break
        except:
            print "Enter 1 or 2"      
   
    while (periodic!=1 and periodic!=2):
        periodic=int(input("Enter 1 for Periodic \nEnter 2 For Non-Periodic\n")) 
                
    #choosing how many blocks for a grid column
    if ((part==1)):
        size=4
    elif(  (part==2)):
        size=4
    else:
        size=8
    
    fun = Simulator.Simulator( part,size,periodic)
    size=fun.size  
    
    #loop for simulations. It runs over the number of simulations chosen.Enter a positive number. not words etc.
    while True:
        try:
            sim = int(input ("Enter The Number Of Runs of the Simulations: "))
            break
        except:
            print "Enter a positive number"      
    
    while (sim<=0):
        print "Enter a positive number"
        sim = int(input ("Enter The Number Of Runs of the Simulations: "))
        
    fun.populate()
    for r in range(sim):
        x,y = [],[]
        
        plt.ion()
        fig = plt.figure()
        print
        print "Simulation: ", r + 1 
        fun.populateDensity()   
        fun.populatePotential() #calls potential
        fun.update()        #update positions, potentials
        for i in range(len(fun.matrix)):
            x.append(fun.matrix[i][1])  
            y.append(fun.matrix[i][2])
        plt.cla()           #clears prev point
        plt.plot(x, y, "*") #each particle is a star
        plt.xlim(xmax=size,xmin=0)  #bounding the grid for x positions
        plt.ylim(ymax=size,ymin=0)        
        plt.show()
        plt.draw()          #draw stars
