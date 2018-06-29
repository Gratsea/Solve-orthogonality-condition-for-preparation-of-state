# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:50:16 2018

@author: kgratsea

Preparation of a desired state over sites {2,n} with n-number of steps. Total final state over sites {1,n+1}. 
Initial state arbitrary superposition over 1 site.
"""

import numpy as np
from scipy.optimize import fsolve


#desired sate over sites {2,n}
#n=3
#target= np.array([[0.6148,0.3492],[0.3493,-0.6148]])

n=5
target=np.array ([[0.52619533, 0.27592482, 0.26093188, 0.2510644 ], [0.2290105 , 0.48391584, 0.3214155 , 0.35487592]]) 

def kate(i,given_state):
    u_i = np.array([[given_state[0][i-1]],[given_state[1][i]]])
    return u_i


def system(d):
    '''
    n=3
    final_state = np.array([[d[0],0.6148,0.3493,0],[0,0.3492,-0.6148,d[1]]])
    '''
    n=5
    final_state = np.array([[d[0],0.52619533, 0.27592482, 0.26093188, 0.2510644,0 ], [0,0.2290105 , 0.48391584, 0.3214155 , 0.35487592,d[1]]]) 

    #system of equations
    equations = []
    #sites 2,3
    #equation1 = np.conj(d[0])* final_state[0][1]+ np.conj(final_state[1][1])*final_state[1][2] + np.conj(final_state[0][1])*final_state[0][2]+np.conj(final_state[1][2])*d[1] 
    
    #condition : the summation is only over the "internal" sites
    s=n-1
    summ = 0
    for i in range (2,s-1,+1) :   
        v_i = kate(i,final_state)
        v_j = kate(n-s+i,final_state)
        summ += np.dot(np.matrix.getH(v_i),v_j)
    
    #the "hole" orthonomality condition, summation over all sites
    summation = summ + np.conj(final_state[1][1])*final_state[1][2]  + np.conj(final_state[0][s-1])*final_state[0][n-1]
    equation1= summation -( np.conj(final_state[0][0])* final_state[0][1]+np.conj(final_state[1][s])*final_state[1][n] )

    #equation1= summ -( np.conj(final_state[0][0])* final_state[0][1]+ np.conj(final_state[1][1])*final_state[1][2]  + np.conj(final_state[0][s-1])*final_state[0][n-1]+np.conj(final_state[1][s])*final_state[1][n] )
    equations.append( equation1 )
    equation2 = np.linalg.norm(final_state)-1
    equations.append( equation2 )
    
    print (d,equation1)
    
    return equations

n=3 #number of sites at the phi state
k=n+1
d0 = fsolve(system,[0.7,0.5]) #real and imag parts of each constant di
print ("d0",d0)
'''
#n=3
final_state1 = np.array([[d0[0]],[0],[0.6148],[0.3492],[0.3493],[-0.6148],[0],[d0[1]]])

final_state2 = np.array([[d0[0],0.6148,0.3493,0],[0,0.3492,-0.6148,d0[1]]])
'''


#n=5
final_state1 =  np.array([[d0[0],0.52619533, 0.27592482, 0.26093188, 0.2510644,0 ], [0,0.2290105 , 0.48391584, 0.3214155 , 0.35487592,d0[1]]]) 


print (final_state1)
