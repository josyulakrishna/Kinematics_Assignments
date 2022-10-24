#!/usr/bin/env python
# coding: utf-8

# In[1]:
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import symbols
from scipy.integrate import solve_ivp
import numpy as np
import scipy.integrate as it
import pdb


# In[2]:


a1,b1,c1,d1,e1,f1 = symbols('a1 b1 c1 d1 e1 f1') 
a2,b2,c2,d2,e2,f2= symbols('a2 b2 c2 d2 e2 f2') 
t = symbols('t')


# In[3]:


eq_pos_t1 = a1*t**5 + b1*t**4 + c1*t**3 + d1*t**2 + e1*t + f1 
eq_vel_t1 = e1 + 2*d1*t + 3*c1*t**2 + 4*b1*t**3 + 5*a1*t**4
eq_acc_t1 = 2*d1 + 6*c1*t + 12*b1*t**2 + 20*a1*t**3

eq_pos_t2 = a2*t**5 + b2*t**4 + c2*t**3 + d2*t**2 + e2*t + f2 
eq_vel_t2 = e2 + 2*d2*t + 3*c2*t**2 + 4*b2*t**3 + 5*a2*t**4
eq_acc_t2 = 2*d2 + 6*c2*t + 12*b2*t**2 + 20*a2*t**3


# In[4]:


theta1_bvp={'pos':[[0,2],[-np.pi/4., np.pi/4.]], 'vel':[[0,2],[0,0]], 'acc':[[0,2],[0,0]]}


# In[5]:


theta2_bvp={'pos':[[0,2],[0., np.pi/2.]], 'vel':[[0,2],[0,0]], 'acc':[[0,2],[0,0]]}


# In[6]:


def trajGen(theta1_bvp, theta2_bvp):
    expr_list_t1 =  []
    for i in range(len(theta1_bvp['pos'][0])): 
        expr_list_t1.append(eq_pos_t1.subs(t,theta1_bvp['pos'][0][i])-theta1_bvp['pos'][1][i])
    for i in range(len(theta1_bvp['vel'][0])): 
        expr_list_t1.append(eq_vel_t1.subs(t,theta1_bvp['vel'][0][i])-theta1_bvp['vel'][1][i])
    for i in range(len(theta1_bvp['vel'][0])): 
        expr_list_t1.append(eq_acc_t1.subs(t,theta1_bvp['acc'][0][i])-theta1_bvp['acc'][1][i])
    expr_list_t2 =  []
    for i in range(len(theta2_bvp['pos'][0])): 
        expr_list_t2.append(eq_pos_t2.subs(t,theta2_bvp['pos'][0][i])-theta2_bvp['pos'][1][i])
    for i in range(len(theta2_bvp['vel'][0])): 
        expr_list_t2.append(eq_vel_t2.subs(t,theta2_bvp['vel'][0][i])-theta2_bvp['vel'][1][i])
    for i in range(len(theta2_bvp['vel'][0])): 
        expr_list_t2.append(eq_acc_t2.subs(t,theta2_bvp['acc'][0][i])-theta2_bvp['acc'][1][i])
    return solve(expr_list_t1), solve(expr_list_t2)


# In[7]:


theta1co, theta2co = trajGen(theta1_bvp, theta2_bvp) #test coeffes fetch


# In[8]:


def trajectory(time1, theta1co, theta2co): 
    theta1co[t]=time1
    theta2co[t]=time1
    q = [eq_pos_t1.subs(theta1co),eq_pos_t2.subs(theta2co)]
    qdot = [eq_vel_t1.subs(theta1co),eq_vel_t2.subs(theta2co)]
    qddot = [eq_acc_t1.subs(theta1co),eq_acc_t2.subs(theta2co)]
    return q, qdot, qddot


# In[9]:


q, qdot, qddot = trajectory(0, theta1co, theta2co) #test trajecgenrator


# In[10]:


# https://www.ijser.org/researchpaper/Kinematic-and-Dynamic-Analysis-of-Two-Link-Robot-Arm-using-PID-with-Friction-Compensator.pdf
# def senseAngle(x, y, l1, l2 ):
#     pass
# use for inverse kinematics calculation


# In[11]:


# th1, th2, u, v = symbols("th1 th2 u v", cls=Function)
# tau1, tau2 = symbols("tau1 tau2")
# k1,k2,k3,k4,k5,k6,k7 = symbols("k1 k2 k3 k4 k5 k6 k7")


# In[12]:


# ddtheta1_d/dt^2
# eq3 = k1*u(t).diff(t) + k2*v(t).diff(t) - k3*u(t)*v(t) - k4*v(t)**2 - tau1 
# ddtheta2_d/dt^2
# eq4 = k5*u(t).diff(t) + k6*v(t).diff(t) - k7*u(t)*v(t) - tau2


# In[13]:

# theta1[-1], theta2[-1], theta1dot[-1], theta2dot[-1], tau1, tau2
def ddTh_dT(t, Y, th1, th2, th1dot, th2dot, Y3, Y4):
    # pass the current values as the initial values of theta, compute torque from PID
    m2 = 2.
    m1 = 3.
    I1 = 2.
    I2 = 1.
    L1 = 1.
    g = 9.8
    r1 = 0.5
    r2 = 0.5

    u = th1dot  # u = dtheta1_dt
    v = th2dot # v = dtheta2_dt
    th1 = th1  # theta1
    th2 = th2  # theta2

    tau1 = Y3  # Torque1
    tau2 = Y4  # Torque2

    k1 = I1 + m1 * r1 ** 2 + m2 * r2 ** 2 + m2 * L1 ** 2 + 2 * m2 * L1 * r2 * np.cos(th1)
    k2 = m2 * r2 ** 2 + m2 * L1 * r2 * np.cos(th2)
    k3 = 2. * L1 * m2 * r2 * np.sin(th2)
    k4 = 1. * L1 * m2 * r2 * np.sin(th2)
    k51 = g * ((m1 * r1 + m2 * L1) * np.cos(th1) + m2 * r2 * g * np.cos(th1 + th2))

    k5 = m2 * r2 ** 2 + m2 * L1 * r2 * np.cos(th2)
    k6 = m2 * r2 ** 2 + I2
    k7 = 1. * L1 * m2 * r2 * np.sin(th2)
    k8 = m2 * r2 * g * np.cos(th1 + th2)

    du_dt = -(((-k51) * k6 + k2 * k8 + k6 * tau1 - k2 * tau2 - k2 * k7 * u ** 2 + k3 * k6 * u * v + k4 * k6 * v ** 2) / (
                          k2 * k5 - k1 * k6))
    dv_dt = -(((-k5) * k51 + k1 * k8 + k5 * tau1 - k1 * tau2 - k1 * k7 * u ** 2 + k3 * k5 * u * v + k4 * k5 * v ** 2) / (
                          (-k2) * k5 + k1 * k6))

    return [Y[0], Y[1], du_dt, dv_dt]


# derivatives = solve_ivp(ddTh_dT, [0,1],[-np.pi/4.,0],args= (0, 0, 0.1, 0.1))


# In[15]:


# theta1_ = it.cumtrapz(derivatives.y[0],initial=-np.pi/4)[0]


# Pseudocode for PID control
# ![image.png](attachment:image.png)
# computing torques as 
# $\tau = $
# ![image-4.png](attachment:image-4.png)
# ![image-3.png](attachment:image-3.png)

# In[16]:




def pidLoopbck(t,q,qdot,qddot):
#     intialise errors by 1 for steady point control
#     theta1_e = [1]
#     theta1_edot = [0]
#     theta1_eddot = [0]

#     theta2_e =[1]
#     theta2_edot = [0]
#     theta2_eddot = [0]

    #list for storing errors n theta1
    l_theta1_e=[1]
    l_theta1_edot=[0]
    l_theta1_eddot=[0]
    #list for storing theta1
    l_theta1_=[np.float(q[0,0])]
    l_theta1_dot=[np.float(qdot[0,0])]
    l_theta1_ddot=[np.float(qddot[0,0])]

    #list for storing errors in theta2
    l_theta2_e=[1]
    l_theta2_edot=[0]
    l_theta2_eddot=[0]
    #list for storing theta2
    l_theta2_=[np.float(q[0,1])]
    l_theta2_dot=[np.float(qdot[0,1])]
    l_theta2_ddot=[np.float(qddot[0,1])]


    Kp1 = 2.5
    Kd1 = 2.3
    Ki1 = 0.3


    Kp2 = 2.4
    Kd2 = 2.7
    Ki2 = 0.3

    dt = 1./len(t)

    for i in range(len(t)):
        print(t[i])
        theta1_e = np.float(q[i,0])-l_theta1_[i]
        theta1_e_dot = np.float(qdot[i,0]) - l_theta1_dot[i]
        theta1_e_ddot = np.float(qddot[i,0])-l_theta1_ddot[i]

#         pdb.set_trace()

        l_theta1_e.append(theta1_e)
        int_theta1_e = it.simpson(l_theta1_e[:i+1])

        tau1 = Kp1*theta1_e+Kd1*theta1_e_dot+Ki1*int_theta1_e


        theta2_e = np.float(q[i,1])-l_theta2_e[i]
        theta2_e_dot = np.float(qdot[i,1])-l_theta2_dot[i]
        theta2_e_ddot =np.float( qddot[i,1])-l_theta2_ddot[i]

        l_theta2_e.append(theta2_e)
        int_theta2_e = it.simpson(l_theta2_e[:i+1])

        tau2 = Kp2*theta2_e+Kd2*theta2_e_dot+Ki2*int_theta2_e

        #function, time, initial states theta1dot and theta2dot
        if i==0:
            theta_dots = solve_ivp(ddTh_dT, [t[i],t[i]+dt], [0.,0.],args= (np.float(q[i,0]), np.float(q[i,1]), tau1, tau2))
        else:
            theta1dot = l_theta1_dot[i-1]-l_theta1_dot[i]
            theta2dot = l_theta2_dot[i-1]-l_theta2_dot[i]
            theta_dots = solve_ivp(ddTh_dT, [t[i],t[i]+dt], [theta1dot,theta2dot],args= (np.float(q[i,0]), np.float(q[i,1]), tau1, tau2))
        theta1 = it.simpson(theta_dots.y[0]) #, initial=l_theta1_[i])
        theta2 = it.simpson(theta_dots.y[1]) #, initial=l_theta1_[i])
        l_theta1_.append(theta1)
        l_theta2_.append(theta1)
        #caluculate theta dot actual
        l_theta1_dot.append(theta_dots.y[0][-1])
        l_theta2_dot.append(theta_dots.y[1][-1])
        #caluclate theta ddot actual
        l_theta1_ddot.append(theta_dots.y[0][-1]-theta_dots.y[0][0])
        l_theta2_ddot.append(theta_dots.y[1][-1]-theta_dots.y[1][0])
    plt.plot(l_theta1_e, "-r")
    plt.plot(l_theta1_, "-b")
    plt.show()


def pidLoop(t, q, qdot, qddot):

    theta1 = [float(q[0,0])]
    theta2 = [float(q[0,1])]
    theta1dot = [float(qdot[0,0])]
    theta2dot = [float(qdot[0,1])]
    theta1ddot = [float(qddot[0,0])]
    theta2ddot = [float(qddot[0,1])]

    #list for storing errors in theta2
    theta1_e = [1];
    theta1_edot = [0];
    theta1_eddot = [0];

    theta2_e = [1];
    theta2_edot = [0];
    theta2_eddot = [0];

    m2 = 2.
    m1 = 3.
    I1 = 2.
    I2 = 1.
    L1 = 1.
    g = 9.8
    r1 = 0.5
    r2 = 0.5

    Kp1 = 0.3
    Kd1 = 0.3
    Ki1 = 3

    Kp2 = 3
    Kd2 = -0.3
    Ki2 = 10
    dt = 1. / len(t)
    tau1_l = []
    tau2_l = []

    for i in range(len(t)):
        print(t[i])
        theta1d = q[i,0]
        theta1dotd = qdot[i,0]
        theta1ddotd = qddot[i,0]

        theta2d = q[i,1]
        theta2dotd = qdot[i,1]
        theta2ddotd = qddot[i,1]

        # calucluate error;

        theta1_e.append(theta1d - theta1_e[-1])
        theta1_edot.append(theta1dotd - theta1_edot[-1])
        theta1_eddot.append(theta1ddotd - theta1_eddot[-1])

        theta2_e.append(theta2d - theta2_e[-1])
        theta2_edot.append(theta2dotd - theta2_edot[-1])
        theta2_eddot.append( theta2ddotd - theta2_eddot[-1])

        u1 = Kp1 * theta1_e[-1] + Kd1 * theta1_edot[-1] + Ki1 * it.simpson(theta1_e)


        u2 = Kp2 * theta2_e[-1] + Kd2 * theta2_edot[-1] + Ki2 * it.simpson(theta2_e)

        inp_pid = np.array([[u1], [u2]])

        td = np.array([[theta1dot[-1]], [theta2dot[-1]]])

        tdd = np.array([[theta1ddot[-1]], [theta2ddot[-1]]])

        M = np.array([[I1 + m1 * r1 ** 2 + m2 * r2 ** 2 + m2 * L1 ** 2 + 2 * m2 * L1 * r2 * np.cos(theta2[-1]),
                       m2 * r2 ** 2 + m2 * L1 * r2 * np.cos(theta2[-1])],
                      [m2 * r2 ** 2 + m2 * L1 * r2 * np.cos(theta2[-1]), m2 * r2 ** 2 + I2]])
        C = np.array([[-m2 * L1 * r2 * theta2dot[-1] * np.sin(theta2[-1]),
                       -m2 * L1 * r2 * theta1dot[-1] * np.sin(theta2[-1]) - m2 * L1 * r2 * theta2dot[-1] * np.sin(
                           theta2[-1])], [m2 * L1 * r2 *theta1dot[-1] * np.sin(theta1[-1]), 0]])

        G = np.array(
            [[g * ((m1 * r1 + m2 * L1) * np.cos(theta1[-1]) + m2 * r2 * g * np.cos(theta2[-1] + theta2[-1]))],
             [ m2 * r2 * g * np.cos(theta2[-1] + theta1[-1])]])

        result1 = M@(tdd+ inp_pid) + C@td + G
        tau1 = result1[0,0]
        tau2 = result1[1,0]
        theta_dots = solve_ivp(ddTh_dT,
                               [t[i], t[i] + dt],
                               [theta1[-1], theta2[-1], theta1dot[-1], theta2dot[-1]],
                               args=(theta1[-1], theta2[-1], theta1dot[-1], theta2dot[-1], tau1, tau2)
                               )
        theta1.append(theta_dots.y[0,-1])
        theta2.append(theta_dots.y[1,-1])
        theta1dot.append(theta_dots.y[2,-1]);
        theta2dot.append(theta_dots.y[3, -1]);

        theta1ddot.append(theta1dot[-1] - theta1dot[-2])
        theta2ddot.append(theta2dot[-1] - theta2dot[-2])

        tau1_l.append(tau1)
        tau2_l.append(tau2)

    plt.plot(q[:,0], "-r", label="reference")
    plt.plot(theta1, "-b", label="actual")
    plt.legend(loc="upper left")
    plt.xlabel("time steps")
    plt.ylabel("theta")
    plt.show()

    plt.plot(q[:,1], "-r", label="error")
    plt.plot(theta2, "-b", label="actual")
    plt.legend(loc="upper left")
    plt.xlabel("time steps")
    plt.ylabel("theta value")
    plt.show()


    plt.plot(tau1_l, "-r", label="tau1")
    plt.plot(tau2_l, "-b", label="tau2")
    plt.legend(loc="upper left")
    plt.xlabel("time steps")
    plt.ylabel("tau value")
    plt.show()

# In[ ]:
q = [] #{theta1,theta2}
qdot = [] #{theta1dot,theta2dot}
qddot = [] #{theta1ddot,theta2ddot}
timesteps = 50
for i in np.linspace(0,2,timesteps):
    q1, qdot1, qddot1 = trajectory(i, theta1co, theta2co)
    q.append(q1)
    qdot.append(qdot1)
    qddot.append(qddot1)


# In[17]:


q=np.array(q)
qdot=np.array(qdot)
qddot=np.array(qddot)

pidLoop(np.linspace(0,2,timesteps),q,qdot,qddot)# integrate.simpson([0])

# # sns.set()
# plt.plot(l_theta1_e, "-r", label="error")
# plt.plot(l_theta1_, "-b", label="actual")
# # plt.legend('theta1 actual, theta1 error', ncol=2, loc='upper left');
# plt.legend(loc="upper left")
# plt.xlabel("time steps")
# plt.ylabel("theta")
# plt.ylim(-100, 100)
# plt.show()
#
# plt.plot(l_theta2_e, "-r", label="error")
# plt.plot(l_theta2_, "-b", label="actual")
# # plt.legend('theta1 actual, theta1 error', ncol=2, loc='upper left');
# plt.legend(loc="upper left")
# plt.xlabel("time steps")
# plt.ylabel("theta value")
# plt.ylim(-100, 100)
# plt.show()
# In[ ]:



qdot


# In[ ]:




