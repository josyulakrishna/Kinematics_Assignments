def PID(t, q, dq, ddq)

    I1=1; I2=1; r1=.5;
    r2=.5; L1=1; L2=1;
    m1=3; m2=2;

    tt = 0:0.1:2;
    theta1 = [q(1,j)];
    theta2 = [q(2,j)];

    tau1_l = [];
    tau2_l = [];

    theta1dot = [dq[1,j]];
    theta2dot = [dq(2,j)];

    theta1ddot = [ddq(1,j)];
    theta2ddot = [ddq(2,j)];

    % errors;
    theta1_e = [1];
    theta1_edot = [0];
    theta1_eddot = [0];

    theta2_e = [1];
    theta2_edot = [0];
    theta2_eddot = [0];

    N=20;
    tf=2;
    timegap = tf / (N - 1);
    for i=1:N
    %intialise variables;
        % actual values

        %desired values ;
        theta1d = q(1,i);
        theta1dotd = dq(1,i);
        theta1ddotd = ddq(1,i);

        theta2d = q(2,i);
        theta2dotd = dq(2,i);
        theta2ddotd = ddq(2,i);
        %calucluate error;

        theta1_e(end+1) = theta1d - theta1_e(end);
        theta1_edot(end+1) = theta1dotd - theta1_edot(end);
        theta1_eddot(end+1) = theta1ddotd - theta1_eddot(end);

        theta2_e(end+1) = theta2d - theta2_e(end);
        theta2_edot(end+1) = theta2dotd - theta2_edot(end);
        theta2_eddot(end+1) = theta2ddotd - theta2_eddot(end);

        Kp1 = -0.2;
        kd1 = 0.0 ;
        ki1 = 0.0 ;

        Kp2 = -0.2;
        kd2 = 0.0 ;
        ki2 = 0.0 ;

        u1 = theta1ddotd + Kp1*(theta1_e(end))+kd1*(theta1_edot(end))+ki1*trapz(theta1_e);

        u2 = theta2ddotd + Kp2*(theta2_e(end))+kd2*(theta2_edot(end))+ki2*trapz(theta2_e);

        u = theta1dot(end);
        v = theta2dot(end);

        g = 9.8;

        k1 = I1+I2+L1^2+m2+m1*r1^2+m2*r2^2+2*L1*m2*r2*cos(theta2(end));
        k2 = I2+m2*r2^2+L1*m2*r2*cos(theta2(end));
        k3 = -L1*m2*r2*sin(theta2(end))*u^2;
        k4 =2*L1*m2*r2*sin(theta2(end))*u*v;
        k5 = -g*m2*r2*sin(theta1(end)+theta2(end))-L1*g*m2*sin(theta1(end))-g*m1*r1*sin(theta1(end));
        k6 =I2+m2*r2^2+L1*m2*r2*cos(theta2(end));
        k7 =I2+m2*r2^2 ;
        k8 = L1*m2*r2*sin(theta2(end))*u^2;
        k9=-g*m2*r2*sin(theta1(end)+theta2(end));

        M=[k1 k2;k6 k7];
        C=[k3 k4;k8 0];
        G=[k5;k9];


        re1  = M*([theta1ddotd;theta2ddotd]+[u1;u2])+ C + G;
        tau1 = re1(1);
        tau2 = re1(2);
        para=@(t,Y)dtheta_dt(t, Y, theta1(end), theta2(end), theta1dot(end), theta2dot(end), tau1, tau2 );
        [t,y]=ode45(para, [timegap*(i - 1), timegap*(i - 1)+0.1], [theta1(end); theta2(end); theta1dot(end); theta2dot(end)]);
    %     tempr = trapz(y);
        theta1(end+1) = y(end,1);
        theta2(end+1) = y(end,2);
        theta1dot(end+1) = y(end,3);
        theta2dot(end+1) = y(end,4);
        theta1ddot(end+1) = theta1dot(end-1) - theta1dot(end);
        theta2ddot(end+1) = theta2dot(end-1) - theta2dot(end);
        j=j+1;
        tau1_l(end+1)=tau1;
        tau2_l(end+1)=tau2;

    end
       plot(tau1_l)
       hold on;
       plot( tau2_l)
       figure
       plot(theta1)

    end