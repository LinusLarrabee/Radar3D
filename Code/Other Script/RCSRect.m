w = 2; %距离向长度
d = 2; %方位向长度
C=3.0e8;  %光速(m/s)
Fc=100e6;  %雷达射频 1.57GHz
Lambda=C/Fc;    %雷达工作波长
k = 2*pi/Lambda;
theta = pi/2*(1e-3:1e-3:1);
fai = pi/12;
deltaCar = pi^(-1)*(k.*d.*w).^2 .*sinc(k.*d.*sin(theta)*cos(fai)).^2 ...
          .*sinc(k*w*sin(theta)*sin(fai)).^2.*0.5.* ...
        ((1-sin(theta)*cos(fai).^2).^(0.5) ...
         +(1-sin(theta)*sin(fai).^2).^(0.5));
     figure(2);
     plot(theta*180/pi,20*log(deltaCar));
     
     ylabel ('RCS \delta');
xlabel ('\theta - 度');