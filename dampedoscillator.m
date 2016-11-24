function dampedoscillator(A,x0,k1,m,b,tmax)
c=[0;0.5;0.5;1];
a=[0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0];
w=[1/6 1/3 1/3 1/6];
F=@(t,s)[s(2);-k1/m*s(1)-b/m*s(2)];
s(:,1)=[x0+A;0];
t(1)=0;
dt=0.01;
i=1;
 while t(i)<tmax
        k=zeros(length(s(:,1)),length(c));
            for j=1:length(c)
                k(:,j)=dt*F(t(i)+dt*c(j),s(:,i)+k*a(j,:)');
                
            end
            s(:,i+1)=s(:,i)+k*w';
            t(i+1)=t(i)+dt;
            i=i+1;
 end
 delete datafile.dat;
 
 %For the animation of damped spring system
% figure(1);
%     for i=1:length(s(1,:))
%     clf;
%     plot(s(1,i),0,'s','markersize',40,'markerfacecolor','b');
%     axis([0 max(s(1,:)) -0.5 1]);
%     hold on;
%     x=linspace(0,s(1,i),20);
%     y(1:2:length(x))=0.04;
%     y(2:2:length(x))=-0.04;
%     plot(x,y)
%     plot([0 s(1,i)], [0 0]);
%     pause(0.01);
%     hold off   
%   end
  
  
 %Basically what this does is it creates a data point of displacement and
 %time in an interval of 0.01, these values are then saved into a
 %datafile.dat(incase we need a backup of the data) and we use this file to read the data through 'dlmread'. It
 %read any numeric data as a matrix and plot it's component to get graph of
 %x(t) Vs. t
figure(2);
for i=1:length(s(1,:))
    t(i)=t(i)+dt;
    xvalue=s(1,i)';
    tvalue=t(i)';
    dataset=[tvalue xvalue];
    save datafile.dat dataset -ascii -append;   
end
dampedMatrix=dlmread('datafile.dat'); %dlmread-read ascii delimited file of numeric data as a matrix
displacement=dampedMatrix(:,2);
time=dampedMatrix(:,1);
plot(time,displacement,'r-','linewidth',2);
xlabel('time(t)');
ylabel('displacement x(t)'); 
title('graph of displacement x(t) Vs. time t in damped harmonic oscillator');
end
