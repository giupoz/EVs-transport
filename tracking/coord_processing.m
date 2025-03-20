format long
clear all

% get txt file with ev coordinates
file=uigetfile(['analyzed/' '*.txt']);

% A is a Nx2 matrix with all the coordinates
A=importdata(file);
xpixel=A(:,1);
ypixel=A(:,2);

%conversion: pixels to um (conversion factors 0.084 or 0.108)
X=xpixel*0.084;
Y=ypixel*0.084;

%scatter plot of coordinates
figure(1)
scatter(X,Y)
hold on
axis equal

%starting point maked with green color and ending point marked with red
plot(X(1),Y(1),'g.','MarkerSize',20)
plot(X(end),Y(end),'r.','MarkerSize',20)
title(file)
%number of coordinates
N=length(X);

%neuron shape approximation: coordinates interpolation with grade 1
%polynomial function
pol=polyfit(X,Y,1);
plot(X,pol(2)+pol(1).*X)

%Variant for exception (trajectory almost vertical)
[Ma,i1]=max(Y);
plot(X(i1),Y(i1),'m.','MarkerSize',20)
[mi,i2]=min(Y);
plot(X(i2),Y(i2),'m.','MarkerSize',20)

Xe=[X(i1) X(i2)];
Ye=[Y(i1) Y(i2)];

pol2=polyfit(Xe,Ye,1);
plot(Xe,pol2(2)+pol2(1).*Xe)

%average between the two angular coefficients
m=(pol(1)+pol2(1))/2;

saveas(figure(1),['coordinates_images/' file '_coor.png']);

%%
%sin() and cos() of the angle formed by the neuron with the x-axis (coeff)
coeff=atan(pol(1)); % pick pol(1), pol2(1) or m 
cosen=cos(coeff);
sen=sin(coeff);

%creation of a vector with the distances between the points
dis_x=zeros(N,1);
dis_y=zeros(N,1);
dis_x(2:N)=(X(2:end)-X(1:end-1));
dis_y(2:N)=(Y(2:end)-Y(1:end-1));

%projection of the distances between the points on the neuron
%tangential distance
dis_t = dis_x.*cosen + dis_y.*sen;
dis_t(isnan(dis_t))=0;

dis_t_woz=[];
j=1;
for i=1:N
    if dis_t(i)~=0
        dis_t_woz(j)=dis_t(i);
        j=j+1;
    end
end

%normal distance
dis_n = dis_y.*cosen - dis_x.*sen;
dis_n(isnan(dis_n))=0;

%istant velocity - tangential direction
ivel_t=dis_t./0.5;
%istantant velocity - normal direction
ivel_n=dis_n./0.5;

%tangential velocities vector clean up:
% if the velocity is less than one pixel per unit time (0.084/0.5=0.168)
% set it to zero and save a vector of non-null velocities
ivel_woz=[];
count_zeros=0;
for i=1:N
    if abs(ivel_t(i))<0.168
        ivel_t(i)=0;
        count_zeros=count_zeros+1;
    else
        ivel_woz=[ivel_woz; ivel_t(i)];
    end
end

%calculate the average speed taking into account the signs
AverageSpeed = ones(N,1).*abs(mean(ivel_t));

%absolute progressive pathlength calculation in tangential direction
abs_dis_t=abs(dis_t);
path=0;
abs_pathlength_t=zeros(N,1);
for k=1:N
    path=path+abs_dis_t(k);
    abs_pathlength_t(k)=path;
end

%Compute the pathlength WITH progressive SIGN in tangential direction
%(result: vector of distances of each point from the initial point where
%the last element of the vector is the distance between the inital and the
%final point

path=0;
pathlength_t=zeros(N,1);
for k=1:N
    path=path+dis_t(k);
    pathlength_t(k)=path;
end

%computing pathlength (with sign in normal dir)
path=0;
pathlength_n=zeros(N,1);
for k=1:N
    path=path+dis_n(k);
    pathlength_n(k)=path;
end

%time vector
seq_end=0.5*(N-1);
timev=[0:0.5:seq_end]';

%absolute progressive average velocity vector in tangential direction
abs_AverageSpeed = abs_pathlength_t./timev;
abs_AverageSpeed(isnan(abs_AverageSpeed))=0;

%moving average tangential velocity on window of 5 elem
MA_tg=zeros(N,1);
MA_tg(1)=mean(ivel_t(1:3));
MA_tg(2)=mean(ivel_t(1:4));
for i=3:N-2
    MA_tg(i)=mean(ivel_t(i-2:i+2));
end
MA_tg(N-1)=mean(ivel_t(N-3:N));
MA_tg(N)=mean(ivel_t(N-2:N));
skewness(MA_tg)

%moving average without zeroes
MA_woz=[];
for i=1:N
    if MA_tg(i)~=0
        MA_woz=[MA_woz; MA_tg(i)];
    end
end


%create table
T=table(X,Y,dis_x,dis_y,dis_t,dis_n,ivel_t,ivel_n,abs_pathlength_t,pathlength_t,timev,abs_AverageSpeed,AverageSpeed,MA_tg);
Name=['all_data/',file,'_data.xlsx'];
writetable(T,Name);

%% PLOTS



% figure(2)
% subplot(1,2,1)
% histogram(MA_tg,'BinWidth',0.05);
% title('Moving average velocity over 5 elements frequency')
% xlabel('v')
% ylabel('frequency')
% subplot(1,2,2)
% histogram(MA_woz,'BinWidth',0.05);
% title('Moving average velocity over 5 elements frequency - no zeros')
% xlabel('v')
% ylabel('frequency')
% 
% figure(3)
% subplot(1,2,1)
% histogram(ivel_t,'BinWidth',0.15);
% title('Istant velocity frequency')
% xlabel('v')
% ylabel('frequency')
% subplot(1,2,2)
% histogram(ivel_woz,'BinWidth',0.15);
% title('Istant velocity frequency - no zeros')
% xlabel('v')
% ylabel('frequency')

figure(4)
plot(timev,pathlength_t)
% hold on
% plot(timev,pathlength_n)
ylim([-3 6])
title(['Trajectory',file])
xlabel('time')
ylabel('position')
saveas(figure(4),['path_images/' file '_path.png']);

% figure(5)
% plot(timev,ivel_t)
% title('Instant velocity over time')
% xlabel('time')
% ylabel('v')
% 
% figure(6)
% plot(timev, MA_tg)
% title('Moving average velocity')
% xlabel('time')
% ylabel('v')


figure(7)
histogram(dis_t,'BinWidth',0.05)
title('tangential displacements frequencies')
xlabel('dx')
ylabel('frequency')


figure(8)
histogram(dis_t_woz,'BinWidth',0.05)
title('tangential displacements frequencies (no zeros)')
xlabel('dx')
ylabel('frequency')





