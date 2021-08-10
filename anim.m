clear all;
clc;
load sckadata
l1 = 0.4;l2 = 0.5;l3 = 0.4;
%Fixed link length = 0.2;
figure
axis square
line([0 0.2],[0 0],'LineWidth',3,'Color','k'); % creating fixed link line
axis([-0.8 0.8 -0.8 0.8]);
N = size(t,2);
link1 = [];link2 = [];link3 = [];Opos = [];Apos = [];Bpos = [];Cpos = [];timedisplay = [];

phi1 = pcoordsall(3,1);phi2 = pcoordsall(6,1);phi3=pcoordsall(9,1);
xO = 0; yO = 0 ;
xA = l1*cos(phi1);yA = l1*sin(phi1);
xB = xA + l2*cos(phi2);yB = yA + l2*sin(phi2);
xC = 0.2 ; yC = 0;
link1 = line([0 xA],[0 yA],'LineWidth',3);
link2 = line([xA xB],[yA yB],'LineWidth',3);
link3 = line([xB xC],[yB yC],'LineWidth',3);
%slider = rectangle('Position',[xB-0.2,yB-0.1,0.4,0.2],'FaceColor','c');
Opos = rectangle('Position',[-0.05,-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor','r'); %creating circle with at origin with red colour
Apos = rectangle('Position',[xA-0.05,yA-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor','r');
Bpos = rectangle('Position',[xB-0.05,yB-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor','r');
Cpos = rectangle('Position',[xC-0.05,yC-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor','r');
timedisplay = text(0.7,0.7,num2str(t(1))); %showing timer at point(1,1)

for i = 1:10:N
    set(link1,'xdata',[0 xA],'ydata',[0 yA]);
    set(link2,'xdata',[xA xB],'ydata',[yA yB]);
    set(link3,'xdata',[xB xC],'ydata',[yB yC]);
    %set(slider,'Position',[xB-0.2,yB-0.1,0.4,0.2]);
    set(Opos,'Position',[xO-0.0150,yO-0.0150,0.03,0.03]);
    set(Apos,'Position',[xA-0.0150,yA-0.0150,0.03,0.03]);
    set(Bpos,'Position',[xB-0.0150,yB-0.0150,0.03,0.03]);
    set(Cpos,'Position',[xC-0.0150,yC-0.0150,0.03,0.03]);
    set(timedisplay,'Position',[0.7 0.7],'string',num2str(t(i)));
    drawnow(); % Draws all commands at a time
    phi1 = pcoordsall(3,i);phi2 = pcoordsall(6,i);phi3=pcoordsall(9,1);
    xA = l1*cos(phi1);yA = l1*sin(phi1);
    xB = xA + l2*cos(phi2);yB = yA + l2*sin(phi2);
    xC = 0.2 ; yC = 0;
    pause(0.001);
end