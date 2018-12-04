clear
close all

G = 6.674e-11; % m^3*kg^-1*s^-2
r = 0.5e11; % m
M_sun = 1.989e30;
m1 = 1*M_sun; % kg
m2 = 2*M_sun; % kg
m3 = 1000; % kg
p1_init = [-2e12,0]; % m
p2_init = [2e12,0]; % m
p3_init = [-2e12-r,0]; % m
v1_init = [0,-5000]; % m/s
v2_init = v1_init/-2; % m/s
v3_init = [0,-5000-sqrt(G*m1/r)]; % m/s
prd = 2*pi*sqrt(r^3/(G*m1));


% calculate orbital parameters
[xcom,ycom,~] = COM([p1_init,0;p2_init,0],[m1;m2]);
r1_init = p1_init-[xcom,ycom];
H1 = cross([r1_init,0],[v1_init,0]);
h1 = norm(H1);
r2_init = p2_init-[xcom,ycom];
H2 = cross([r2_init,0],[v2_init,0]);
h2 = norm(H2);

mu = G*(m1+m2);
mu1 = mu*(m2/(m1+m2))^3;
mu2 = mu*(m1/(m1+m2))^3;

e1 = 1/mu1*((norm(v1_init)^2-mu1/norm(r1_init))*r1_init-dot(r1_init,v1_init)*v1_init);
e2 = 1/mu2*((norm(v2_init)^2-mu2/norm(r2_init))*r2_init-dot(r2_init,v2_init)*v2_init);
e = norm(e1); % = norm(e2);
ax1 = h1^2/(mu1*(1-e^2));
ax2 = h2^2/(mu2*(1-e^2));
ax = ax1+ax2;
rp = ax*(1-e);
ra = ax*(1+e);
T = sqrt(4*pi^2*ax^3/mu);
T_yrs = T/(3600*24*365);
fprintf('Eccentricity is %.4f \n',e);
fprintf('Period is %.2f years \n',T_yrs);
fprintf('Periapsis Radius is %.0f km \n',rp/1e3);
fprintf('Apoapsis Radius is %.0f km \n',ra/1e3);

dt = 15000;
max_t = 2400000000;
L = round(max_t/dt)+1;

a1 = zeros(L,2); a2 = zeros(L,2); a31 = zeros(L,2); a32 = zeros(L,2); a3 = zeros(L,2);
v1 = zeros(L,2)+v1_init; v2 = zeros(L,2)+v2_init; v3 = zeros(L,2)+v3_init;
p1 = zeros(L,2)+p1_init; p2 = zeros(L,2)+p2_init; p3 = zeros(L,2)+p3_init;
max_dist = 0;
min_dist = realmax;

for i = 1:L
    
    if i > 1
        v1(i,:) = v1(i-1,:)+a1(i-1,:)*dt;
        v2(i,:) = v2(i-1,:)+a2(i-1,:)*dt;
        v3(i,:) = v3(i-1,:)+a3(i-1,:)*dt;
        p1(i,:) = p1(i-1,:)+v1(i,:)*dt;
        p2(i,:) = p2(i-1,:)+v2(i,:)*dt;
        p3(i,:) = p3(i-1,:)+v3(i,:)*dt;
    end
    r = p2(i,:)-p1(i,:);
    
    [a1(i,1),a1(i,2)] = accel(m2,r);
    [a2(i,1),a2(i,2)] = accel(m1,-r);
    [a31(i,1),a31(i,2)] = accel(m1,p1(i,:)-p3(i,:));
    [a32(i,1),a32(i,2)] = accel(m2,p2(i,:)-p3(i,:));
    a3(i,:) = a31(i,:)+a32(i,:);
    
    d = norm(r);
    if d > max_dist
        max_dist = d;
        imax = i;
    end
    if d < min_dist
        min_dist = d;
        imin = i;
    end

end



% tol = 1000*dt;
% cond1 = (p1(:,1) < p1_init(1) + tol) & (p1(:,1) > p1_init(1) - tol);
% cond2 = (p1(:,2) < p1_init(2) + tol) & (p1(:,2) > p1_init(2) - tol);
% findT = find(cond1 & cond2);
% T = findT(2)*dt; %/(3600*24*365);

% ax = ((mu*T^2)/(4*pi^2))^(1/3);
% ax1 = ax*m2/(m1+m2);
% ax2 = ax*m1/(m1+m2);


% comet(p1(:,1),p1(:,2))
% comet(p2(:,1),p2(:,2))
% plot(p1(:,1),p1(:,2),'r')
% plot(p2(:,1),p2(:,2),'b')
% scatter(p1(1,1),p1(1,2),'ro')
% scatter(p2(1,1),p2(1,2),'bo')

figure()
hold on
axis equal
g = animatedline('Color','r');
h = animatedline('Color','b');
%g = animatedline('MaximumNumPoints',1,'Marker','o','MarkerSize',6,'MarkerFaceColor','r');
%h = animatedline('Color','b','MaximumNumPoints',1,'Marker','o','MarkerSize',10,'MarkerFaceColor','b');
j = animatedline('Color','g');

for k = 1:length(p1)
    addpoints(g,p1(k,1),p1(k,2));
    addpoints(h,p2(k,1),p2(k,2));
    %addpoints(j,p3(k,1),p3(k,2));
    drawnow limitrate
end
scatter(p1(imax,1),p1(imax,2),'*k');
scatter(p2(imax,1),p2(imax,2),'*k','HandleVisibility','off');
scatter(p1(imin,1),p1(imin,2),'ok');
scatter(p2(imin,1),p2(imin,2),'ok','HandleVisibility','off');
% legend('Orbit 1','Orbit 2','Maximum Distance','Minimum Distance')
fprintf('Maximum Distance is %.0f km \n',max_dist/1e3);
fprintf('Minimum Distance is %.0f km \n',min_dist/1e3);
fprintf('Starting speed is %.2f km/s \n',norm(v3_init)/1e3);
fprintf('Ending speed is %.2f km/s \n',norm(v3(end,:))/1e3);





function [ax,ay] = accel(m,r)
G = 6.674e-11; % m^3*kg^-1*s^-2
R = norm(r);
ax = G*m*r(1)/R^3;
ay = G*m*r(2)/R^3;
end