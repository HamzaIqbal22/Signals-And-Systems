% % A.1
% % Set component values:
% R = [1e4, 1e4, 1e4]; C = [1e-6, 1e-6];
% % Determine the coefficients for characteristic equation:
% A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% % Determine characteristic roots:
% lambda = roots(A);
%  
% %The poly command takes in the matrix of roots and  returns the original
% %polynomial equation.
% poly(lambda);

% A.2

%  R = [1e4, 1e4, 1e4]; 
%  C = [1e-6, 1e-6];
%  A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
%  lambda = roots(A);
%  poly(lambda);
% 
% t = [0:0.0005:0.1]; 
% u = @(t) 1.0 * (t>=0);
% h = @(t) (C(1).* exp(lambda(1).* t) + C(2).* exp(lambda(2).*t)).*(u(t));
% 
% plot(t,h(t));
% xlabel ("t");
% ylabel ("h(t)");
% title ('Plot of h()');
% 
% grid;

% A.3
% function[lambda] = lab2(R,C)
%  
% % Determine the coefficients for characteristic equation:
% A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
%  
% % Determine characteristic roots:
% lambda = roots(A);



% B.1
% CH2MP4.m : Chapter 2, MATLAB Program 4
% Script M-file graphically demonstrates the convolution process.figure(1) 
% Create figure window and make visible on screen
% u = @(t) 1.0*(t>=0);
% x = @(t) 1.5*sin(pi*t).*(u(t)-u(t-1));
% h = @(t) 1.5*(u(t)-u(t-1.5))-u(t-2)+u(t-2.5)
% dtau = 0.005;
% tau = -1:dtau:4;ti = 0; 
% tvec = -.25:.1:3.75; 
% y = NaN*zeros(1,length(tvec));
% % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end

% % B.2
% % CH2MP4.m : Chapter 2, MATLAB Program 4
% % Script M-file graphically demonstrates the convolution process.figure(1) 
% % Create figure window and make visible on screen
% u = @(t) 1.0*(t>=0);
% x = @(t) 1.0.*(u(t)-u(t-2));
% h = @(t) 1.0.*(t+1).*(u(t+1)-u(t));
% dtau = 0.005;
% tau = -1:dtau:4;ti = 0; 
% tvec = -.25:.1:3.75; 
% y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end



% % B.2
% % CH2MP4.m : Chapter 2, MATLAB Program 4
% % Script M-file graphically demonstrates the convolution process.figure(1) 
% % Create figure window and make visible on screen
% u = @(t) 1.0*(t>=0);
% x = @(t) 1.0.*(u(t)-u(t-2));
% h = @(t) 1.0.*(t+1).*(u(t+1)-u(t));
% dtau = 0.005;
% tau = -1:dtau:4;ti = 0; 
% tvec = -.25:.1:3.75; 
% y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end

% % B.3
% %Convolution of A
% % CH2MP4.m : Chapter 2, MATLAB Program 4
% % Script M-file graphically demonstrates the convolution process.figure(1) 
% % Create figure window and make visible on screen
% u = @(t) 1.0*(t>=0);
% x = @(t) (u(t-4)-u(t-6));
% h = @(t) 2.0.*(u(t+5) - u(t+4));
% dtau = 0.005;
% tau = -15:dtau:15;ti = 0; 
% tvec = -15:1:15; 
% y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end



% % B.3 (Convolute B:)
% % CH2MP4.m : Chapter 2, MATLAB Program 4
% % Script M-file graphically demonstrates the convolution process.figure(1) 
% % Create figure window and make visible on screen
% u = @(t) 1.0*(t>=0);
% x = @(t) 1.0.*(u(t-3)-u(t-5));
% h = @(t) 2.0.*(u(t+5)-u(t+3));
% dtau = 0.005;
% tau = -15:dtau:15;ti = 0; 
% tvec = -15:1:15; 
% y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end


% %c1
% h1 = @(t) exp(t/5).*u(t);
% h2 = @(t) 4*exp(-t/5).*u(t);
% h3 = @(t) 4*exp(-t).*u(t);
% h4 = @(t) 4*(exp(-t/5) - exp(-t)).*u(t);
% 
% u = @(t) 1.0.* (t>=0);
% t = [-1:0.001:5];
%  
% plot(t, h1(t));
% 
% ylabel("h(t)");
% xlabel("t");
% 
% 
% hold on;
% plot(t,h2(t));
% plot(t,h3(t));
% plot(t,h4(t));
% 
% legend ("h1", "h2", "h3", "h4");
%  
% hold off;

% % C.3
% % h2(t)
% % CH2MP4.m : Chapter 2, MATLAB Program 4
% % Script M-file graphically demonstrates the convolution process.
% % Create figure window and make visible on screen
% u = @(t) 1.0.*(t>=0);
% x = @(t) sin(5.*t).*(u(t)-u(t-3));
% h = @(t) 4*exp(-t/5).*u(t);
% dtau = 0.005; 
% tau = 0:dtau:20;ti = 0; 
% tvec =[0:0.1:20];
% y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
% for t = tvec,
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],"edgecolor","none");
%     xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
%     xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
%     axis([tau(1) tau(end) -1.0 2.0]); grid;
%     pause;
% end






% C.3
% h3(t)
% CH2MP4.m : Chapter 2, MATLAB Program 4
% Script M-file graphically demonstrates the convolution process.
% Create figure window and make visible on screen
u = @(t) 1.0.*(t>=0);
x = @(t) sin(5.*t).*(u(t)-u(t-3));
h = @(t) 4*exp(-t).*u(t);
dtau = 0.005; 
tau = 0:dtau:20;ti = 0; 
tvec =[0:0.1:20];
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec,
    ti = ti+1; % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau); 
    % Trapezoidal approximation of convolution integral
    subplot(2,1,1),plot(tau,h(tau),"k-",tau,x(t-tau),"k--",t,0,"ok");
    axis([tau(1) tau(end) -2.0 2.5]);
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
        [.8 .8 .8],"edgecolor","none");
    xlabel("\tau"); title("h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]");
    c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
    subplot(2,1,2),plot(tvec,y,"k",tvec(ti),y(ti),"ok");
    xlabel("t"); ylabel("y(t) = \int h(\tau)x(t-\tau) d\tau");
    axis([tau(1) tau(end) -1.0 2.0]); grid;
    pause;
end






