% % a.3
% function [D1,D2,D3]=lab3(d,n)
%  
% D1 = (1/2)*(abs(n)==3)+(1/4)*(abs(n)==1);
% D2 = (1/n*pi).*(sin(pi*n/2));
% D3 = (1/n*pi).*(sin(pi*n/4));
% 
% if (d == 1) 
% D1 = (1/2)*(abs(n)==3)+(1/4)*(abs(n)==1);
% end    
% 
% if (d == 2) 
% D2 = (1/n*pi).*(sin(pi*n/2));
% end
% 
% if (d == 3) 
% D3 = (1/n*pi).*(sin(pi*n/4));
% end
% 
% end


% a.4: (a)
% x1(t)
% clf;
% n = (-5:5);
%  
% D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi ))+(1./pi.*n).*sin((3+n).*pi)+(1./(2.*n.*pi).*sin((1+n).*pi))+(1./(2.*n.*pi).*sin((1-n).*pi));
%  
% subplot(1,2,1);
% stem(n,abs(D_n),'.k');
% xlabel('n'); 
% ylabel('|D_n|');
%  
% subplot(1,2,2); 
% stem(n,angle(D_n),'.k');
% xlabel('n'); 
% ylabel('\angle D_n [rad]');


% %x2(t)
% clf;
% n = (-5:5);
%  
% D_n = (1 ./(n.*pi).*sin((n.*pi)./2));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% 
% %x3(t)
% clf;
% n = (-5:5);
%  
% D_n = (1./(n.*pi).*sin((n.*pi)./4));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% %a.4
% %x1(t)
% clf;
% n = (-20:20);
%  
% D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi ))+(1./pi.*n).*sin((3+n).*pi)+(1./(2.*n.*pi).*sin((1+n).*pi))+(1./(2.*n.*pi).*sin((1-n).*pi));
%  
% subplot(1,2,1);
% stem(n,abs(D_n),'.k');
% xlabel('n'); 
% ylabel('|D_n|');
%  
% subplot(1,2,2); 
% stem(n,angle(D_n),'.k');
% xlabel('n'); 
% ylabel('\angle D_n [rad]');
% 
% %x2(t)
% clf;
% n = (-20:20);
%  
% D_n = (1 ./(n.*pi).*sin((n.*pi)./2));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% 
% %x3(t)
% clf;
% n = (-20:20);
%  
% D_n = (1./(n.*pi).*sin((n.*pi)./4));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% %a.4 (c)
% %x1(t)
% clf;
% n = (-50:50);
%  
% D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi ))+(1./pi.*n).*sin((3+n).*pi)+(1./(2.*n.*pi).*sin((1+n).*pi))+(1./(2.*n.*pi).*sin((1-n).*pi));
%  
% subplot(1,2,1);
% stem(n,abs(D_n),'.k');
% xlabel('n'); 
% ylabel('|D_n|');
%  
% subplot(1,2,2); 
% stem(n,angle(D_n),'.k');
% xlabel('n'); 
% ylabel('\angle D_n [rad]');
% 
% 
% %x2(t)
% clf;
% n = (-50:50);
%  
% D_n = (1 ./(n.*pi).*sin((n.*pi)./2));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% 
% %x3(t)
% clf;
% n = (-50:50);
%  
% D_n = (1./(n.*pi).*sin((n.*pi)./4));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% %a.4 (d)
% %x1(t)
% clf;
% n = (-500:500);
%  
% D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi ))+(1./pi.*n).*sin((3+n).*pi)+(1./(2.*n.*pi).*sin((1+n).*pi))+(1./(2.*n.*pi).*sin((1-n).*pi));
%  
% subplot(1,2,1);
% stem(n,abs(D_n),'.k');
% xlabel('n'); 
% ylabel('|D_n|');
%  
% subplot(1,2,2); 
% stem(n,angle(D_n),'.k');
% xlabel('n'); 
% ylabel('\angle D_n [rad]');
% 
% %x2(t)
% clf;
% n = (-500:500);
%  
% D_n = (1 ./(n.*pi).*sin((n.*pi)./2));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% %x3(t)
% clf;
% n = (-500:500);
%  
% D_n = (1./(n.*pi).*sin((n.*pi)./4));
%  
% subplot(1,2,1); stem(n,abs(D_n),'.k');
% xlabel('n'); ylabel('|D_n|');
% subplot(1,2,2); stem(n,angle(D_n),'.k');
% xlabel('n'); ylabel('\angle D_n [rad]');
% 
% 
%%A5
% function [D] = lab3(Dn)
% n=-500:500;
% D=0.25*sinc(n/4);
% t=[-300:1:300];
% w=pi*0.1;
% x=zeros(size(t));
% for i = 1:length(n)
%  x=x+D(i)*exp(j*n(i)*w*t);
%  't'
% end
% 
% figure(5);
% plot(t,x,'k')
% xlabel('t(sec)');
% ylabel('x(t)');
% axis([-300 300 -1 2]);
% title('Fourier Coefficients Reconstructed');
% grid;




