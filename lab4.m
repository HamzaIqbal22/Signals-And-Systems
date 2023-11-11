%a.1
% u = @(t) 1.0*(t>=0);
% x = @(t) 1.*(u(t)-u(t-10));
% h = @(t) 1.*(u(t)-u(t-10)); % same as x(t)
% dtau = 0.005; 
% tau = -1:dtau:25;ti = 0; 
% tvec = -.25:.1:25;
% y = NaN*zeros(1,length(tvec)); 
% % Pre-allocate memory
% for t = tvec
%     ti = ti+1; % Time index
%     xh = x(t-tau).*h(tau); 
%     lxh = length(xh);
%     y(ti) = sum(xh.*dtau); 
%     % Trapezoidal approximation of convolution integral
%     subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
%     axis([tau(1) tau(end) -2.0 2.5]);
%     patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
%         [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
%         [.8 .8 .8],'edgecolor','none');
%     xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
%     c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
%     subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
%     xlabel('t'); ylabel('z(t) = \int h(\tau)x(t-\tau) d\tau');
%     axis([tau(1) tau(end) -1.0 20.0]); grid;
%     pause
% end

%a.2
%  Xf = fft (x);
%  Zf= Xf.*Xf; 


% %a.3
% N = 100; PulseWidth = 10;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
%  
% stairs (t,x); grid on; axis([-10,110,-0.1,1.1])
% 
% f = [-(N/2):1:(N/2)-1]*(1/N);
%  
% subplot(211); plot(f, fftshift( abs(Zf))); title ("Magnitude "); xlabel ("w"); ylabel ("|z(w)|"); grid on;
% subplot(212); plot(f, fftshift( angle (Zf))); title ("Angle "); xlabel ("w"); ylabel (" <z(w)"); grid on;

 
% %a.4
% N = 100; PulseWidth = 10;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% f = [-(N/2):1:(N/2)-1]*(1/N);
% Xw = fft(x);
% Zw = Xw.*Xw;
% Zt = ifft(Zw);
% w = 2.*pi.*f;
% figure(1);
% subplot(211)
% plot([0,t+1],[0,Zt]); axis([-1 25 -1.0 15]); grid on;
% title('Z(t) = x(t)*x(t) in Time Domain');
% xlabel('t'); ylabel('Z(t)');
% subplot(212); plot(w, fftshift(Zw)); grid on;
% title('Z(jw) = X(jw)X(jw) in Frequency Domain');
% xlabel('w'); ylabel('Z(jw)');
% figure(3);
% %fftshift to centre the double sided spectrum on the freq origin
% %to display magnitude
% subplot(211); plot(w, fftshift(abs(Zw))); grid on;
% title('Magnitude'); xlabel('w'); ylabel('|z(jw)|');
% %display phase
% subplot(212); plot(w, fftshift(angle(Zw))); grid on;
% title('Angle'); xlabel('w'); ylabel('<z(jw)');


%a.5
% N = 100; PulseWidth = 5;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% Xf = fft (x);
% f = [-(N/2):1:(N/2)-1]*(1/N);
% subplot(211); plot(f,fftshift( abs(Xf)));title('Magnitude'); xlabel('w'); ylabel('|z(jw)|'); grid on;
% subplot(212); plot(f,fftshift(angle(Xf)));title('Angle');xlabel('w'); ylabel('<z(jw)'); grid on;
 

% N = 100; PulseWidth = 25;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% Xf = fft (x);
% f = [-(N/2):1:(N/2)-1]*(1/N);
% subplot(211); plot(f,fftshift( abs(Xf)));title('Magnitude'); xlabel('w'); ylabel('|z(jw)|'); grid on;
% subplot(212); plot(f,fftshift(angle(Xf)));title('Angle');xlabel('w'); ylabel('<z(jw)'); grid on;
 

%a.6 (w+)
% N = 100; PulseWidth = 10;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% w = x.*exp(1i*(pi/3)*t);
% Wf = fft (w);
% f = [-(N/2):1:(N/2)-1]*(1/N);
% subplot(211); plot(f,fftshift( abs(Wf)));title ("Magnitude "); xlabel ("w"); ylabel ("|z(w)|"); grid on;
% subplot(212); plot(f,fftshift(angle(Wf))); title ("Angle "); xlabel ("w"); ylabel ("<z(w)");grid on;

%a.6 (w-)
% N = 100; PulseWidth = 10;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% w = x.*exp(-1i*(pi/3)*t);
% Wf = fft (w);
% f = [-(N/2):1:(N/2)-1]*(1/N);
% subplot(211); plot(f,fftshift( abs(Wf)));title ("Magnitude "); xlabel ("w"); ylabel ("|z(w)|"); grid on;
% subplot(212); plot(f,fftshift(angle(Wf))); title ("Angle "); xlabel ("w"); ylabel ("<z(w)");grid on;

% a.6 (wc) 
% N = 100; PulseWidth = 10;
% t = [0:1:(N-1)];
% x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
% w = x.*cos((pi/3)*t);
% Wf = fft (w);
% f = [-(N/2):1:(N/2)-1]*(1/N);
% subplot(211); plot(f,fftshift( abs(Wf)));title ("Magnitude "); xlabel ("w"); ylabel ("|z(w)|"); grid on;
% subplot(212); plot(f,fftshift(angle(Wf))); title ("Angle "); xlabel ("w"); ylabel ("<z(w)");grid on;
 
% %B
% 
%  
% o1 = osc(2000, 80000);
% y1 = xspeech .* o1;
% y2 = conv (y1, hChannel);
%  
% figure (1);
% subplot (211); MagSpect (xspeech); title ( ' xspeech Audio Signal ');
% subplot (212); MagSect (hChannel); title ('hChannel');
% 
% 
% figure (2);
% MagSpect (y2);
% title ('Transmitted Signal over Channel');
%  
% o2 = osc(2000, 80810);
% y3 = y2 .* o2;
% y4 = conv(y3, hLPF2000);
%  
% figure (3);
% MagSpect (y4);
% title ('Decoded Signal from Channel Output');
% xaudio = real(y4);
%  
% figure (1);
% subplot (211); MagSpect(xspeech);grid on;  title('xspeech'); 
%  
% subplot (412); MagSpect(hLPF2000);grid on;  title('hLPF2000'); 
% axis ([-1.55*10^4 1.55*10^4 -150 50 ]);
%  
% subplot (413); MagSpect(hLPF2500);grid on;  title('hLPF2500'); 
% axis ([-1.55*10^4 1.55*10^4 -150 50 ]);
%  
% subplot (412); MagSpect(hChannel);grid on;  title('hChannel'); 
% axis ([-1.55*10^4 1.55*10^4 -150 50 ]);
