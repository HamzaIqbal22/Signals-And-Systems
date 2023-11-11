% %A.1

% %Creating the function
% f = @(t) exp(-t).*cos(2*pi*t);
% 
% %Setting the value of t
% t= (-2:2);
% 
% plot (t,f(t));
% xlabel('t');
% ylabel('f(t)');
% grid;
% 
% title('Figure 1.46');


% %Creating the function
% f = @(t) exp(-t).*cos(2*pi*t);
%  
% %Setting the value of t
% t = (-2:0.01:2);
% 
% plot(t,f(t));
% xlabel('t');
% ylabel('f(t)');
% grid;
% title('Figure 1.47');


% %A.2
% 
% %Creating the function
% f = @(t) exp(-t);
% 
% %Setting the value of t
% t = (-2:2);
% 
% plot(t,f(t));
% xlabel('t');
% ylabel('f(t)');
% grid;
% title('A.2. Graph');

% %B.1
% 
% p = @(t) 1.0.*((t>=0)&(t<1));
% t = (-1:0.01:2);
% plot(t,p(t));
% xlabel('t');
% ylabel('p(t) = u(t) - u(t-1)');
% axis ([-1 2 -.1 1.1]);
% title ('Figure 1.50');

%B.2
 
% %Creating the function
% 
% u = @(t) 1.0.*(t>=0);
% p = @(t) u(t)-u(t-1);
% r = @(t) t.*p(t);
% n = @(t) r(t) + r(-t+2);
% 
% %Setting the value of t
% t = (-2:0.01:2);
% plot(t,r(t));
% xlabel('t');
% ylabel('r(t) = t*p(t)');
% axis ([-1 2 -.1 1.1]);
% 
% plot(t,n(t));
% xlabel('t');
% ylabel('n(t) = r(t) + r(-t +2)');

%B.3
 
% u = @(t) 1.0.*(t>=0);
% p = @(t) u(t)-u(t-1);
% r = @(t) t.*p(t);
% n = @(t) r(t) + r(-t+2);
% n1 = @(t) n(0.5*t);
% n2 = @(t) n1(t+0.5);
% t = (-1:0.01:5);
% plot(t,n1(t),'-k',t,n2(t),':k');
% xlabel('t');
% ylabel('n1(t) & n2(t)');
% axis ([-1 5 -.1 1.1]);

%B.4
 
% u = @(t) 1.0.*(t>=0);
% p = @(t) u(t)-u(t-1);
% r = @(t) t.*p(t);
% n = @(t) r(t) + r(-t+2);
% 
% n1 = @(t) n(0.5*t);
% n2 = @(t) n1(t+0.5);
% n3 = @(t) n(t+0.25);
% n4 = @(t) n3(0.5.*t);
% 
% t = (-1:0.01:5);
% 
% plot(t,n3(t),'-k',t,n4(t),':k');
% xlabel('t');
% ylabel('n3(t) & n4(t)');
% axis ([-1 5 -.1 1.1]);


% %C.1
% f = @(t)exp(-2*t).*cos(4*pi*t);
% g = @(t) f(t).* u(t);
% u = @(t) 1.0.* (t>=0);
% t = (-2:0.01:2);
% plot(t,g(2*t+1));
% xlabel('t');
% ylabel('g(2+1)');
% grid;
% title ('Figure 1.51');
 
 
 
% C.2 
% s = @(t) g(t+1);
% t = (0:0.01:4);
% plot (t, s(t));
% grid;
% xlabel ('t');
% ylabel ('s(t) = g(t+1)');

% C.3
% u = @(t) 1.0 .* (t>=0);
% t = (0:0.01:4);
% 
% for alpha = 1:2:7
%     s = @(t) exp (-2).* exp(-alpha.*t).*cos(4*pi*t).*u(t);
%     plot (t, s(t));
%     grid;
%     xlabel ('t');
%     ylabel ('s(t)');
%     hold on;
% end   
%     
% legend ('Alpha 1', 'Alpha 3', 'Alpha 5', 'Alpha 7')


%D.2.a)
% load('ELE532_Lab1_Data.mat')
% num_rows = size(B,1);
% num_cols = size(B,2);
%  
% for i=1:1:num_rows
%     for j=1:1:num_cols
%         if(abs(B(i,j)) < 0.01)
%             B(i,j)=0;
%         end
%     end
% end
 
%D.2.b)
% load('ELE532_Lab1_Data.mat')
% B([abs(B) >= 0.01])= 0;
%  
% %D2.c) i)
% tic
% load('ELE532_Lab1_Data.mat')
%  
% num_rows = size(B,1);
% num_cols = size(B,2);
%  
% for i=1:1:num_rows
%     for j=1:1:num_cols
%         if(abs(B(i,j)) < 0.01)
%             B(i,j)=0;
%         end
%     end
% end
%  
% fprintf('\nD2:\nFor part a:\n');
% toc
%  
% %D2.c) ii)
%  
% tic
% load('ELE532_Lab1_Data.mat')
% B([abs(B) >= 0.01])= 0;
% fprintf('\nD2: \nFor part b: Elapsed time is :\n');
% toc
%  
% %D.3:
load('ELE532_Lab1_Data.mat')
 
%Copying data array x_audio into audio.
audio = x_audio;
 
num_rows = size(audio,1);
num_cols = size(audio,2);
 
number_of_zeros = 0;
 
for i = 1: num_rows
    for j = 1: num_cols
        if(abs(audio(i,j) == 0))
            number_of_zeroes = number_of_zeroes +1;
        end
    end
end
    
    
 % # of elements in 0
 fprintf("\n" + number_of_zeroes);
 
 % Play sound
 sound(audio,8000)

