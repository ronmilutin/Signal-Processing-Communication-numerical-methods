clear 
%q1c%
figure
[bii,aii]= butter(8,0.86,'s');
[Hii,wii]= freqs(bii,aii,1024);
[bblt,ablt]= butter(7,0.92,'s');
[Hblt,wblt]= freqs(bblt,ablt,1024);
plot(wii,abs(Hii),"red")
hold on
plot(wblt,abs(Hblt),"blue")
legend('II','BLT')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
axis([0 2 0 1.2])
hold off
%q1d,e%
[bzii,azii]= impinvar(bii,aii,1)
[bzblt,azblt]= bilinear(bblt,ablt,1)
[Haii,waii]= freqz(bzii,azii,1024);
[Hablt,wablt]= freqz(bzblt,azblt,1024);
%q1f%
figure
plot(waii,20*log(abs(Haii(1:length(waii)))),"red")
hold on
plot(wablt,20*log(abs(Hablt(1:length(wablt)))),"blue")
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
legend('II','BLT')
%% Q2a,b %%
wrect=rectwin(22)
whamming=hamming(22)
whann=hann(22)
wblackman=blackman(22)
wbartlett=bartlett(22)

brect_fir=fir1(21,0.3,wrect)
bhamming_fir=fir1(21,0.3,whamming)
bhann_fir=fir1(21,0.3,whann)
bblackman_fir=fir1(21,0.3,wblackman)
bbartlett_fir=fir1(21,0.3,wbartlett)

[H_rect,w_rect]= freqz(brect_fir, 1, 1024);
[H_hamming,w_hamming]= freqz(bhamming_fir, 1, 1024);
[H_hann,w_hann]= freqz(bhann_fir, 1, 1024);
[H_blackman,w_blackman]= freqz(bblackman_fir, 1, 1024);
[H_bartlett,w_bartlett]= freqz(bbartlett_fir, 1, 1024);

figure(), subplot(2,1,1)
plot(w_rect, 20*log(abs(H_rect)),"LineWidth",1)
title('rect')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
grid on
subplot(2,1,2)
plot(w_rect, angle(H_rect),"LineWidth",1)
xlabel('\Omega [rad]')
ylabel('Phase [rad]')

figure(), subplot(2,1,1)
plot(w_hamming, 20*log(abs(H_hamming)),"LineWidth",1)
title('hamming')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
grid on
subplot(2,1,2)
plot(w_hamming, angle(H_hamming),"LineWidth",1)
xlabel('\Omega [rad]')
ylabel('Phase [rad]')

figure(), subplot(2,1,1)
plot(w_hann, 20*log(abs(H_hann)),"LineWidth",1)
title('hann')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
grid on
subplot(2,1,2)
plot(w_hann, angle(H_hann),"LineWidth",1)
xlabel('\Omega [rad]')
ylabel('Phase [rad]')

figure(), subplot(2,1,1)
plot(w_blackman, 20*log(abs(H_blackman)),"LineWidth",1)
title('blackman')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
grid on
subplot(2,1,2)
plot(w_blackman, angle(H_blackman),"LineWidth",1)
xlabel('\Omega [rad]')
ylabel('Phase [rad]')

figure(), subplot(2,1,1)
plot(w_bartlett, 20*log(abs(H_bartlett)),"LineWidth",1)
title('bartlett')
xlabel('\Omega [rad]')
ylabel('Amplitude [DB]')
grid on
subplot(2,1,2)
plot(w_bartlett, angle(H_bartlett),"LineWidth",1)
xlabel('\Omega [rad]')
ylabel('Phase [rad]')

ideal = zeros(1,512);
for i =1:512
    if i<(153)
        ideal(i) = 1;
    end
end
figure(), subplot(2,1,1)
plot(0:0.006:3.066,ideal);
ylim([0 pi]);
title('ideal')
ylabel('Magnitude [DB]')
xlabel('\Omega [rad]')
grid on
subplot(2,1,2)
ideal_phase = angle(ideal);
x_axis = 0:0.006:pi;
plot(x_axis(1:512),ideal_phase)
ylabel('Phase [rad]')
xlabel('\Omega [rad]')
%Q2c%
wvtool(wblackman)
wvtool(wbartlett)
wvtool(whamming)
wvtool(whann)
wvtool(wrect)
%Q2d%
%%
filters = [H_rect H_hamming H_hann H_blackman H_bartlett];
eror = zeros(5,1);
for i = 1:5
    minus = (abs(ideal - filters(i))).^2;
    eror(i) = (1/(2*pi))*(trapz(0:0.006:3.066,minus))
end
%%
 %parameters of kaiser filter

M = 29; %number of samples
beta = 3.5035;
wc = 19/50*pi;
wp = 0.3*pi;
ws = 0.46*pi;
 %generate kaiser filter

 w_kaiser = kaiser(M+1,beta);
 n = 0:M;
 figure()
 plot(n,w_kaiser)
 xlabel('Sample index n');
ylabel('Amplitude');
title('Kaiser window');
xa=impinvar(w_kaiser,n);

% calculate of hd[n] - ideal filter in the time domian 

hd = sin(wc*(n-M/2))./(pi*(n-M/2));
figure
plot(n, hd);
xlabel('Sample index n');
ylabel('Amplitude');
title('hd');

% calculate & generate h[n] in the time domaim
figure;
h_n = hd.*w_kaiser'; % multiple two different vectors
plot(n,h_n,"r");
xlabel('sample index n');
ylabel('Amplitude');
title('h[n] kaiser windowing');

% generate abs(h(w)) and angle of h(w) in DB

N = 512; % Number of frequency points
[H_kaiser, w] = freqz(h_n, 1, N); % moving h[n] to the frequency domain

figure;
plot(w,20*log10(abs(H_kaiser)));
xlabel('Frequency [rad]');
ylabel('Magnitude [dB]');
title('Magunitude Kaiser');
figure;
H_phase = angle(H_kaiser);
plot(w,H_phase, 'LineWidth', 2);
xlabel('Frequency [rad]');
ylabel('Phase [dB]');
title('Phase Kaiser');

% generate the error approximation
figure;
E_w = (1 - abs(H_kaiser)).*(w<wp) + abs(H_kaiser).*(w>ws);
plot(w, E_w, 'LineWidth', 2);
hold on;
xlabel('Frequency [rad]');
ylabel('Error');
title('approximation error - Kaiser');
























