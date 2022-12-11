%% ELE 632 LABORATORY 5: SMAPLING AND DISCRETE FOURIER TRANSFORM

% Ashwin Vimalanathan

%% PART A: DISCRETE FOURIER TRANSFORM AND ZERO PADDING

%% QUESTION 01

No = 10;
Nz = 490;
Nz2 = 100;
Nz3 = 400;
n = (0:1:(No-1));
r = 0:No-1;
r = r- (No/2);
fr = r/No;
df1 = 0:1:Nz-1;
df2 = 0:1:Nz2-1; 
df3 = 0:1:Nz3-1;
x = exp(1i.*2*pi*n.*0.3) + exp(1i.*2*pi*n.*0.33);
x2 =cos(2.*pi.*n.*0.3)+0.5.*cos(2.*pi.*n.*0.4);
X_r = fft(x); 
X_r2 = fft(x2); 

figure;
stem(fr,fftshift(abs(X_r)));
title('Problem A.1: Magnitude of DFT of signal x1[n] when N=10');
xlabel('\omega');
ylabel('|X(\omega)|');
grid on;

figure;
stem(fr,fftshift(abs(X_r2)));
title('Problem A.1: Magnitude of DFT of signal x2[n] when N=10');
xlabel('\omega');
ylabel('|X(\omega)|'); 
grid on;
 
disp('A.1i) The signal x2[n] has a symmetric spectrum and x1[n] does not have a symmetric spectrum. The reasoning for this is because time-domain signal x2[n] contains only cosine terms which are all even and all real valued terms. The x1[n] signal has complex terms and that caused the spectrum to not be spectrum.');

disp('A.1ii) In signal x1[n], the frequency resolution with value 0.1 intervals is not able to allow the visualization of two components placed at 0.3 and 0.33 because they are relatively closely spaced. As a result of this, it is not possible to distinguish both frequency components in the plotted figures. But, in the case of the signal x1[n], the two components in frequency spectrum can be seen distingushed at 0.3 and 0.4 as the components are separated by a spacing of 0.1 units away from each other.');

disp('A.1iii) The other frequency components are visible in the plotted figure for the DFT of the signal because there was not an integer multiple of periods available as well as the components possible occuring near the region of the abrupt corner caused by the sampling which the truncated signal. when a sharp corner occurs a stream of additional components occurs similar to the rect fucntion.');

%% QUESTION 02

X_rr = fft(x,490); 
X_rr2 = fft(x2,490);

figure;
stem(df1,fftshift(abs(X_rr)));
title('Problem A.2: Magnitude of DFT of signal x1[n] padded with 490 zeros when N=10');
xlabel('\omega');
ylabel('|X(\omega)|');
grid on;

figure;
stem(df1,fftshift(abs(X_rr2)));
title('Problem A.2: Magnitude of DFT of signal x2[n] padded with 490 zeros when N=10');
xlabel('\omega');
ylabel('|X(\omega)|');
grid on;

disp('A.2) For the DFT of the signals x1[n] and x2[n], there are only sampling of 10 taken, thoguht the zero paddding creates a length of 500. This indicates there are 500 frequency components evenly spaced between the frequency interval of -0.5 and 0.5. For signal x1[n], there is still not enough information to distinguish between the closely spaced components at 0.3 and 0.33. This is expects because the padding only alters the width of the frequency interval but not the information held behidn these components.');
 
 %% QUESTION 03
N = 100;
n=(0:N-1);
r = n;
f_r = r/N;

x = exp(1i.*2*pi*n.*0.3) + exp(1i.*2*pi*n.*0.33);
X1 = fft(x);
figure;
hold on;
subplot(2,1,1);
stem(f_r - 0.5, fftshift(abs(X1)));
title('Problem A.3: Magnitude of DFT of signal x1[n] when N=100');
xlabel('f_r');
ylabel('|X_1 (f_r)|');

x2 =cos(2.*pi.*n.*0.3)+0.5.*cos(2.*pi.*n.*0.4);
X2 = fft(x2);
subplot(2,1,2);
stem(f_r - 0.5, fftshift(abs(X2)));
title('Problem A.3: Magnitude of DFT of signal x2[n] when N=100');
xlabel('f_r');
ylabel('|X_1 (f_r)|');
hold off;

disp('A.3) The DFT of the signal remains symmetric and there are now two components located at the points 0.3 and 0.4 present in the DFT of the signal');
 
 %% QUESTION 04
 
N=100;
n=(0:99);
x = exp(1i.*2*pi*n.*0.3) + exp(1i.*2*pi*n.*0.33);
x1_padded = [x,zeros(1,400)];
X1_padded = fft(x1_padded);
f_r = (0:length(x1_padded)-1)/length(x1_padded);
figure;
hold on;
subplot(2,1,1);
stem(f_r - 0.5, fftshift(abs(X1_padded)));
title('Problem A.4: Magnitude of DFT of signal x1[n] padded with 400 zeros when N = 10');
xlabel('f_r');
ylabel('|X_1padded (f_r)|');
axis([-0.5 0.5 0 110]);

n=(0:99);
x2 =cos(2.*pi.*n.*0.3)+0.5.*cos(2.*pi.*n.*0.4);
x2_padded = [x2,zeros(1,400)];
X2_padded = fft(x2_padded);
f_r = (0:length(x2_padded)-1)/length(x2_padded);
subplot(2,1,2);
stem(f_r - 0.5, fftshift(abs(X2_padded)));
title('Problem A.4: Magnitude of DFT of signal x2[n] padded with 400 zeros when N = 10');
xlabel('f_r');
ylabel('|X_2padded (f_r)|');

disp('A.4) There is an integer multiple of periods available currently, so therefore the DFT of the signals is definitely improved.');

 %% PART B: SAMPLING
 
 %% QUESTION 01
 
clear;
close all;

load chirp.mat;
filename =  'chirp.wav';
audiowrite(filename,y,Fs);
clear y Fs;
[y,fs] = audioread(filename);
N_0 = length(y); %Sample length
T_0 = N_0/fs;%Total time 
T = 1/fs;%Sampling Interval

%% QUESTION 02

t = linspace(0,T_0,N_0);
figure;
plot(t,y);
title('Problem B.2: Chirp Audio Signal plotted with respect to time');
xlabel('Time(s)');
ylabel('Audio Signal');
%% QUESTION 03

Y = fft(y);

f = (-N_0/2:N_0/2-1)*100/N_0;

figure;
subplot(2,1,1);
plot(f,abs(fftshift(Y)));
title('Problem B.3: The DFT of signal');
xlabel('f_r');
ylabel('|Y(f_r)|');
subplot(2,1,2);
plot(f,angle(fftshift(Y)));
xlabel('f_r');
ylabel('\angle Y(f_r)');

%% QUESTION 04

y1 = downsample(y,2);
N_01 = length(y1); %Sample length
T_01 = N_01/fs;%Total time 
T1 = 1/fs;%Sampling Interval

%% QUESTION 05

t = linspace(0,T_01,N_01);
figure;
plot(t,y1);
title('Problem B.5: Audio signal sub-sampled at rate of 2');
xlabel('Time(s)');
ylabel('Audio Signal');

%% QUESTION 06

Y1 = fft(y1);
f = (-N_01/2:N_01/2-1)*100/N_01;
figure;
subplot(2,1,1);
plot(f,abs(fftshift(Y1)));
title('Problem B.6: The DFT of signal');
xlabel('f_r');
ylabel('|Y1(f_r)|');
subplot(2,1,2);
plot(f,angle(fftshift(Y1)));
xlabel('f_r');
ylabel('\angle Y1(f_r)');

disp('B.6i) The spectrum of y is a rather larger frequency range, 8kHz, in comparison to the spectrum of y1 which has a range of 4kHz because of the sub-sampling thats occured. In addition, the order of the frequency components are re-arranged due to the aliasing.');

disp('B.6ii) To compare the enlarged parts of the signals, its quite evident that the spectrum of the original signal is definitely more closely and tightly packed than in comparison to the sub-sampled signal.');

%% QUESTION 07

sound(y,fs);

disp('B.7) From listening to the sub-sampled ad original signal, the sub-sampled signal sounds to closely resemble a bird chirping, although there is more distortion occurring, this is possible due to the frequency components overlapping due to aliasing.');
%% QUESTION 08

y2 = downsample(y,5);
N_02 = length(y2); 
T_02 = N_02/fs;
T2 = 1/fs;
t = linspace(0,T_02,N_02);

figure;
plot(t,y2);
title('Problem B.8: Audio Signal subsampled at rate of 5');
xlabel('Time(s)');
ylabel('Audio Signal');
Y2 = fft(y2);
f = (-N_02/2:N_02/2-1)*100/N_02;

figure;
subplot(2,1,1);
plot(f,abs(fftshift(Y2)));
title('Problem B.8: The DFT of the signal');
xlabel('f_r');
ylabel('|Y2(f_r)|');
subplot(2,1,2);
plot(f,angle(fftshift(Y2)));
xlabel('f_r');
ylabel('\angle Y2(f_r)');

disp('B.8) The audio signal is not distinguishable. The spectrum of the signal has a much smaller frequency range compared to the original signal spectrum and lots of aliasing is occurring which results in a heavily distorted sound.');

%% PART C: FILTER DESIGN

%% QUESTION 01

clear;

file = 'handel.wav';
[y,Fs] = audioread(file);
audio = y; 
DFT_audio = fftshift(fft(audio)); 
half = Fs/2; 
t = 0:1:length(DFT_audio)-1; % Time x-axis
t = t/10000;
f = linspace(0, Fs, length(y)); % Frequency x-axis
f = f-half;
H = abs(f) < 2000;
H = transpose(H);
filtered_audio = H.*DFT_audio;

figure;
subplot(2,1,1);
plot(t,audio);
title('Audio Signal in Time Domain');
xlabel('Time(s)');
grid on;   
subplot(2,1,2);
plot(f,abs(DFT_audio)); 
title('Audio Signal in Frequency Domain');
xlabel('Frequency(Hz)');
grid on;

figure;
plot(f,abs(H));
title('2kHz Lowpass Filter');
xlabel('Frequency(Hz)');
grid on;
    
figure;
subplot(2,1,1);
plot(f,abs(filtered_audio));
title('Audio Signal in Frequency-Domain from -2kHz to 2kHz');
xlabel('Frequency(Hz)');
grid on;
subplot(2,1,2);
plot(t,real(ifft(fftshift(filtered_audio))));
title('Filtered Audio Signal in Time-Domain');
xlabel('Time(s)');
grid on;
    
%% QUESTION 02

sound(real(ifft(fftshift(filtered_audio))),Fs);

%% QUESTION 03


H2 = ~(abs(f) >= 16 & abs(f) <= 256);% Bandpass filtering out between 16-256
H2 = transpose(H2);
filtered_audio2 = DFT_audio.*H2;% Bass frequencies filtered out

figure;
plot(f,abs(H2));
title('Bass Filter between 16-256 Hz');
xlabel('Frequency(Hz)');
grid on;

figure;
subplot(2,1,1);
plot(f,abs(filtered_audio2));
title('Filtered Audio Signal in the Frequency-Domain');
xlabel('Frequency(Hz)');
grid on;
subplot(2,1,2);
plot(t,real(ifft(fftshift(filtered_audio2))));
title('Filtered Audio Signal in the Time-Domain');
xlabel('Time(s)');
grid on;

sound(real(ifft(fftshift(filtered_audio2))),Fs);    

disp('C.3) The sound is distorted and does not represent the signal properly. This system is not realizable because of the extremely narrow range of frequencies(16-256) have been removed and the filter becomes quite unstable.');

%% QUESTION 04


H3 = abs(f) >= 2048 & abs(f) <= 16384; % 2048-16384 
H3 = transpose(H3);
H3 = H3.*0.25;
filtered_audio3 = DFT_audio+(DFT_audio.*H3);

figure;
plot(f,real(H3));
title('Treble Filter between 2048-16384 Hz');
xlabel('Frequency(Hz)');
grid on;

figure;
subplot(2,1,1);
plot(f,abs(filtered_audio3));
title('Amplified Audio Signal in the Frequency-Domain');
xlabel('Frequency(Hz)');
grid on;
subplot(2,1,2);
plot(t, real(ifft(filtered_audio3)));
title('Amplified Audio Signal in the Time-Domain');
xlabel('Time(s)');
grid on;

sound(real(ifft(fftshift(filtered_audio3))),Fs); 

%% QUESTION 05

disp('The scaling property of the DFT was utilized.');

 