function [noiseimg_, purenoise] = noise_mix(videotmp, static_noise, percent, randseed)
rng(randseed)
% noise mixing function
[height1,width1,counts]=size(videotmp);

image = videotmp;%-mean(videotmp(:));
imagefft = fft2(image);
imageamp = abs(imagefft);
imageph = angle(imagefft);
rmsc = 23;%std2(image);

% rmsc = rmsc*.75;
rand_samp = rand(height1,width1,counts);
randphase = 2*pi*rand_samp-pi;

if static_noise
    randphase1 = repmat(randphase(:,:,end), 1, 1, counts);
else
    randphase1 = randphase;
end
Sphase = percent.*sin(imageph)+(1-percent).*sin(randphase1);
Cphase = percent.*cos(imageph)+(1-percent).*cos(randphase1);
noisetmp = atan(Sphase./Cphase);
noisetmp(Cphase<0)=noisetmp(Cphase<0)+pi;
noisetmp( (Sphase < 0)&(Cphase > 0) )=noisetmp( (Sphase < 0)&(Cphase > 0) )+2.*pi;

noiseimg = real(ifft2(imageamp.*exp(sqrt(-1)*noisetmp)));
noiseimg2 = noiseimg.*sqrt((rmsc.^2)/(std2(noiseimg)).^2);
noiseimg_ = imresize(uint8(noiseimg2+127.5), [height1*2, width1*2]);

Sphase2 = sin(randphase1);
Cphase2 = cos(randphase1);
noisetmp2 = atan(Sphase2./Cphase2);
noisetmp2(Cphase2<0)=noisetmp2(Cphase2<0)+pi;
noisetmp2( (Sphase2 < 0)&(Cphase2 > 0) )=noisetmp2( (Sphase2 < 0)&(Cphase2 > 0) )+2.*pi;

noiseimgp = real(ifft2(imageamp.*exp(sqrt(-1)*noisetmp2)));
noiseimgp2 = noiseimgp.*sqrt((rmsc.^2)/(std2(noiseimgp)).^2);
purenoise = imresize(uint8(noiseimgp2+255/2), [height1*2, width1*2]);