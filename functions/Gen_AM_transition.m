% This script creates envelope of AM stimulus in van den Berg*, Wong* et al
% the final output is the vector total_mag, with timing in TT (or x)

%% Generation of envelope

% display 
startT = -0.2; endT = 2.2;
Fs = 1/5.12e-6;         % sampling rate
showGates = 1; % whether the different envelopes/gates are shown
showGates = 0; % whether the different envelopes/gates are shown


% stimulus parameter
    Md = 1;  % ; modulation depth
    Mf = 16; % Hz; modulation frequency
    stimDur = 2; %s

    % rise and fall gates of overall stimulus
    RiseTime = 5e-3; %s
    FallTime = 5e-3; %s
    % transition between UM and AM
    transStart = 1;%s; transition start
    transDur = 0.02; %s; transition duration

    % (cosine)phase for the AM transition start
    ModPower    = 1+0.5*Md^2;
    phi = 2*pi-acos(sqrt(ModPower)-1); % phi in radian [0,pi] 

TT = startT:1/Fs:endT;
x = TT;
preIdx = x< transStart ;
transIdx = x>=transStart & x <= transStart+transDur;
postIdx = x > transStart+transDur;
onsetIdx = x >= 0 & x <= RiseTime;
offsetIdx = x >= stimDur - FallTime & x <= stimDur;


% AM envelope
AM_mag = nan(size(x)); AM_mag(transIdx|postIdx) = 1+Md*cos(2*pi*Mf*(x(transIdx|postIdx)-transStart)+phi);

% AM transition "gate"
AM_gate = zeros(size(x)); 
AM_gate(transIdx)= sin((x(transIdx)-1)/transDur *pi/2);
AM_gate(postIdx) = 1;

% UM "envelope"
UM_mag = NaN(size(x)); UM_mag(preIdx|transIdx) = sqrt(ModPower);

% UM transition "gate"
UM_gate = zeros(size(x)); 
UM_gate(preIdx) = 1; 
UM_gate(transIdx) = cos((x(transIdx)-1)/transDur *pi/2); 
UM_gate(postIdx) = 0;

% overall cos2 gating
onOff_gate = ones(size(x));
onOff_gate(x<0) = 0; % before stim
onOff_gate(x>stimDur) = 0; % after stim
onOff_gate(onsetIdx) = 0.5-0.5*cos(x(onsetIdx)/ RiseTime * pi);
onOff_gate(offsetIdx) = 0.5-0.5*cos((x(offsetIdx)-stimDur)/ FallTime * pi);

% combining the results
total_mag = nan(size(x));
total_mag(preIdx) = UM_mag(preIdx);
total_mag(postIdx) = AM_mag(postIdx);
total_mag(transIdx) = sqrt( (AM_mag(transIdx).*AM_gate(transIdx)).^2 + (UM_mag(transIdx).*UM_gate(transIdx)).^2);
total_mag = total_mag .* onOff_gate;

%%
fig = figure;
fig.Position = [100,50,750,450];


for ss = 1:3
switch ss
    case 1
        subplot(2,2,[1,2]);
    case 2
        subplot(2,2,3);
    case 3
        subplot(2,2,4);
end
hold on;
if showGates; p4 = plot(x,UM_mag,'-b','LineWidth',2);end
if showGates; p2 = plot(x,AM_mag,'-r','LineWidth',2);end
if showGates; p6 = plot(x,onOff_gate,':','Color',[0.5,0.5,0.5],'LineWidth',3);end
p1 = plot(x,total_mag,'-k','LineWidth',1);
if showGates; p5 = plot(x,UM_gate,'--','Color',[0,0,0.7],'LineWidth',1.5);end
if showGates; p3 = plot(x,AM_gate,'--','Color',[0.7,0,0],'LineWidth',1.5);end
switch ss
    case 1
        xlim([startT,endT]);
    case 2
        xlim([0.95,1.10]);
    case 3
        xlim([1.90,2.05]);
end

xlabel('Time (s)');
end


fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi

%% LOCAL FUNCTIONS
function    [mainSnd,sideSnd,TT,FF]    =   gennoise(Par, guardBand,NoiseRamp)
% function to generate different kind of noise

%% Basic parameters
FSam        =   Par.Fs;             % Hz; Sampling Rate
Dur         =   Par.SigDur / 1000;         % ms->s; Duration of Signal
FreqBand    =   Par.FreqBand;
useLogDensity = Par.logdensity;
SidePower    =  1; % power

if nargin < 5|| isempty(guardBand)||length(FreqBand) == 1 ; guardBand = 0;end % Hz
if nargin < 6|| isempty(NoiseRamp)||length(FreqBand) == 1 ; NoiseRamp = 0;end % Hz

%% Specify Frequency Bands
if length(FreqBand) == 1 
    f1      =   (FreqBand(1,1)*1000);      %kHz -> Hz lower bound
    f2      =   (FreqBand(1,1)*1000);      %kHz -> Hz upper bound
else
    f1      =   (FreqBand(1,1)*1000);      %kHz -> Hz lower bound
    f2      =   (FreqBand(1,2)*1000);      %kHz -> Hz upper bound
end

MaskFreqs   = [f1 - guardBand,...
             f1 + NoiseRamp + guardBand,...
             f2 - NoiseRamp - guardBand,...
             f2 + guardBand];

if (max(MaskFreqs) > FSam / 2)
   warning('Nyquist is violated!!!') 
end

%% Derived parameters
nSamp       =   round(FSam*Dur);           % Number of samples in signal
dF          =   1/Dur;              % frequency resolution

%% Use Calibration
% DACmax      =   Ref(1,3);
% RefdB       =   Ref(1,2);

% Generate Signal
%% Select frequncy band
FF          =   0:dF:dF*(nSamp-1);  % Freq Axis

mainIdx     =   (FF >= f1) & (FF <= f2);
sideIdx     =   (FF >= MaskFreqs(1)) & (FF < MaskFreqs(2))...
                    | (FF > MaskFreqs(3)) & (FF <= MaskFreqs(4));

mainN       =   sum([mainIdx]);           % number of freq samples in band
sideN       =   sum([sideIdx]);       % number of freq samples in band

%% Generate random phased spectrum

mainXX          =   zeros(1,nSamp);   % initialize with zeros
mainXX(mainIdx) =   exp(2*pi*rand(1,mainN)*1i); % euler form - flat spectrum

if NoiseRamp > 0
    rampIdx1 = (FF >= f1) & (FF < f1 + NoiseRamp);
    rampIdx2 = (FF > f2 - NoiseRamp) & (FF <= f2);
    mainXX(rampIdx1) = mainXX(rampIdx1) .* sqrt(linspace(0,1,sum(rampIdx1)));
    mainXX(rampIdx2) = mainXX(rampIdx2) .* sqrt(linspace(1,0,sum(rampIdx2)));
end

sideXX          =   zeros(1,nSamp);   % initialize with zeros
sideXX(sideIdx) =   sqrt(SidePower) * exp(2*pi*rand(1,sideN)*1i); % euler form - flat spectrum

mainN       =   sum(abs(mainXX).^2);  % number of freq samples in band
sideN       =   sum(abs(sideXX).^2);  % number of freq samples in band
totalN      =   mainN + sideN;  % number of freq samples in band

%% log vs linear power density scaling
if useLogDensity > 0
    rawMS = rms(sideXX+mainXX)^2;

    logDenIdx   = FF >  useLogDensity; % anything above useLogDensity will be 1/f
    flatIdx     = FF <= useLogDensity; % anything <= useLogDensity will be flat
    
    % apply 1/f scaling
    sideXX_logden(logDenIdx)  =   sideXX(logDenIdx) ./ sqrt(FF(logDenIdx));
    mainXX_logden(logDenIdx)  =   mainXX(logDenIdx) ./ sqrt(FF(logDenIdx));
    
    % apply flat scaling
    sideXX_logden(flatIdx)  =   sideXX(flatIdx) ./ sqrt(useLogDensity);
    mainXX_logden(flatIdx)  =   mainXX(flatIdx) ./ sqrt(useLogDensity);

    % rescale total power (RMS)
    newMS       =   rms(sideXX_logden+mainXX_logden)^2;
    scale       =   sqrt(rawMS / newMS);
    sideXX      =   sideXX_logden .* scale;
    mainXX      =   mainXX_logden .* scale;
end
%% apply calibration
% % ToneSPL		=	Int - 10 * log10(totalN);	%-- Each component contributes Lvl - 10*log10(# components) to the overall level --%
% % -- scale so that main sound is the target intensity --
% ToneSPL		=	Int - 10 * log10(mainN);	%-- Each component contributes Lvl - 10*log10(# components) to the overall level --%
% 
% sideXX(sideIdx)      =   sideXX(sideIdx).*getamp(Gain,FF(sideIdx),ToneSPL,RefdB,DACmax);
% mainXX(mainIdx)      =   mainXX(mainIdx).*getamp(Gain,FF(mainIdx),ToneSPL,RefdB,DACmax);

%% generate t-domain signal

mainSnd         =   fft(mainXX);  % fft;
sideSnd         =   fft(sideXX);  % fft;
TT              =   [0:nSamp-1]./FSam;    % time axis vector
end

function    [outSnd,modSnd,sideSnd,scramSide]    =   local_scramble(FSam, InSnd, Mf, Md)
% function to generate different kind of noise
% input InSnd should be complex
if isempty(Mf) || isempty(Md) || Mf == 0 || Md == 0; outSnd = InSnd; return;end

%% Basic parameters
nSamp       =   length(InSnd);
TT          =   [1:nSamp]./FSam;

%% [to be added] check Mf to be integer multiple of 1/Dur

%% Modulate
phi = 2*pi-acos(sqrt((1+0.5*Md^2))-1); % phi in radian [0,pi] 
    % Note: "instantaneous power" of starting phase is matched to scrambled version 
modSnd         =   InSnd .* (1+Md*cos(2*pi*(Mf*TT) + phi));

% modSnd      =   InSnd .* (1+Md*cos(2*pi*(Mf*TT)));

sideSnd     =   modSnd - InSnd;
sideXX      =   ifft(sideSnd);

%% Scramble phase of sideband
scramXX     =   sideXX .* exp(2*pi*1i*rand(size(sideXX)));

%% generate t-domain signal
scramSide   = fft(scramXX);
outSnd      = InSnd + scramSide;

end
