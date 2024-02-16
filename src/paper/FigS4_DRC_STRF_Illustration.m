% This script re-creates Figure S4 from van den Berg*, Wong*, Houtak,
% Williamson, Borst (2024) iScience

load('.\Data\Neuro\DRC\M30_S04_DRC.mat','Stm') % stimulus matrix
load('.\Data\Neuro\DRC\M30_S04_DRC_SpkT.mat','SpkT') % Trial-aligned spike times
load('.\Data\Neuro\DRC\M30_DRC_PSTH.mat','PSTH') % PSTH
 
close all

FontSize = 16;
%% Get stimulus mask (spectrogram) and generate waveform
Mask = Stm.Mask{1};
Onset = Stm.Onset{1};
Durs = Stm.DurMask{1};

K = size(Mask,2);
minFreq = 2000; maxFreq = 64000;
Freq = logspace(log10(2000),log10(64000),K);
freq = Freq./1000;
TimeStep = 5; %ms
Gain = [Freq',ones(K,2)];
Ref = [1,70,100];
Fs = 192000;
Spk = 1;

% generate DRC waveform (whole)
[Snd,Pul] = genDRC(TimeStep,Freq,Onset,Durs,Mask,Fs,Spk,Gain,Ref);

% get a concatenated version of the waveform (first 10 seconds)
snipDur = 10; %s
Sound = Snd(1:(snipDur*Fs));
tt = (1:length(Sound))./Fs;

% ---  uncomment the following line to save the waveform as a .wav file ---
% audiowrite('DRC_sample.wav',Sound./max(abs(Sound)),Fs)
% -------------------------------------------------------------------------
%% Figure 8A-D Plot Spectrogram, Waveform, Raster, PSTH
fig = figure;
fig.Position = [50,100,1000,750];

% Stimulus spectrogram
ax1 = subplot(5,1,[1,2]);
imagesc(Mask',[20,70])
set(gca,'YDir','normal')
set(gca,'ytick',1:12:K,'yticklabels',freq(1:12:K));
ylabel('Freq (kHz)')
colormap(sqrt(flipud(gray)))
% colormap([1,1,1;jet])
% cb = colorbar('Position',[.92,.45,.015,.45]);
cb = colorbar('Position',[.92,.65,.015,.25]);
cb.Limits = [25,70];
cb.Label.String = 'Intensity (dB SPL)';
set(ax1,'FontSize',FontSize)
set(ax1,'xtick',(0:0.5:2).*1000./TimeStep,'xticklabels',(0:0.5:2));
% xlabel(ax1,'Time (s)')

xlim([0,2000/TimeStep])

% Stimulus waveform
ax2 = subplot(5,1,3);
plot(tt,Sound,'k');
set(ax2,'FontSize',FontSize)
set(ax2,'xtick',(0:0.5:2),'ytick',[]);
% xlabel(ax2,'Time (s)')
xlim([0,2])

u1 = 1;
u2 = 11;

ax3 = subplot(5,1,4);

[~,~,~] = plotraster(ax3, SpkT(:,u1),(1:20)', [0,0,1],0,1);hold on;
[~,~,~] = plotraster(ax3, SpkT(:,u2),(1:20)', [1,0,0],0,1);hold on;

set(ax3,'FontSize',FontSize)
xlim(ax3,[0,2])
set(ax3,'xtick',(0:0.5:2));
% xlabel(ax3,'Time (s)')
ylim(ax3,[.5,20.5]);ylabel('Rep #');

ax4 = subplot(5,1,5);
% errorbar(mean(PSTH{1,u},2),std(PSTH{1,u},[],2));
stairs(mean(PSTH{1,u1},2),'Color','b'); hold on;
stairs(mean(PSTH{1,u2},2),'Color','r');
set(ax4,'FontSize',FontSize)
set(ax4,'xtick',(0:0.5:2).*1000./TimeStep,'xticklabels',(0:0.5:2));
xlabel(ax4,'Time (s)')
ylabel(ax4,'r (spk/bin)');
xlim([0,2000/TimeStep])
ylim([0,2]);

fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'.\Figures\Links\FigS4A-D.pdf')


% %% Stimulus scaling (input non-linearity)
% stimulus = Mask;
% stimulus=stimulus-max(max(stimulus));
% stimulus=10.^(stimulus./20);
% 
% figure;
% uInt = unique(Mask);
% uStim = uInt-max(uInt);
% uStim=10.^(uStim./20);
% plot(uInt,uStim,'-o','LineWidth',2,'MarkerSize',9);
% xlabel('intensity (dBSPL)')
% ylabel('''stimulus strength''')
% set(gca,'FontSize',16)
% xlim([25,70])
% axis square
% 
% %% Stimulus characteristics
% % histogram of durations
% figure('Position',[500,400,500,200])
% histogram(Durs*TimeStep,20:5:220)
% xlim([0,Inf])
% ylabel('Num. tones');
% ylim([0,2000])
% set(gca,'FontSize',16)
% xticks([0:50:220])
% xlabel('Tone duration (ms)');
% 
% % histogram of intensities
% figure('Position',[500,400,500,200])
% histogram(Onset,[22.5:5:72.5])
% ylabel('Num. tones');
% set(gca,'FontSize',16)
% xticks([20:10:75]);
% xlabel('Tone intensity (dB)');

%% STRF, PRF/CGF
B = load('.\Data\Neuro\DRC\ASD_6_0x5_0x5\Results_M30_S04_J8_M10_ASD_amp.mat'); % baseline data
S = load('.\Data\Neuro\DRC\ASD_6_0x5_0x5\Results_M30_S10_J8_M10_ASD_amp.mat'); % salicylate data


%% Figure 8E-J plot PRF & CGF baseline

units = [u1,u2];

fig = figure;
fig.Position = [1120,100,600,750];

for uu = 1:2
    prf = flip(B.results(units(uu)).full_rank_sparse_rep.wtf',2);
    J = size(prf,2); K = size(prf,1);
    strf = reshape(B.results(units(uu)).strf.ww(2:end),J,K)';
    cgf = B.results(units(uu)).full_rank_sparse_rep.wtauphi';
    M = size(cgf,2); N = (size(cgf,1)-1)/2;
    cgfmask = ones(size(cgf));
    cgfmask(N+1,1) = 0;

    cmax = max(abs([strf(:)]));%max(abs([strf(:);prf(:)]));
    subplot(3,2,uu)
    imagesc(fliplr(strf),[-cmax,cmax]); colormap(jet);
    colorbar
    setPRFAxis(gca,'reverse')
    set(gca,'FontSize',FontSize)
    title('STRF')

    cmax = max(abs([prf(:)]));%max(abs([strf(:);prf(:)]));
    subplot(3,2,2+uu)
    imagesc(fliplr(prf),[-cmax,cmax]); colormap(jet);
    colorbar
    setPRFAxis(gca,'reverse')
    set(gca,'FontSize',FontSize)
    title('PRF')

    cmax = max(abs([cgf(:)]));
    subplot(3,2,4+uu)
    imagesc((cgf), ...
        'AlphaData',cgfmask, ...
        [-cmax,cmax]); colormap(jet); 
    colorbar
    setCGFAxis(gca,'reverse')
    set(gca,'FontSize',FontSize)
    title('CGF')
    
end

fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'.\Figures\Links\FigS4E-J.pdf')
%% functions

function [Snd,Pul] = genDRC(TimeStep,freqs,Onset,Durs,Mask,Fs,Spk,Gain,Ref)
    disp('Generating DRC signal...');
    
    lgth = size(Onset,1);
    
    TimeStepSamp = round(TimeStep*Fs/1000);
    maxDur = max(Durs,[],'all');
    NSamp = maxDur*TimeStepSamp;
    BlockDurSamp = (lgth+maxDur-1) * TimeStepSamp;
    GateTime = 5;
 
    % Generate pulse (for first frequency)
    Pul = reshape(repmat((Mask(:,1) > -Inf)',TimeStepSamp,1) , 1, []);

    % Generate sound
    Snd = zeros(1,BlockDurSamp);
    % staggered generation
    for j = 1:maxDur
        tempSnd = zeros(1,BlockDurSamp);
        for i = j:maxDur:lgth
%             % -- debugging --
%             if (~mod(i,50)); disp([num2str(i) ' / ' num2str(lgth) ' timestep']); end
%             % ---------------
            offset = (i - 1) * TimeStepSamp;
            levels = Onset(i,:);
            durations = Durs(i,:);
            UDur = unique(durations);
            UDur = UDur(UDur~=0);
            chord = zeros(1,NSamp);
            for d = UDur
                idx = durations == d;
                newChord = genChord(d*TimeStep,freqs(idx),levels(idx),Fs, GateTime, Spk, Gain, Ref);
                if length(newChord) < NSamp
                    newChord(end+1:NSamp) = 0;
                end
                chord = chord + newChord;
            end
            tempSnd(offset+(1:NSamp)) = chord;
        end
        Snd = Snd + tempSnd;
    end
    
end
    
function chord = genChord(Dur,freqs,levels,Fs, GateTime, Spk, Gain, Ref)
% generate multi-tone signal of Dur with sampling frequency Fs
% 
% INPUTS:
%   Dur         = duration in ms
%   freqs       = vector of frequencies
%   levels      = vector of levels
%   Fs          = sampling frequency
%   GateTime    = duration of cosine gate
%   Gain, Ref   = gain to apply for setup


RefdB		=	Ref(1,2);
DACmax		=	Ref(1,3);

sel			=	Gain(:,3) == Spk;
Gain		=	Gain(sel,:);

%-- Carrier time axis --%
cSamp		=	round( (Dur/1000)*Fs );		%-- # samples		--%
cTime		=	( (1:cSamp)-1 ) / Fs;		%-- Time axis [sec]	--%


Amp		=	getamp(Gain,freqs,levels,RefdB,DACmax);

chord	=	Amp * sin(2*pi*freqs'*cTime);



if( size(chord,2) > 2 )
    Nenv			=	round( GateTime * 1e-3 * Fs );
    chord			=	envelope(chord',Nenv)';
end

end


function Amp = getamp(GainTable,Freq,Lvl,RefdB,DACmax)

    %-- Select interpolate gain from the calibration table --%
    uSpk = unique(GainTable(:,3));
    nSpk = length(uSpk);

    Gain2D = nan(nSpk,length(Freq));
    for s = 1:nSpk
        Spk = uSpk(s);
        GT = GainTable(GainTable(:,3) == Spk,:);
        Gain2D(s,:) = interp1(GT(:,1),GT(:,2),Freq(:));
    end

    Gain = mean(Gain2D,1); % if both speakers 

    Amp		=	(10.^((Lvl-RefdB)/20)) ./ Gain;

    if ( Amp >= DACmax )
        if( idx+1 > size(GainTable,1) )
            Gain	=	min([ mean( GainTable(idx-1:idx,2) ) 9]);
        else
            Gain	=	min([ mean( GainTable(idx-1:idx+1,2) ) 9]);
        end
        Amp		=	(10.^((Lvl-RefdB)/20)) ./ Gain;
        warning(['Used adjacent gains to calibrate f= ' num2str(round(Freq)) ' Hz & clipped to 9 V if necessary.'])
    end
    
end


function Sig = envelope(Sig, NEnv)
% Create a smooth on- and offset envelope for an auditory signal
%
% function SIG = ENVELOPE (SIG, NENV)
%
% .. Dr. P ...

if (length(NEnv) == 1); NEnv = [NEnv,NEnv];end

SigLen = size(Sig,1);

if (SigLen < 2*NEnv)

  disp ('-- ERROR: Envelope length greater than signal');

else

  Env1 = ( sin(0.5*pi*(0:NEnv(1))/NEnv(1)) ).^2;
  Env2 = flip( sin(0.5*pi*(0:NEnv(2))/NEnv(2)) ).^2;
  head = 1:(NEnv(1)+1);
  tail = (SigLen-NEnv(2)):SigLen;

  for i=1:size(Sig,2)
    Sig(head,i) = Env1' .* Sig(head,i);
    Sig(tail,i) = Env2' .* Sig(tail,i);
  end
end
end


%% LOCAL Functions
function setPRFAxis(h,xdir)
    if nargin < 2; xdir = 'normal';end
    K = 61; J = 8;
    freq = logspace(log10(2),log10(64),K);
    delay = 5.*(0:J);
    set(h,'YDir','normal','XDir',xdir)
    set(h,'ytick',1:12:K,'yticklabel',freq(1:12:K));
    xlabel('Delay (ms)');
    set(h,'xtick',1:2:J+1,'xticklabel',delay(1:2:J+1));
    ylabel('Freq (kHz)');
end

function setCGFAxis(h,xdir)
    if nargin < 2; xdir = 'normal';end
    N = 12; M = 10;
    phi = (-N:N)./12;
    tau = 5.*(0:M);
    set(h,'YDir','normal','XDir',xdir)
    set(h,'ytick',1:N/2:2*N+1,'yticklabel',phi(1:N/2:end));
    set(h,'xtick',1:5:M+1,'xticklabel',tau(1:5:M+1));
    xlabel('Tau (ms)');
    ylabel('Phi (oct)');
%     rectangle('Position',[0.5,N+0.5,1,1],'EdgeColor',[.5 .5 .5],'FaceColor',[.5 .5 .5],'LineWidth',0.05)
%     rectangle('Position',[0.5,N+0.5,1,1],'EdgeColor','none','FaceColor',[.5 .5 .5])
end