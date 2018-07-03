%% QUEST threshold estimation with dynamic and static expressions.
% 2017-11-23 Junpeng Lao @FRIBOURG,CH
clear;
clc;
sca
%% set up and parameters
ScreenX=60;
ScreenY=34;
PixelX=1920;
PixelY=1080;
% setting the size of the stimuli at a distance of 60cm
% scale=3;
scale=(3*ScreenY/ScreenX)*(PixelX/PixelY);
static_noise=1;

Screen('Preference', 'SkipSyncTests', 1);
% setup experiment
subInitials = input('participant''s name: ','s');
subage = input('participant''s age: ','s');
subgender = input('participant''s gender: (1 Female 2 Male) ','s');
newexp = str2double(input('is it a new participant (1) or a continue participant (2)?: ','s'));
randtype = input('type of random: (0 Uniform <--> 5 Half-Mixture) ','s');
mixing = str2double(randtype)/10;
filename = sprintf('%s',subInitials);
% create text file for data and parameters recording
datatxt=[filename '_g' subgender '_a' subage '_r' randtype '_DyExpthr.csv'];
fileID = fopen(datatxt, 'w');
formatSpec = '%d %d %d %d %d %d %d %d %d %d %s %s %s %d\n';
%% Video stimuli related
% video scale
fps1=30; % refresh rate for the movie

% load videos
folder=cd;
VideoMat=load('JOV2013ExpressionStimuli-3.mat');
VideoMat=rmfield(VideoMat,'randomsequences');

expressions=fieldnames(VideoMat);
itemall=eval(['fieldnames(VideoMat.', expressions{1}, ')']);
% counter of stimulus presentation
CounterM = zeros(2, length(expressions), length(itemall));

% load Mask.mat
% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Reseed the random-number generator for each expt.
rng(sum(100*clock));
%% random table
stimulitype = {'Static'; 'Dynamic'};
Nrep = 8; % Each identity repeated Nrep times, meaning that for each 
           % expression and each stimuli type (dynamic or static) subject
           % see 8*Nrep times. 
           % Currently this adds up to 6 exp * 8 id * 8 rep * 2 stim = 768
Nblock = 10;
Nexpress = length(expressions);
Nstimtyp = length(stimulitype);
Nitem = length(itemall);
ii=0;
Nexpall = Nstimtyp*Nexpress*Nitem*Nrep;
exptable = cell(Nexpall,4);
for iexp=1:Nexpress
    for iitem=1:Nitem
        for irep=1:Nrep
            randseed = randi(100000);
            for istim = 1:Nstimtyp
                ii=ii+1;
                exptable(ii,1)=stimulitype(istim);
                exptable(ii,2)=expressions(iexp);
                exptable(ii,3)=itemall(iitem);
                exptable{ii,4}=randseed;
            end
        end
    end
end
seedlist = unique([exptable{:,4}])';
seedlist(:,2)=NaN;
% same pdf of both static and dynamic
pdf_x = linspace(0,1,100);

Rsp_allexp = NaN(Nexpress, Nexpall/Nexpress, 2);

if newexp==1
    randseq=randperm(Nexpall);
    ii2=0;
    save([filename,'.mat'],'ii2','randseq')
else
    try
        load([filename,'.mat']);
    catch
        randseq=randperm(Nexpall);
        ii2=0;
        save([filename,'.mat'],'ii2','randseq')
    end
end
%% set screen
screenNumber=max(Screen('Screens'));
[mainw, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2);
ifi=Screen('GetFlipInterval', mainw);
waitframe=round(1/fps1/ifi);% at refresh rate 60HZ
% hide mouse cursor
HideCursor;
    
Screen('FillRect', mainw, [255/2 255/2 255/2 0]);% black screen
Screen('Flip', mainw);

% % second screen for optimal synchronization purpose.
% white=WhiteIndex(1);
% [w2, wRect2]=Screen('OpenWindow',1, white);%%% ,[0 0 800 600]
%% reaction key
% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

recal=KbName('z');
afk=KbName('p');% Fear
ank=KbName('c');% Anger
dik=KbName('d');% Disgust
hak=KbName('j');% Happy
sak=KbName('t');% Sadness
suk=KbName('s');% Surprise
skipk=KbName('l');% Skip
breakcode=KbName('ESCAPE');% breakkey

% screen size
screenx=0.5*wRect(3);
screeny=0.5*wRect(4);
%% introduction
instru_v = 'The experiment is about to start...';
% Write instruction message for subject, nicely centered in the
% middle of the display, in white color. As usual, the special
% character '\n' introduces a line-break:
Screen('TextSize', mainw, 50);
DrawFormattedText(mainw, instru_v, 'center', 'center', WhiteIndex(mainw));
Screen('Flip', mainw);

% Wait for mouse click:
GetClicks(mainw);

% Clear screen to background color (our 'white' as set at the
% beginning):
Screen('Flip', mainw);

% Wait a second before starting trial
WaitSecs(1.000);
%% main experiment
while ii2<Nexpall
    ii2=ii2+1;
    %         DrawFormattedText(w2, num2str(itrial), 'center', 'center', BlackIndex(w2));
    %         Screen('Flip', w2);
    
    [KeyIsDown, vert, KeyCode]=KbCheck;
    Screen('TextSize', mainw, 50);
    DrawFormattedText(mainw, '+', 'center', 'center', WhiteIndex(mainw)); % fixation
    [VBLTimestamp, startrt2] = Screen('Flip', mainw);
    
    while (GetSecs - startrt2)<.5 %presentation until respond
        [KeyIsDown, endrt2, KeyCode]=KbCheck;
        if KeyCode(breakcode)==1
            break
        end
    end
    
    eval(['videotmp1=VideoMat.' exptable{randseq(ii2),2} '.' exptable{randseq(ii2),3} ';'])
    expid = strcmp(expressions, exptable{randseq(ii2),2});
    recorder = find(isnan(Rsp_allexp(expid,:,1)));
    
    [height1,width1,counts]=size(videotmp1);
    
    if strcmp(exptable{randseq(ii2),1}, stimulitype(1)) == 1 % static
        videotmp = repmat(videotmp1(:,:,end),[1,1,counts]);
    else
        videotmp = videotmp1;
    end
    % videotmp(mask==0)=255/2;
    randseed = exptable{randseq(ii2),4};
    
    % random percentage
    seed_loc = find(seedlist(:,1)==randseed);
    if ~isnan(seedlist(seed_loc,2))
        % this seed was used before, reading the presented percentage so
        % that static and dynamic are using the same noise pattern at the
        % same noise intensity level
        percent = seedlist(seed_loc,2);
    else
        percent = rand();
        if recorder(1)>20 % Uniform until the 20th trials
            history_resp = squeeze(Rsp_allexp(expid,1:(recorder(1)-1),:));
            % smooth
            [b,dev,stats] = glmfit(history_resp(:,1),history_resp(:,2),'binomial','link','logit');
            yfit = glmval(b, pdf_x, 'probit');
            % mixture
            cdf_tosamp = mixing*(yfit-min(yfit))/(max(yfit)-min(yfit))+(1-mixing)*pdf_x';
            % get precentage
            idxnew = find(cdf_tosamp(2:end)>percent);
            if ~isempty(idxnew)
                percent = pdf_x(idxnew(1));
            end
        end
        seedlist(seed_loc,2) = percent;
    end
    
    % noise mixing
    [noisemixedstim, purenoise] = noise_mix(videotmp, static_noise, percent, randseed);
    
    % load first frame
    width2=width1*scale;
    height2=height1*scale;
    
    r1=[screenx,screeny];
    cRect1=SetRect(r1(1)-width2/2, r1(2)-height2/2, r1(1)+width2/2, r1(2)+height2/2);
    
    debtrial=GetSecs;
    currenttime=debtrial;
    oldtime=debtrial;
    echantillon=0;
    tvbl = Screen('Flip', mainw);
    % Playback loop: Runs until end of movie
    for iframe=1:counts
        imdata1=squeeze(noisemixedstim(:,:,iframe));
        tex1=Screen('MakeTexture', mainw, imdata1);
        % Draw the new texture immediately to screen:
        Screen('DrawTexture', mainw, tex1, [], cRect1);
        %             % Draw border
        %             Screen('FrameRect', mainw, [255 255 255],cRect1, 5);
        % Update display:
        tvbl = Screen('Flip', mainw, tvbl + ifi*(waitframe-0.5));

        % Release texture:
        Screen('Close', tex1);
    end
    tex2=Screen('MakeTexture', mainw, squeeze(purenoise(:,:,end)));
    % Draw the new texture immediately to screen:
    Screen('DrawTexture', mainw, tex2, [], cRect1);
    DrawFormattedText(mainw, '?', 'center', 'center', WhiteIndex(mainw)); % fixation
    [VBLTimestamp, startrt2] = Screen('Flip', mainw);
    ACC=0;
    while (GetSecs - startrt2) %presentation until respond
        [KeyIsDown, endrt2, KeyCode]=KbCheck;
        if KeyCode(afk)==1
            resp=1;
            ACC=strcmp(exptable{randseq(ii2),2},'Fear');
            break
        elseif KeyCode(ank)==1
            resp=2;
            ACC=strcmp(exptable{randseq(ii2),2},'Anger');
            break
        elseif KeyCode(dik)==1
            resp=3;
            ACC=strcmp(exptable{randseq(ii2),2},'Disgust');
            break
        elseif KeyCode(hak)==1
            resp=4;
            ACC=strcmp(exptable{randseq(ii2),2},'Happiness');
            break
        elseif KeyCode(sak)==1
            resp=5;
            ACC=strcmp(exptable{randseq(ii2),2},'Sadness');
            break
        elseif KeyCode(suk)==1
            resp=6;
            ACC=strcmp(exptable{randseq(ii2),2},'Surprise');
            break
        elseif KeyCode(skipk)==1
            resp=NaN;
            ACC=0;
            break
        end
    end
    
    Rsp_allexp(expid,recorder(1),:) = [percent, ACC];
    
    rt=endrt2-startrt2;
    % save subject respond
    datatowrite=[ii2 resp ACC rt percent cRect1(1) cRect1(2) cRect1(3) cRect1(4) randseq(ii2) exptable(randseq(ii2),:)];
    fprintf(fileID,formatSpec,datatowrite{:});
    save([filename,'.mat'], 'ii2', 'randseq');
    
    if mod(ii2, round(Nexpall/Nblock))==0
        instru_v = 'The next block will start soon...';
        % Write instruction message for subject, nicely centered in the
        % middle of the display, in white color. As usual, the special
        % character '\n' introduces a line-break:
        Screen('TextSize', mainw, 50);
        DrawFormattedText(mainw, instru_v, 'center', 'center', WhiteIndex(mainw));
        Screen('Flip', mainw);
        
        % Wait for mouse click:
        GetClicks(mainw);
        
        % Clear screen to background color (our 'white' as set at the
        % beginning):
        Screen('Flip', mainw);
        
        % Wait a second before starting trial
        WaitSecs(1.000);
    end
end

%% End experiment
Screen('TextSize', mainw, 50);
KbCheck;
WaitSecs(0.1);
endword = 'Please wait...';
% Write instruction message for subject, nicely centered in the
% middle of the display, in black color. As usual, the special
% character '\n' introduces a line-break:
DrawFormattedText(mainw, endword, 'center', 'center', WhiteIndex(mainw));

% Update the display to show the instruction text:
Screen('Flip', mainw);

% Wait 1 sec
WaitSecs(1);

%%close window
fclose('all');
Screen('Closeall');
ShowCursor;
