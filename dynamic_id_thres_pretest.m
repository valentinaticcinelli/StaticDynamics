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
filename = sprintf('%s',subInitials);
% create text file for data and parameters recording
datatxt=[filename '_DyIdthr_pretest.csv'];
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
Nrep = 1;  % Each identity repeated Nrep times, meaning that for each 
           % expression and each stimuli type (dynamic or static) subject
           % see 8*Nrep times. 
           % Currently this adds up to 6 exp * 8 id * 8 rep * 2 stim = 768
Nblock = 1;
Nexpress = length(expressions);
Nstimtyp = length(stimulitype);
Nitem = length(itemall);

logger = zeros(Nexpress, Nitem);
ii2=0;
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

% name list
idname = {'Antoine', 'Sophie', 'Didier', 'Fanny', 'Gregoire', 'Helene', 'Joseph', 'Katia'};
% Male1, Female1, Male2, Female2, Male3, Female3, Male4, Female4
ak=KbName('a');% Antoine
sk=KbName('s');% Sophie
dk=KbName('d');% Didier
fk=KbName('f');% Fanny
gk=KbName('g');% Gregoire
hk=KbName('h');% Helene
jk=KbName('j');% Joseph
kk=KbName('k');% Katia

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
while sum(logger(:))<(Nexpress*Nitem)
    if ii2>(Nexpress*Nitem)*2
        break
    end
    [a,b]=find(logger~=1);
    randseq=randperm(length(a));
    t_exp = expressions{a(randseq(1))};
    t_id = itemall{b(randseq(1))};
    t_stim = stimulitype{(rand<.5)+1};
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
    
    eval(['videotmp1=VideoMat.' t_exp '.' t_id ';'])
    idexp = strcmp(itemall, t_id);
    
    [height1,width1,counts]=size(videotmp1);
    
    if strcmp(t_stim, stimulitype(1)) == 1 % static
        videotmp = repmat(videotmp1(:,:,end),[1,1,counts]);
    else
        videotmp = videotmp1;
    end
    % videotmp(mask==0)=255/2;
    randseed = randi(100000);
    percent = .75;
    
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
        if KeyCode(ak)==1
            resp=1;
            ACC=strcmp(t_id,'Male1');
            break
        elseif KeyCode(sk)==1
            resp=2;
            ACC=strcmp(t_id,'Female1');
            break
        elseif KeyCode(dk)==1
            resp=3;
            ACC=strcmp(t_id,'Male2');
            break
        elseif KeyCode(fk)==1
            resp=4;
            ACC=strcmp(t_id,'Female2');
            break
        elseif KeyCode(gk)==1
            resp=5;
            ACC=strcmp(t_id,'Male3');
            break
        elseif KeyCode(hk)==1
            resp=6;
            ACC=strcmp(t_id,'Female3');
            break
        elseif KeyCode(jk)==1
            resp=7;
            ACC=strcmp(t_id,'Male4');
            break
        elseif KeyCode(kk)==1
            resp=8;
            ACC=strcmp(t_id,'Female4');
            break
        elseif KeyCode(skipk)==1
            resp=NaN;
            ACC=0;
            break
        end
    end
    if ACC==0
        logger(:, b(randseq(1))) = 0;
        DrawFormattedText(mainw, 'X', 'center', r1(2), [255, 0, 0]);
        DrawFormattedText(mainw, idname{idexp}, 'center', r1(2)+50, WhiteIndex(mainw));
        Screen('Flip', mainw);
        WaitSecs(1)
        Screen('Flip', mainw);
    else
        logger(a(randseq(1)), b(randseq(1))) = 1;
    end
        
    rt=endrt2-startrt2;
    % save subject respond
    datatowrite=[ii2 resp ACC rt percent cRect1(1) cRect1(2) cRect1(3) cRect1(4) 0. {t_exp}, {t_id}, {t_stim}, {randseed}];
    fprintf(fileID,formatSpec,datatowrite{:});
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
