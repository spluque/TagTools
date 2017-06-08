function [newcue, wtype, lev, freqr, stype] = whistle_class(tag, cues, doisave)
%routine to manually classify whistles, mark the type, and indicate relative level
%inputs are:
%   tag     tag name string eg 'gm09_219a'
%   cues    a 2 column vector of [starttimes duration] in seconds since tagon for events in
%           the tagout.  If they are not all whistles (eg some are start of
%           clicking, etc.) it is ok - you will have the chance to delete them from
%           the whistle count during processing.
%   doisave 1 if you want to save output to a mat file, 0 if not.
%
%outputs are:
%   newcue  an adjusted vector of whistle times in cst [start end dur]
%   wtype   a vector of numbers to indicate whistle types.  For gm09_219a,
%           whistle types are 1A-C - short flat whistle followed by upsweep, sounds
%           like a redwing blackbird, sometimes with upsweeps clustered (1A),
%           sometimes not(B), sometimes preceded by a click burst (C); 2, simple tonal
%           curve, maybe with harmonics; 3, extended, very FM whistle; 4, short
%           chirp, usually fermata-shaped; 5, pulses, grunts and squeaks; 6,
%           ultrasonic; 7, other
%   lev     indicator of relative level of whistle - 1 = faint, 2 = quite
%           clearly audible, 3 = high level
%   freqr    frequency range of the whistle including harmonics - [min max]
%            in Hertz
%   stype   a descriptor of the sound
%           SDR, August 2009

%**************************************************************************
%  Begin taogout-specific information
%**************************************************************************
path = 'E:\tag\data'; %base path for tag data
settagpath('audio',path,'cal',[path '\cal'], 'raw',[path '\raw'],'prh','E:\tag\data\prh'); %set dtag path structure
[dummy,afs,rcue] = tagwavread(tag,100,0.1); %to read in afs
%**************************************************************************
%  End taogout-specific information
%**************************************************************************

%**************************************************************************
%  Preallocate space for outputs
%**************************************************************************
newcue = zeros(length(cues),2);
wtype = zeros(length(cues),2);
lev = zeros(length(cues),1);
freqr = zeros(length(cues),2);
stype = cell(length(cues),1);
load (['whistleclass_' tag(1:2) tag(6:end)], 'newcue', 'wtype', 'lev', 'freqr');%get previous data to add on to if classification was started before
%**************************************************************************
%  Extract a Whistle clip
%**************************************************************************
amp = 150; %amplification factor for playing sound
bk = 2048;%spectrogram block size
st = cues(:,1); %start times of marked events
et = cues(:,1) + cues(:,2); %end times of marked events
for k = 233:length(st)
    %load in whistle data for one whistle
    [w,afs,rcue] = tagwavread(tag,st(k)-0.1,et(k) - st(k)+0.3); %read in whistle audio
    figure(1); clf; spectrogram(w(:,1),hamming(bk),bk/2,bk,afs, 'yaxis'); ylim([0,80000]);
    pl = input('Play sound? 1= yes; 0 = no. ');
    if pl == 1
        sound(w*amp,afs);
    end
    keepthisone = input('Is this a whistle, or another sound you want included in the whistle classification? 0=no 1=yes   ');
    stype{k} = input('Enter a word to describe this sound (e.g. call, soc, etc...): ', 's');
    if keepthisone == 1
    disp('click on: whistle start time, whistle end time, and lowest & highest frequencies that are part of the whistle');
    inps = ginput(4); sortrows(inps,1); hp = min(inps(:,2)); lp = max(inps(:,2)); wst = floor(min(inps(:,1))*afs); wed = ceil(max(inps(:,1))*afs);
    wd = w(wst:wed,:);
    newcue(k,1:2) = [wst./afs + st(k) - 0.15, wed./afs + st(k) - 0.15];%start and end in cst
    freqr(k,1:2) = [hp,lp];
    wtype(k,1) = input('Which whistle type? Enter 1-7:  ');
    if wtype(k,1) == 1
        wtype(k,2) = input('Which whistle subtype? Enter 1-3 (blackbird/pod, clean, clicks):  ');
    end
    lev(k) = input('What is the level? Enter 1-3 (faint, audible, loud):  ');
    if doisave == 1
        save(['whistleclass_' tag(1:2) tag(6:end)], 'newcue', 'wtype', 'lev', 'freqr'); 
    end
    clear w inps 
    end
end
newcue(newcue(:,1) == 0,:) = []; %get rid of entries that are not whistles
wtype(wtype(:,1) ==0,:) = [];
lev(lev==0) = [];
freqr(freqr(:,1)==0,:) = [];
if doisave == 1
    save(['whistleclass_' tag(1:2) tag(6:end)], 'newcue', 'wtype', 'lev', 'freqr'); 
end
    