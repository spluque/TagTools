function [lunges, lunge_times] = detect_lunges(Aw, fs, fc, ref, depth, speed, plot)
% Detect the occurences of a lunges from triaxial acceleration data.
%
% Inputs:
%   Aw is the nx3 acceleration matrix with columns [ax ay az] within the 
%       frame of the whale. Acceleration can be in either g or m/s^2.
%   fs is the sampling rate in Hz of the acceleration signals.
%   fc (optional) specifies the cut-off frequency of a low-pass filter. The
%       filter cut-off frequency is with respect to 1=Nyquist frequency.
%       Filtering adds no group delay. If fc is missing, the default value
%       will be 0.3 Hz.
%   ref is the gravitational field strength in the same units as Aw. The
%       default value is 9.81 which assumes that A is in m/s^2. Use ref=1 
%       if the units of Aw is g.
%   depth is a vector of the depth of the whale in meters.
%   speed is a the speed of the whale in m/s. This can either be a single
%       number, representing a constant speed over time, or a vector of
%       numbers, representing changing speeds over time.
%   plot is a conditional statement in which 'true' returns a plot of the
%       Awx column-vector with the peaks of detected lunge marked with a 
%       black dot and 'false' does not return a plot. Default is 'true'.
% 
% Outputs:
%   lunges is a vector of cues to lunges from Aw.
%   lunge_times is a vector of times in seconds since start of recording.


if nargin == 0
    help detect_lunges
end

if nargin < 2
    error("input for fs is required")
end

if nargin < 3 || isempty(fc)
    fc = .3;
end

if nargin < 3 || isempty(ref)
    ref = 9.81;
end

if nargin < 4
    error("input for depth is required")
end

if nargin < 5
    error("input for speed is required")
end

if nargin < 6 || or isempty(plot)
    plot = 'true';
end

%calculate pitch and roll
[pitch, roll , v] = a2pr(Aw, fc);

%calculate EXA (excess x-axes acceleration)
if ref == 0 || isempty(ref)
    g = 9.81;
elseif ref == 1
    g = 1;
end
EXA = Aw(:,1) - g * sin(pitch);

%calculate njerk of EXA signal and Aw signal
j_EXA = njerk(EXA, fs);
j_Aw = njerk(Aw, fs);

%low-pass filter all signals to remove noise from signal
Aw_filt = fir_nodelay(Aw, round(8 / fc), fc, 'low');
EXA_filt = fir_nodelay(EXA, round(8 / fc), fc, 'low');
j_EXA_filt = fir_nodelay(j_EXA, round(8 / fc), fc, 'low');
j_Aw_filt = fir_nodelay(j_Aw, round(8 / fc), fc, 'low');
speed_filt = fir_nodelay(speed, round(8 / fc), fc, 'low');

%chunk signals into one second bins
for a = 1:round(size(Aw_filt, 1) / fs)
    A_bin = cell(round(size(Aw_filt, 1) / fs), 1);
    bin = 1:fs:(size(Aw_filt, 1));
    b = 1;
    while b < length(bin)
        %for b = 1:length(bin)
        %while  length(bin)
        %    break
        %end
        A_bin{a} = Aw_filt((bin(b):bin(b+1)),1);
        b = b + 1;
    end
end

%chunk signals into 1 second bins
%separate bins into shallow vs deep depths
%find bins with peak values
%compare different signals looking for matching parameters characteristic
%   of a lunge
%mark these lunge detected bins and then plot the signal with the bins
%   marked

end