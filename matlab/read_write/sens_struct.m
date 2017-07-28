function    X = sens_struct(data,fs,depid,type,name)

%    X = sens_struct(data,fs,depid,type)  % regularly sampled data
%    or
%    X = sens_struct(data,T,depid,type)   % irregularly sampled data
%    or
%    X = sens_struct(data,fs,depid,type,name)  % regularly sampled data
%    or
%    X = sens_struct(data,T,depid,type,name)   % irregularly sampled data
%
%    Generate a sensor structure from a sensor data vector or matrix.
%
%    Inputs:
%    data is a vector or matrix of sensor data. This can be in any unit
%     and frame but the metadata generated by sens_struct assumes a default
%     unit and frame. Change these manually in X after running sens_struct
%     to the correct unit and frame.
%    fs is the sampling rate of the sensor data in Hz (samples per second).
%    T is the time in seconds of each measurement in data for irregularly
%     sampled data. The time reference (i.e., the 0 time) should be with
%     respect to the start time of the deployment.
%    depid is a string containing the deployment identifier for these data.
%    type is a string containing the first few letters of the sensor type,
%     e.g., acc for acceleration. These will be matched to the list of
%     sensor names in the sensor_names.csv file. If more than one sensor 
%     matches type, a warning will be given. type can be in upper or lower 
%     case.
%    name is the optional name to give the variable, e.g., T_EXT. If a name
%     is not given, a default value will be selected.
%
%    Returns:
%    X is a sensor structure with metadata fields pre-populated from the 
%     sensor_names.csv file. Change these as needed to the correct values.
%
%    Example:
%     A = sens_struct(Aw,fs,'md13_134a','acc')
%     A.frame = 'animal' ;     % change frame indication
%
%    Valid: Matlab, Octave
%    markjohnson@st-andrews.ac.uk
%    Last modified: 10 July 2017

if nargin<4,
   help sens_struct
   [S,hdr]=read_csv('sensor_names.csv') ;
   fprintf(' Predefined sensor types:\n') ;
   for k=1:size(S,1),
      fprintf('  %s\n',S(k).name);
   end
   X = [] ;
   return
end

if length(fs)==1,    % regularly sampled data
   X.data = data ;
	X.sampling = 'regular' ;
	X.sampling_rate = fs ;
	X.sampling_rate_unit = 'Hz' ;

else                 % irregular data
   if length(fs) ~= size(data,1),
      fprintf('Error: number of sampling times does not match number of samples\n') ;
      return
   end
   X.data = [fs(:),data] ;
	X.sampling = 'irregular' ;
	X.sampling_time = 'column 1' ;
	X.sampling_time_unit = 'second' ;
   fs = [] ;
end

X.depid = depid ;
X.creation_date = datestr(now) ;
X.history = 'sens_struct' ;

% read in sensor names database and compare against type
[S,hdr]=read_csv('sensor_names.csv',1) ;
k = strncmpi(type,{S.name},length(type)) ;
if all(k==0),
   fprintf('Warning: unknown sensor type %s. Set metadata manually\n', type) ;
   X.name = type ;
   return ;
end

if sum(k)>1,
   fprintf('More than one sensor type match "%s". Retry with a longer type name.\n',type) ;
   X= [] ;
   return
end

k = find(k) ;
nc = str2double(S(k).axes) ;
if size(data,2)~=nc,
   fprintf('Warning: size of data does not match number of columns (%d) expected for %s\n',nc,S(k).name) ;
end

if nargin<5,
   X.name = S(k).abbrev ;
else
   X.name = name ;
end

X.full_name = S(k).name ;
X.description = S(k).description ;
X.unit = S(k).def_units ;
X.unit_name = S(k).def_unit_name ;
X.unit_label = S(k).def_label ;
X.start_offset = 0 ;
X.start_offset_units = 'second' ;

if ~isempty(S(k).def_cols),
   if isempty(fs),
      X.column_name = strcat('time,',strip_quotes(S(k).def_cols)) ;
   else
      X.column_name = strip_quotes(S(k).def_cols) ;
   end
end

if ~isempty(S(k).def_frame),
   X.frame = S(k).def_frame ;
end

if ~isempty(S(k).def_axes),
   X.axes = S(k).def_axes ;
end
