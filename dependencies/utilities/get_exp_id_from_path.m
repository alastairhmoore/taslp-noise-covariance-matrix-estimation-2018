function[out_str] = get_exp_id_from_path(path_to_working_dir,timestamp)
%get_exp_id_from_path returns a string used to uniquely identify an experiment
%in the format
% <create_date>_<user>_<id>_<run_date>_<run_time>_<tag>
%where
% <create_date>_<user>_<id>:
%   uniquely defines the set of scenarios which will be generated
% <create_date>:
%   yyyymmdd
% <user>:
%   initials of person setting up the experiment.
% <id>:
%   2 digit integer - distinguish between multiple tests setup by a single
%   user on the same day
% <run_time>:
%   yyyymmdd_HHMMSS indicates the date and time at which the code was run
%   - see datestr for fields
% ssssssss:
%   short description of experiment for easy identification
%

if nargin<2
    timestamp=now;
end

[path_to_parent_dir,this_dir,~] = fileparts(path_to_working_dir);
[~,parent_dir,~] = fileparts(path_to_parent_dir);

idc_split = strfind(this_dir,'_');
if numel(idc_split)<2
    error('parent folder name should be in format yyyymmdd_user_id_desc')
else
    s1_date = this_dir(1:idc_split(1)-1);
    s2_user = this_dir(idc_split(1)+1:idc_split(2)-1)
    if numel(idc_split)<3
        s3_id = this_dir(idc_split(2)+1:end);
        s4_desc = '';
    else
        s3_id = this_dir(idc_split(2)+1:idc_split(3)-1);
        s4_desc = this_dir(idc_split(3)+1:end)
    end
end


%light touch validation - convert strings to numbers
if ~isdatestr(s1_date) || ~isinteger(s3_id)
    error('Directory name is badly formed')
end

out_str = sprintf('%s_%s_%s_%s_%s',...
    s1_date, s2_user, s3_id, datestr(timestamp,'yyyymmdd_HHMMSS'),s4_desc);

function[bool] = isdatestr(in)
bool = 1;
if isempty(str2double(in)) %check its a number
    bool = 0;
end
if length( in) ~=8 %with correct number of digits
    bool = 0;
end
% hard code a sensible date range
mindate = datenum('20170102','yyyymmdd');
maxdate = datenum('20201231','yyyymmdd');
thisdate =  datenum(in,'yyyymmdd');
if thisdate > maxdate || thisdate < mindate
    bool=0;
end

function[bool] = isinteger(in)
bool = 1;
try
    in = str2double(in);
catch
    bool = 0;
end
if mod(in,1)~=0
    bool = 0;
    
end
if in>99
    bool = 0;
end