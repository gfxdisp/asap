function [] = ahelp(c)
%AHELP Open browser and show documentation or examples.

if nargin==0
  c = 0;
end

s = which('ainit');
s = strrep(s,'ainit.m','');
s = fullfile(s,'html',filesep);

if c>0
  web([s 'example' num2str(c) '.html'])
elseif c==0
  web([s 'documentation.html'])
elseif c==-1
  web([s 'examples.html'])
end