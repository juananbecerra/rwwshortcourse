function figureName(name)
% Function that makes the current figure the one that has name "name"

fh = findobj( 'Type', 'Figure', 'Name', name);
if length(fh)==0
    figure('Name',name);
else
    figure(fh(1));
end