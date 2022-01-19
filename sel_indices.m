function indices = sel_indices(xdpd,ydpd,perc)

% Seleccion del trozo de señal para modelar
%perc =  0.015;
[maxy, indy] = max(abs(ydpd));
indmodinf = floor(length(ydpd)*perc*floor((indy/length(ydpd))/perc))+1;
indmodsup = ceil(length(ydpd)*perc*ceil((indy/length(ydpd))/perc))-1;
if(indmodsup>length(ydpd))
    indmodsup = length(ydpd);
    indmodinf = length(ydpd)-floor(perc*length(ydpd));
end

indices = indmodinf:indmodsup;