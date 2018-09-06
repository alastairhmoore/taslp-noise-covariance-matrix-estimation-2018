function[cm] = get_voicebox_thermliny_linecolors(N)
%GET_VOICEBOX_THERMLINY_LINECOLORS is a handy wrapper to get nicely scaled
%colors for plotting lines
cm = v_colormap('v_thermliny',[],N+1);
cm = flipud(cm(1:end-1,:));