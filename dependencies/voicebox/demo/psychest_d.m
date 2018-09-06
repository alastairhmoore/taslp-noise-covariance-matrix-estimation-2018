%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   this script uses psycest to estimate the parameters of a simulated psychometric function  %
%                                                                                             %
%   It plots graphs after the trial and then asks the user to enter the number of trials to   %
%   skip before plotting graphs again.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
usemat = ''; % blank string [normal] or else file name of previous probe sequence
LOG=fopen('psycest.log','wt'); % open log file
q0=[];              % clear parameter structure with algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   optional algorithm parameters   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q0.la=1;            % look 2-ahead flag
% q0.cs=1             % weighting of slope relative to SRT
% q0.op=0.1;          % outlier probability
% q0.qm=1;            % model: 1=logistic, 2=gaussian
% q0.cf=2;            % cost function type: 1=variance, 2=entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     define ground truth model   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prob @ threshold, threshhold (dB), slope at threshold (prob/dB), p(miss), p(guess), function type]
qq0 = [0.5 0 0.1 0.01 0 1]'; % ground truth model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     define model initialization sent to psycest  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [thresh, p(miss), p(guess), SNR min, SNR max, Slope min, Slope max (prob/dB)]
p0 = [qq0(1) 0.01 0 -20 20 0 0.5]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     initialize the model          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(usemat)
    % initialize a model
    [xx,ii,m,v]=psycest(-1,p0,q0,[],LOG); % enable logging
    nt=200; % number of trials
else
    load(usemat); % restore previous values of qq0, p9, q9, msr9
    [xx,ii,m,v]=psycest(-1,p9,q9,[],LOG);
    res0=msr9;
    nt=size(res0,1); % number of trials
end
srtint=qq0(2)+[-qq0(1) 0 1-qq0(1)]/qq0(3); % SRT and slope intercepts
srtax=linspace(srtint*[1.6; 0; -0.6],srtint*[-0.6; 0; 1.6],200); % SRT axis for ground truth

psycest(0); % output model info to log file
[p1,q1]=psycest(0);
savs=zeros(nt,4);
ihalt=1;
inchalt0=1;
i=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     loop for each probe           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while i<nt
    i=i+1;
    if isempty(usemat)
        rr=psychofunc('r',qq0,xx);
    else
        xx=res0(i,2);
        rr=res0(i,3);
    end
    savs(i,1:2)=[xx,rr];
    [xx,ii,m,v]=psycest(1,xx,rr);
    savs(i,3:4)=m(:,1,1)';
    if i==ihalt
        [xx,ii,m,v,mr,vr]=psycest(1) % Determine robust outputs in addition [takes longer]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    print errors on console            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'    SRT: True=%.2f, Est=%.2f(+-%.2f), Err=%.2f\n    Slope:  True=%.3f, Est=%.3f(+-%.3f), Err=%.3f\n', ...
            qq0(2),m(1),sqrt(v(1)),m(1)-qq0(2),qq0(3),m(2),sqrt(v(3)),m(2)-qq0(3));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    plot the ground truth psychometric function      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(100);
        plot(srtax,psychofunc('',qq0,srtax),'--k'); % plot the ground truth graph
        hold on 
        plot(srtax,psychofunc('',[p1(1) mr(1:2) p1(2:3)' q1(10)]',srtax),'-r'); % plot the robust mean graph
        plot(srtax,psychofunc('',[p1(1) m(1:2) p1(2:3)' q1(10)]',srtax),'-b'); % plot the mean graph     
        plot(srtint([1 3 3]),[0 1 0],'--g',[srtax(1) qq0(2) qq0(2)],[qq0(1) qq0(1) 0],'--g');  
        hold off
        axisenlarge([-1 1]);
        legend('True','Robust','Mean','Location','NorthWest');
        title('Ground Truth Psychometric Function');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    plot psycest pdf                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(101);
        psycest(1,'p');
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        hold on
        plot([qq0([2 2]) xlim'],[ylim' log(qq0([3 3]))],'-g');
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    plot next probe cost               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(102);
        psycest(1,'x');
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        hold on
        if q1(18)  % if look 2-ahead
            plot(xlim,[1;1]*srtint,'-g',[1;1]*srtint,ylim,'-g');
            title('Expected cost after next 2 probes [green=true SRT+intercepts]');
        else
            plot([1;1]*srtint,ylim,'-g');
            title('Expected cost after next probe [green=true SRT+intercepts]');
        end
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    plot probe history                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(103);
        psycest(1,'h');
        hold on
        plot([1 i],[1;1]*srtint,':r')
        hold off
         axisenlarge([-1 -1.03]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    plot cost evolution                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(104);
        psycest(1,'c');
        axisenlarge([-1 -1.03]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    make all the figures visible and request user input     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tilefigs;
        if i<nt
            inchalt=input(sprintf('Enter number of trials to skip or 0=quit, -1=Inf [%d]:',inchalt0));
            if isempty(inchalt)
                ihalt=min(i+inchalt0,nt);
            elseif inchalt==0
                nt=i;
            elseif inchalt<0
                ihalt=nt;
            else
                ihalt=min(i+inchalt,nt);
                inchalt0=inchalt;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    save the probe sequnce to allow repeats   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p9,q9,msr9]=psycest(0); % save model parameters and results
save('pqmsr.mat', 'qq0', 'p9','q9','msr9');
fclose(LOG);
tilefigs;
