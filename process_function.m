%Implementation of Nature protocol
%Hongbo Jia, Nathalie L Rochefort1, Xiaowei Chen & Arthur Konnerth1 "In
%vivo two-photon imaging of sensory-evoked dendritic calcium signals in cortical neurons"
%
%Implementation copyright Petros Xanthopoulos 2013-2014
%usage: signalout=process_function(signalin,t_0,t_1,t_2)
% where
% input: signalin is the raw signal 
%t_0,t_1,t_2 are the parameters described in Nature protocol paper
%comments: for a 30Hz imaging systems the following parameter setup is
%recommended (empirical note on Nature paper): 
%t_0= 0.2;
%t_1=0.75;
%t_2=3;


function signalout=process_function(signalin,t_0,t_1,t_2)

F_0=[];

Fs=30; %sampling frequency

t_0_s=floor(t_0*Fs);
t_1_s=floor(t_1*Fs);
t_2_s=floor(t_2*Fs);

F_sm = smooth(signalin,t_1_s);

for i=(t_2_s+1):length(signalin)
    F_0=[F_0 min(F_sm(i-t_2_s:i))];
end

R_0=(F_sm((t_2_s+1):end)-F_0')./F_0';

R_0_sm = EWMA(R_0,t_0_s);


signalout=R_0_sm;
%plot((1:length(R_0_sm))/30,RawIntDen1((t_2_s+1):end))
%hold on
%plot((1:length(R_0_sm))/30,R_0_sm,'r')


