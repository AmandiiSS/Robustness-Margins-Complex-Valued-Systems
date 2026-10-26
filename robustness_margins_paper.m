% Gain, phase, delay margins and sensitivity peak for complex systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num and den are the complex numerator and denominator polynomials
% (set in the workspace):
% num = N_ol(s)
% den = D_ol(s)
% Where N(s)/D(s) is the open-loop transfer function of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREVIOUS CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the coef of the N(jw) and D(jw)
numjw=num; num_length=length(num);
for k=1:1:num_length; numjw(k)=numjw(k)*1i^(num_length-k); end
numjw_conj=conj(numjw);
denjw=den; den_length=length(den);
for k=1:1:den_length; denjw(k)=denjw(k)*1i^(den_length-k); end
denjw_conj=conj(denjw);
%
epsx=1e-12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE MARGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux_num=conv(numjw,numjw_conj); aux_den=conv(denjw,denjw_conj);
aux_num_ext=[zeros(1,length(aux_den)-length(aux_num)) aux_num];
w_roots_pm=roots(aux_num_ext-aux_den);
w_i_pm=[]; pm_candidates=[];
for k = 1:1:length(w_roots_pm)
    if(abs(imag(w_roots_pm(k)))<epsx)
        w_i_pm = [w_i_pm, w_roots_pm(k), w_roots_pm(k)];
        phase = rad2deg(angle(polyval(numjw,w_roots_pm(k))/polyval(denjw,w_roots_pm(k))));
        if(phase >= 0)
            new_phase=phase-180; new_phase_opposite=360+new_phase;
        else
            new_phase=180 + phase; new_phase_opposite=new_phase-360;
        end
        pm_candidates=[pm_candidates, new_phase, new_phase_opposite];
    end
end
%
pm_pos=Inf; pm_neg=Inf;
for k = 1:1:length(pm_candidates)
    if(pm_candidates(k) > 0 && pm_candidates(k) < pm_pos)
        pm_pos=pm_candidates(k);
    end
    if(pm_candidates(k) < 0 && abs(pm_candidates(k)) < abs(pm_neg))
        pm_neg=pm_candidates(k);
    end
end
pm=[pm_neg pm_pos]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELAY MARGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dm_candidates = deg2rad(pm_candidates)./w_i_pm;
dm_pos=Inf; dm_neg=Inf;
for k = 1:1:length(dm_candidates)
    if(dm_candidates(k)>0 && dm_candidates(k)<dm_pos) dm_pos=dm_candidates(k); end
    if(dm_candidates(k)<0 && abs(dm_candidates(k))<abs(dm_neg)) dm_neg=dm_candidates(k); end
end
%
if(dm_pos == Inf) dm_pos=0; end
if(dm_neg == Inf) dm_neg=0; end
dm=[dm_neg*1000 dm_pos*1000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIN MARGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_roots_gm=roots(imag(conv(numjw,denjw_conj)));
w_i_gm=[]; gm_candidates=[];
for k = 1:1:length(w_roots_gm)
    if(abs(imag(w_roots_gm(k)))<epsx && abs(w_roots_gm(k))>epsx && real(polyval(numjw, w_roots_gm(k))/polyval(denjw, w_roots_gm(k)))<-epsx)
        w_i_gm=[w_i_gm, w_roots_gm(k)];
        gain=abs(polyval(numjw, w_roots_gm(k))/polyval(denjw, w_roots_gm(k)));
        new_gain_margin=1/gain;
        gm_candidates=[gm_candidates, new_gain_margin];
    end
end
%
gm_pos=Inf; gm_neg=0;
for k = 1:1:length(gm_candidates)
    if(gm_candidates(k)<gm_pos && gm_candidates(k)>1) gm_pos=gm_candidates(k); end
    if(gm_candidates(k)>gm_neg && gm_candidates(k)<1) gm_neg=gm_candidates(k); end
end
gm = [gm_neg gm_pos]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENSITIVITY PEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum sensitivity for complex systems
%d=denjw; n=numjw; % Polynomials in omega
%ds=denjw_conj; ns=numjw_conj;
dp=zeros(1,den_length-1); np=zeros(1,num_length-1);
for k=1:den_length-1; dp(k)=denjw(k)*(den_length-k); end
for k=1:num_length-1; np(k)=numjw(k)*(num_length-k); end
%
% Maximum sensitivity
%poly=conv(conv(denjw_conj+[zeros(1,length(denjw_conj)-length(numjw_conj)),numjw_conj],denjw_conj),[0, conv(dp,numjw)]-conv(denjw,np));
poly=conv(conv(denjw_conj+[zeros(1,length(denjw_conj)-length(numjw_conj)),numjw_conj],denjw_conj), conv(dp,numjw)-conv(denjw,np));
oms=roots(poly+conj(poly)); noms=length(oms);
%
ms=0;
for k=1:noms
    if abs(imag(oms(k)))<epsx
        om=real(oms(k));
        fresp=polyval(num,1i*om)/polyval(den,1i*om);
        som=abs(1/(1+fresp));
        if som>ms; ms=som; end
    end
end
ms