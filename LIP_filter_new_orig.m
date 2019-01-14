function rq = LIP_filter_new_orig(data)%, method='nelder-mead', xtol=1E-8):
%data is a standard dp struct

nEvents = length(data.x_corrected);
nPulses = size(data.x_corrected , 1);
rq  = init_lip_rqs (nEvents, nPulses);
data = find_errors(data);
for i = 1: nEvents
   % i
    rq = fillS1RQs(i, data, rq);
    rq.eventWithinAcq(i) = data.event_number(i);
    rq.file_number(i) = data.file_number(i); % this is the file in suffix of the acq
    rq.luxstamp(i) = data.luxstamp_samples(i); % save the luxstamp 

    if rq.nS1(i) > 0 % if there is no s1 we don't care about the event (and no corrected pulses)
        rq = fillCorrS2PARQs(i, data, rq); %
        %maybe have condition about not meeting filter?
        rq = fillCorrS2AngleRQs(i, data, rq);
        rq = fillCorrS2FitRQs(i, data, rq); %, method=method, xtol=xtol)
        rq.CorrS1_keV(i) = 13.7/1000 * 1/0.117 .* rq.largestS1PA(i);  
        rq.CorrS2_keV(i) = 13.7/1000 * 1/12.1 .* rq.sumCorrS2PA_xyz_corr(i);
        rq.total_Corr_keV(i) = rq.CorrS1_keV(i) + rq.CorrS2_keV(i);
    end
    
    pulse_filter = data.pulse_classification(:,i) == 2 | data.pulse_classification(:,i) == 4;% | data.s2_like_class5(: , i) == 1 ;
    
    rq.raw_S2_total_keV(i) = sum(data.pulse_area_phe( pulse_filter , i));
end
rq.raw_total_keV= 13.7/1000 * 1/0.117 .* rq.largestS1PA; 
rq.raw_total_keV = rq.raw_total_keV + rq.raw_S2_total_keV; 

end

function data = find_errors(data)

z_err =ones(1, length(data.z_drift_samples));
data.z_err = z_err;

sys = 0.1; %1mm in Run03
ri = data.sd_radius_inf;
rs = data.sd_radius_sup;
sys_err = ones(1, length(ri)) *sys;
sigma_r = sqrt( (0.5*(rs + ri)./sqrt(2.3)).^2. + sys_err.^2.);
data.sigma_r = sigma_r;
sd_phiXR = data.sd_phiXR;
sigma_tot = sqrt(sigma_r.^2. + sd_phiXR.^2. + z_err.^2.);
data.sigma_tot = sigma_tot ;
end

function rq = fillCorrS2FitRQs(i, data, rq) %, xtol=1E-8, method='nelder-mead'):
   % """ Do fits and fill RQs from fitting, return 1 if the event passes the filter, 
    %    return 0 fit the event doesn't pass the filter """

    pulsetypes = data.pulse_classification(: , i);
    x_corr_cm = data.x_corrected(: , i);
    y_corr_cm = data.y_corrected(: , i);
    sigma_tot = data.sigma_tot(: , i);
    pulse_area_phe = data.pulse_area_phe(:,i);

    drift_us = 10. * 0.001 * double(data.z_drift_samples(: , i));
    drift_to_cm = 1.51/10.;% #cm/us
    z_cm = 50. - drift_us * drift_to_cm  ;
    pulse_length = data.pulse_length_samples(:, i );
    %z_cm = data.z_cm(: , i) #50. - drift_us * drift_to_cm  
    pulse_filter = (pulsetypes == 2 | pulsetypes == 4 )  & (x_corr_cm > -99.);
    rq.nCorrS2PulsesBeforeCut(i) = sum(pulse_filter);
    pulse_len_per_samp = double(pulse_area_phe) ./ double(pulse_length);
    plps_cut = pulse_len_per_samp > 1;
    pulse_filter = pulse_filter .* plps_cut; 
    pulse_filter = pulse_filter > 0; % reconvert to logical
    
    if sum(pulse_filter) > 4 % was >2

        xarr = x_corr_cm(pulse_filter);
        yarr = y_corr_cm(pulse_filter);
        zarr = z_cm(pulse_filter);
        %#Fill RQs about the points #############
        %###################################################
        [std_d, avg_d] = getPointsDistanceInfo(xarr,yarr,zarr);
        rq.stdPointDistS2Corr(i) = std_d;
        rq.avgPointDistS2Corr(i) = avg_d;
        
        data_vec = [xarr(end)-xarr(1), yarr(end)-yarr(1), zarr(end)-zarr(1)];
        rq.CorrS2Ldata(i) = norm(data_vec);
        data_dir = data_vec ./ norm(data_vec);

      %   # Fill RQs From track fitting points #############
     %   ########################################################
        tot_err = sigma_tot(pulse_filter); %#np.sqrt(sigma_r[pulse_filter]**2. + sd_phiXR[pulse_filter]**2. + z_err[pulse_filter]**2.)
        
        svd_params = SVD_params(xarr,yarr,zarr);
        svd_params = double(svd_params);
        for idx=1:length(svd_params) 
            rq.CorrS2SVD_params(idx ,i)=svd_params(idx);
        end
             

       % #fit with the magnitude-wise + angle*distance penalty
        %
        %{
        res = scipy.optimize.minimize(chi2_magwise_angle_penalty, svd_params, args=(xarr,yarr,zarr,tot_err,False), method=method, options={'xtol': xtol, 'disp':False})
        if res.success:
            params_opt = res.x;         
            %#get the chi2 of the fit result
            [chi2, chi2_dof] = chi2_magwise_angle_penalty(params_opt, xarr,yarr, zarr,tot_err,True);
            rq.CorrS2Chi2_anglePen(i) = chi2;
            rq.CorrS2Chi2_anglePen_dof(i) = chi2_dof;
            for idx=1:length(params_opt)
                rq.CorrS2TrackFit_anglePen_params(idx , i) = params_opt(idx);
            end
                
            track_dir = params_opt(4:6);
            track_dir = track_dir./norm(track_dir); 
            rq.CorrS2TrackTheta_anglePen(i) = acos(abs(track_dir(3)))*180./pi;
            rq.CorrS2TrackPhi_anglePen(i) = atan2(track_dir(2), track_dir(1))*180./pi;
            rq.CorrS2DataToFitAng_anglePen(i) = acos(dot(data_dir, track_dir)) * 180. / pi;    
            
            [L, PE, PL] = get_tracklength(params_opt);
            if L > -1
                rq.CorrS2Lfit_anglePen(i) = L;
                for idx= 1:length(PE)
                    rq.CorrS2TrackEnter_anglePen(idx , i) = PE(idx);
                end
                for idx=1:length(PL)
                    rq.CorrS2TrackExit_anglePen(idx , i) = PL(idx);
                end
            end
        end
          %}      
        %#fit with the magnitude-wise chi2
        %***************
       % res = scipy.optimize.minimize(chi2_magwise, svd_params, args=(xarr,yarr,zarr,tot_err,False), method=method, options={'xtol': xtol, 'disp':False})   
        params_opt=fminsearch(@chi2_magwise,svd_params);%,xarr,yarr,zarr,xerr,yerr,zerr);%, options);
        [Chi2] = chi2_magwise(params_opt);%, xarr,yarr,zarr, xerr, yerr, zerr);
        
        dof_3d_line = 4 ;%#2 for direction on the unitsphere, 3 for the point - 1 b/c you can translate the point
        chi2_dof = Chi2 / (length(xarr) - dof_3d_line);  
        
        if 1
           % params_opt = res.x;
            %#get the chi2 of the fit result
           % [chi2, chi2_dof] = chi2_magwise(params_opt, xarr,yarr, zarr,data,True);
            rq.CorrS2Chi2(i) = Chi2;
            rq.CorrS2Chi2_dof(i) = chi2_dof;
            for idx = 1:length(params_opt)
                rq.CorrS2TrackFit_params(idx , i) = params_opt(idx);
            end
            
            track_dir = params_opt(4:6); %want the slopes of the parameterized lines
            track_dir = track_dir./norm(track_dir); 
            rq.CorrS2TrackTheta(i) = acos(abs(track_dir(3)))*180./pi;
            rq.CorrS2TrackPhi(i) = atan2(track_dir(2), track_dir(1))*180./pi;
            rq.CorrS2DataToFitAng(i) = acos(dot(data_dir, track_dir)) * 180. / pi;   
            
            [L, PE, PL] = get_tracklength(params_opt);
            if L > -1
                rq.CorrS2Lfit(i) = L;
                for idx =1:length(PE)
                    rq.CorrS2TrackEnter(idx , i) = PE(idx);
                end
                for idx = 1:length(PL)
                    rq.CorrS2TrackExit(idx , i) = PL(idx);
                end
            end
        %return 1; %need some var = 'true'
     %   else
                
     %   return 0
        end
    end
    
    function [chi2_event] = chi2_magwise(params)% , xarr,yarr,zarr,xerr,yerr,zerr)
%""" Return the magnitude-wise claculated component-wise chi2 and chi2/dof of events for minimization. 
% Chi2 procedure does not define any compoent (x,y,z) to be dependent on
%        any other component. """
%      #params are the input line parameters in [x0,y0,z0,a,b,c] format
%    #xarr,yarr,zarr are the x,y,z poisitons of all the hits
%    #xerr,yerr,zerr are their errors - we need to have them spanning
%    multiple functions for the minimization to work (because here matlab is dumb).
 chi2_event = 0;

for j = 1:length(xarr)
        x_val = xarr(j);
        y_val = yarr(j);
        z_val = zarr(j);
    
        d = distance(params, x_val,y_val,z_val);
        
       % #temp = d/err
        chi2 = norm(d).^2 / tot_err(j).^2;
        chi2_event = chi2_event + chi2;
        

        %dof_3d_line = 4 ;%#2 for direction on the unitsphere, 3 for the point - 1 b/c you can translate the point
        %#is there a way to do this fit for events with only 3 hits? 4hits?
      %  chi2_dof = chi2_event / (length(xarr) - dof_3d_line); 

end
end



function d = distance(params,x,y,z)
%    """ Return the orthogonal distance vector from the point (x,y,z)
%        to the line described by params """
    
%    # the line follows vector equation: x = x0 + v*t 
    x_init = [params(1), params(2), params(3)];
    v = [params(4), params(5) ,params(6)];
    v = v ./ norm(v); % #make v a unit vector
    
    P = [x,y,z]; % #data point P
   
   % #get a vector pointing to the closest point 
   %#on the line from the point P
   % #first get r_vec = a vector from P to a point on
   % #the line
    r_vec = (P-x_init) ;
  %  #then get r_proj, the projection of r_vec along the line
    r_proj = dot(r_vec, v);
   % #proceed along the line in the i
   % #direction of the line, v
   % #for magnitude r_proj
   % #this describes the point on the 
   % #line r_fin that is closest to P
    r_fin = x_init + v.* r_proj;

    % #d is the vector between P and the point i
    % #on the line closest to it, r_fin
    d = P - r_fin;
  
end

end


function rq = fillCorrS2PARQs(i, data, rq)
   % """ Fill Pulse Area RQs for Corrected S2 pulses."""
   % #Fill Corr Angle Info - Require pulses to be good corrected values for S2 only

    pulsetypes = data.pulse_classification(: , i);
    x_corr_cm = data.x_corrected(: , i);

    pulse_area_phe = data.pulse_area_phe(: , i);
    xyz_corr_pulse_area_phe = data.xyz_corrected_pulse_area_all_phe(: , i);
    xyz_corr_bot_pulse_area_phe = data.xyz_corrected_pulse_area_bot_phe(: , i);

    pulse_filter = (pulsetypes == 2 | pulsetypes == 4  | s2_like_class5 ==1)  & (x_corr_cm > -99.);
    
    pulse_length = data.pulse_length_samples(:, i );
    pulse_len_per_samp = double(pulse_area_phe) ./ double(pulse_length);
    plps_cut = pulse_len_per_samp > 1;
    pulse_filter = pulse_filter .* plps_cut; %PAT added 0818
    pulse_filter = pulse_filter > 0; % reconvert to logical
    if sum(pulse_filter) > 4 % was > 2. now want 5 pulses %#can calculate DeltaTheta & DeltaPhi for n=3 events
        data.luxstamp_samples(i)
        rq.nCorrS2FitPoints(i) = sum(pulse_filter);
        
       % ## Fill RQs about Pulse Area ###################
       % ###################################################
        
       %uncorrected energies
        mypulseareas = pulse_area_phe(pulse_filter);
        rq.stdCorrS2PA(i) =  std(mypulseareas);    
        [mx, ~] = getMaxAndMinZScore(mypulseareas);
        rq.maxZScoreCorrS2PA(i) = mx       ;

        cs =  cumsum(mypulseareas);
        for idx=1:length(cs)
            rq.cumSumCorrS2PA(idx , i) = cs(idx);
        end
                   
        rq.sumCorrS2PA(i) =  sum(mypulseareas);
        rq.firstCorrS2PA(i) = mypulseareas(1);
        rq.lastCorrS2PA(i)  = mypulseareas(end);
        rq.maxCorrS2PA(i) =  max(mypulseareas);
        rq.minCorrS2PA(i)  =  min(mypulseareas); 
        rq.first2CorrS2PA(i) =  sum(mypulseareas(1 : 2));
        rq.last2CorrS2PA(i)  =  sum(mypulseareas(end-1 : end));
       % ########
       % all corrected energy
        mypulseareas = xyz_corr_pulse_area_phe(pulse_filter);
        rq.moyalTest_xyzCorrPA(i) = moyal_test(mypulseareas);
        rq.stdCorrS2PA_xyz_corr(i) =  std(mypulseareas);        
        [mx, ~] = getMaxAndMinZScore(mypulseareas);
        rq.maxZScoreCorrS2PA_xyz_corr(i) = mx;  
         
        cs =  cumsum(mypulseareas);
        for idx=1:length(cs)
            rq.cumSumCorrS2PA_xyz_corr(idx , i) = cs(idx);
        end
        
        rq.sumCorrS2PA_xyz_corr(i) =  sum(mypulseareas);
        rq.firstCorrS2PA_xyz_corr(i) = mypulseareas(1);
        rq.lastCorrS2PA_xyz_corr(i)  = mypulseareas(end);
        rq.maxCorrS2PA_xyz_corr(i) =  max(mypulseareas);
        rq.minCorrS2PA_xyz_corr(i)  =  min(mypulseareas);
        rq.first2CorrS2PA_xyz_corr(i) =  sum(mypulseareas(1:2));
        rq.last2CorrS2PA_xyz_corr(i)  =  sum(mypulseareas(end-1:end));
        %########
        % corrected bottom energy
        mypulseareas = xyz_corr_bot_pulse_area_phe(pulse_filter);            
            
        rq.stdCorrS2PA_xyz_corr_bot(i) =  std(mypulseareas);        
        [mx, ~] = getMaxAndMinZScore(mypulseareas);
        rq.maxZScoreCorrS2PA_xyz_corr_bot(i) = mx ;

        cs =  cumsum(mypulseareas);
        for idx=1:length(cs)
            rq.cumSumCorrS2PA_xyz_corr_bot(idx,i) = cs(idx);
        end
        
        rq.sumCorrS2PA_xyz_corr_bot(i) =  sum(mypulseareas);
        rq.firstCorrS2PA_xyz_corr_bot(i) = mypulseareas(1);
        rq.lastCorrS2PA_xyz_corr_bot(i)  = mypulseareas(end);
        rq.maxCorrS2PA_xyz_corr_bot(i) =  max(mypulseareas);
        rq.minCorrS2PA_xyz_corr_bot(i)  =  min(mypulseareas);
        rq.first2CorrS2PA_xyz_corr_bot(i) =  sum(mypulseareas(1:2));
        rq.last2CorrS2PA_xyz_corr_bot(i)  =  sum(mypulseareas(end-1:end));
    end
    
    sprintf('Done with track finding')
end  


function rq = fillCorrS2AngleRQs(i, data, rq)
   % """ Fill Pulse Area RQs about angles between poins."""
    pulsetypes = data.pulse_classification(: , i);
    x_corr_cm = data.x_corrected(: , i);
    y_corr_cm = data.y_corrected(: , i);
    pulse_area_phe = data.pulse_area_phe(: , i);

    drift_us = 10. * 0.001 * double(data.z_drift_samples(: , i)); %this line and the line below were commed out and data.z_cm was in
    drift_to_cm = 1.51/10.;% #cm/us
    z_cm = 50. - drift_us .* drift_to_cm  ;
   % z_cm = data.z_cm(: , i)% #50. - drift_us * drift_to_cm  
    pulse_filter = (pulsetypes == 2 | pulsetypes == 4  )  & (x_corr_cm > -99.);
    pulse_length = data.pulse_length_samples(:, i );
    pulse_len_per_samp = double(pulse_area_phe) ./ double(pulse_length);
    plps_cut = pulse_len_per_samp > 1;
    pulse_filter = pulse_filter .* plps_cut; %PAT added 0818
    pulse_filter = pulse_filter > 0; % reconvert to logical
    if sum(pulse_filter) > 4 % was >2
       % #Fill RQs about Angles between points #############
        %###################################################
        xarr = x_corr_cm(pulse_filter);
        yarr = y_corr_cm(pulse_filter);
        zarr = z_cm(pulse_filter);
        
        [dth, dph, ath, aph] = getPointsCombAngleInfo(xarr,yarr,zarr);

        rq.stdCombThetaS2Corr(i) = dth;   
        rq.stdCombPhiS2Corr(i) =  dph;
        rq.avgCombThetaS2Corr(i) = ath;
        rq.avgCombPhiS2Corr(i) = aph;
    
        [dth, dph, ath, aph] = getPointsPairAngleInfo(xarr,yarr,zarr);

        rq.stdPairThetaS2Corr(i) = dth;   
        rq.stdPairPhiS2Corr(i) =  dph;
        rq.avgPairThetaS2Corr(i) = ath;
        rq.avgPairPhiS2Corr(i) = aph;
    end
end


function [dth, dph, ath, aph] = getPointsPairAngleInfo(xarr,yarr,zarr)
    %""" return standard deviation and average of angles between each consecutive pair of points """
    event = [xarr,yarr,zarr] ;% #row is (x,y,z)_i
    phis = zeros(length(xarr-1), 1);
    thetas = zeros(length(xarr-1), 1);

    for idx = 1:(length(xarr)-1)
        %#print ind
        vec = event(idx+1, :)-event(idx, :);
        %#print 'vec = ', vec
        %#print 'vec = ', vec
        vec = vec ./norm(vec); %#unit vector in direction of propogation
        phi = atan2(vec(2),vec(1)) .* 180 ./ pi;
        theta = acos(vec(3)) .* 180 ./pi; %#arccos(z-component / length(=1))
        phis(idx)= phi;
        thetas(idx) = theta;
    end
    %phis = np.array(phis)
    %thetas = np.array(thetas)
    dth = std(thetas);
    dph = std(phis);
    ath = mean(thetas);
    aph = mean(phis);
   
    %return dth, dph, ath, aph
end

function [dth, dph, ath, aph] = getPointsCombAngleInfo(xarr,yarr,zarr)
    %""" return standard deviation and average of angles between all combinations of N choose 2 points """
    
    event = [xarr,yarr,zarr] ;% #row is (x,y,z)_i
    indices = combnk(1:length(xarr), 2);
    phis = zeros( length(indices), 1);
    thetas = zeros( length(indices), 1);
    
 
    for k = 1:length(indices)
        ind = indices(k, :);
        %#print ind
        vec = event(ind(2), :)-event(ind(1), :);
        %#print 'vec = ', vec
        vec = single(vec);
        vec = vec ./norm(vec); %#unit vector in direction of propogation
        phi = atan2(vec(2),vec(1)) .* 180 ./ pi;
        theta = acos(vec(3)) .* 180 ./pi; %#arccos(z-component / length(=1))
        phis(k)= phi;
        thetas(k) = theta;
    end
   % phis = np.array(phis)
    %thetas = np.array(thetas)
    dth = std(thetas);
    dph = std(phis);
    ath = mean(thetas);
    aph = mean(phis);
   
   % return dth, dph, ath, aph
end


function [std_d, avg_d] = getPointsDistanceInfo(xarr,yarr,zarr)
   % """ return standard deviation and average of distances between each consecutive pair of points """
    event = [xarr,yarr,zarr]' ;% #row is (x,y,z)_i
    distances = zeros(1, length(xarr -1));

    for idx = 1:(length(xarr)-1)
        vec = event(:,idx+1)-event(:,idx);
        d = norm(vec) ;%#length of vector in direction of propogation from point i to point i+1
        distances(idx) = d;
    end
        
    %distances = np.array(distances)

    std_d = std(distances);
    avg_d = mean(distances);
   
    %return std_d, avg_d
end



function  params = SVD_params(xx,yy,zz)
%    """ Get the Singular Value Decomposition of the Event. 
%        This returns a line that represents the average behavior of the event, 
%        like PCA """

    evt = [xx,yy,zz];

    datamean = mean(evt);
    x_init=datamean(1);
    y_init=datamean(2);
    z_init=datamean(3);

  %  # Do an SVD (pca) on the mean-centered data.
    v = pca(evt - datamean);
   
    x_slope=v(1,1);
    y_slope=v(2,1);
    z_slope=v(3,1);
    
    params = [x_init,y_init,z_init,x_slope,y_slope,z_slope];
end

function [L, PE, PL] = get_tracklength(params) 
    %this is all written by me to avoid convhull

  top_cm = 50; 
  bot_cm = 0; 
  rmax =double(24.5); % cm this is the max radius of lux assuming circular shap
        
   x0=params(1);
   y0=params(2);
   z0=params(3);
   a=params(4);
   b=params(5);
   c=params(6);
   
   syms t; 
   eqn = x0.^2 + y0.^2 + 2*t.*(x0.*a + y0.*b) + (a.^2 + b.^2).*t.^2 == rmax.^2;
   t = solve(eqn, t);
   t = double(t); % t values are the parameter positions where the particle hists the boundary
   
   x = x0 + a.* t;
   y = y0 + b.*t;
   z = z0 + c.*t;

   for ii=1:length(t)
        if z(ii) > top_cm
            z(ii) = top_cm;
        end
        if z(ii) < bot_cm
            z(ii) = bot_cm;
        end
   end

    if isempty(z) %i.e no solution so only passes through the top
        z(1) = bot_cm;
        z(2) = top_cm;
    end

    if length(z) == 1 %it hits just at a corner thus one solution, this is unlikely. 
        S = 0;
    else
        
        S = sqrt((x(1) - x(2)).^2 + (y(1) - y(2)).^2 + (z(1) - z(2)).^2);
    end

    if S==0
        S=NaN;
    end
    
    L = S;
    PE = [x(1) y(1) z(1)];
    PL = [x(2) y(2) z(2)];
end


function rq = fillS1RQs(i, data, rq)
  %  """ Get some S1 RQs from the data make make LIP-rqs for event i """
   % #get some S1 info
    pulsetypes = data.pulse_classification(: , i);
    pulse_area_phe = data.pulse_area_phe(: , i);

    pulse_filter =( pulsetypes == 1) == 1 ;
    n = sum(pulse_filter);
    rq.nS1(i) = n;

    if n > 0
        s1pas = pulse_area_phe(pulse_filter);
        rq.firstS1PA(i) = s1pas(1);
        rq.sumAllS1PA(i) = sum(s1pas);
        rq.largestS1PA(i) = max(s1pas);
        
        %#get the cumulative area before the first 10 S1s
        %#(maybe useful for cutting e-burps)
        pulse_indices_S1s = find(pulse_filter ==1);
        for ii = 1: length(pulse_indices_S1s)
            ind = pulse_indices_S1s(ii);
            if ii > 10
                break;
            end
            if ind == 1
                rq.PAbeforeS1(ii) = 0;
            else
                rq.PAbeforeS1(ii) = sum(pulse_area_phe(1:ind));
            end
        end
    end
end

function [max_z_score, min_z_score] = getMaxAndMinZScore(myarray)
    %""" return the max and min measure of how many standard deviations each element in input array is from mean """
    std_data = std(myarray);
    avg = mean(myarray);
    
    std_arr = zeros(length(myarray), 1);
    std_arr(:) = std_data;
    avg_arr = zeros(length(myarray),1);
    avg_arr(:) = avg;
    zscores = (myarray - avg_arr) ./ std_arr;

    max_z_score = max(zscores);
    min_z_score = min(zscores); %#note: this can be negative
    
    %return max_z_score, min_z_score
end
