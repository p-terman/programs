function output = Chi2Cal(sim_list)


for q =1:length(sim_list)
output(q).fitval=[];
    dp_name = strcat('C:\Users\Paul\Desktop\matfiles\', sim_list{q});
   % dp=load(dp_name);
dp = load_multi;
%starting with initalizing a bunch of vars where we can collect info along
%the way. Most are not output(q) at end of program. 
sys = 0.1; %1mm 
nEvents = size(dp.x_corrected,2);
%output(q).SumAllPulseArea = [];
%output(q).SumS2PulseArea = [];

output(q).Chi2_array = []; % #np.array([])
output(q).Chi2_dof_array =  []; % #np.array([])

output(q).SumAllPA = []; % #np.array([])
output(q).SumallS2PA = []; % #np.array([])
output(q).SumfitS2PA = [];% #np.array([])

output(q).firstS2PA = [];
output(q).lastS2PA = [];
output(q).x_int_array= [];
output(q).x_slope_array = [];
output(q).y_int_array= [];
output(q).y_slope_array = [];
%output(q).Lfit = [] ;% #given the fit, what is the projected path length through the detector
output(q).Ldata = [] ; %#what is the path length of the actual data

%NGoodPulses = [];
output(q).trackdot = [];

output(q).DeltaTheta = []; % #std of theta for each N set of vectors for N-S2 events
output(q).DeltaPhi = []; % #std of phi

output(q).DeltaE = []; % #std of S2s in fit
output(q).Delta_allS2 = []; % #std of all S2 in event
output(q).Max_Zscore_S2_Area = [];

output(q).firstS1PA = [];
output(q).SumAllS1PA = [];
output(q).largestS1PA = [];
output(q).nS1 = [];
output(q).oneS1mask = [];
output(q).fitmask = [];

output(q).SumAllPA = [];
output(q).SumfitS2PA = []; 

output(q).SumallS2PA = []; 

output(q).luxstamps = [];

output(q).DeltaTheta_3S2 = []; % #std of theta for each N set of vectors for N-S2 events
output(q).DeltaPhi_3S2 = []; % #std of phi

output(q).DeltaTheta = []; % #std of theta for each N set of vectors for N-S2 events
output(q).DeltaPhi = []; % #std of phi

output(q).DeltaE = [];
output(q).Delta_allS2 = [];
output(q).Max_Zscore_S2_Area = [];

nS2 = [];

sigma_x = [];
sigma_y = [];
sigma_z = [];
sigma_R = [];

output(q).sum_total_E_keV = [];
output(q).track_lengths = [];

numfit = 0;
count = 1;
figure
for i = 1 :nEvents
    %#get x,y,z
    
   dp=which_class5_type(dp);
   % dp = post_dpf_s1s2paring_corrections(dp);
    x_corr_cm = double(dp.x_corrected(:,i));
    y_corr_cm = double(dp.y_corrected(:,i));
    
    drift_to_cm = 1.51/10 ; % cm/us
    drift_us = 10 .*  0.001 .* double(dp.z_drift_samples(:,i));
    z_cm = 50 - drift_us .* drift_to_cm;  
    
    sd_radius_inf = dp.sd_radius_inf(:,i);
    sd_radius_sup = dp.sd_radius_sup(:,i);
    sd_phiXR = dp.sd_phiXR(:,i);
    z_err = ones(length(z_cm));
    
    sigma_r = [];
    sigma_r = sqrt( (0.5*(sd_radius_sup + sd_radius_inf)/sqrt(2.3)).^2 + sys.^2) ;
    
    
    phi_corr = atan2(y_corr_cm, x_corr_cm);
    r_corr_cm = sqrt(x_corr_cm.^2 + y_corr_cm.^2);
    
    pulse_area_phe = dp.pulse_area_phe(:,i);
    pulsetypes = dp.pulse_classification(:,i);
    
    
   % #get some S1 info
    pulse_filter_s1 = pulsetypes == 1 | dp.s1_like_class5(:,i) ==1;
    n = sum(pulse_filter_s1);
    output(q).nS1(i) = n;
    if n == 1
        output(q).oneS1mask(i)=true;
    else
        output(q).oneS1mask(i)= false;
    end

    if n > 0
        s1pas = pulse_area_phe(pulse_filter_s1);
        output(q).firstS1PA(i)=s1pas(1);
        output(q).SumAllS1PA(i) = sum(s1pas);
        output(q).largestS1PA(i) = max(s1pas);
    else
        output(q).firstS1PA(i) = 0;
        output(q).SumAllS1PA(i) = 0;
        output(q).largestS1PA(i) = 0;
    end
    
   % #get some S2 info
    pulse_filter_s2 = pulsetypes == 2  | pulsetypes == 4 | dp.s2_like_class5(:,i) == 1;
    n = sum(pulse_filter_s2);
   
    output(q).nS2(i) = n;
    if n > 0
        output(q).SumallS2PA(i) = sum(pulse_area_phe(pulse_filter_s2));
    else
        output(q).SumallS2PA(i)=(0);
    end
    
    
 %   #only plot S2s & Elses>10phd and SEs
    %pulse_filter = ((pulsetypes == 2) | (pulsetypes == 4) | ((pulsetypes == 5) & (pulse_area_phe > 10.))) &((x_corr_cm ~= -500.) & (y_corr_cm ~= -300.) )

    %only plot S2s 
    %pulse_filter = (pulsetypes == 2)  &((x_corr_cm != -500.) & (y_corr_cm != -300.) )
    pulse_filter = pulse_filter_s2  & (x_corr_cm > -99.);

    if sum(pulse_filter) > 2 %#can calculate DeltaTheta & DeltaPhi for n=3 events
        %#print 'ngood ', ngood[i]
        xarr = double(x_corr_cm(pulse_filter));
        yarr = double(y_corr_cm(pulse_filter));
        zarr = z_cm(pulse_filter);
        event = [xarr,yarr,zarr]'; % #row is (x,y,z)_i
        phis = [];
        thetas = [];
        
        %indices = comb_index(length(xarr),2);
        indices = combnk(1:length(xarr), 2);
        for k = 1:length(indices)
            ind = indices(k, :);
            %#print ind
            vec = event(:,ind(2))-event(:,ind(1));
            %#print 'vec = ', vec
            vec = vec ./norm(vec); %#unit vector in direction of propogation
            phi = atan2(vec(2),vec(1)) .* 180 ./ pi;
            theta = acos(vec(3)) .* 180 ./pi; %#arccos(z-component / length(=1))
            phis(k)= phi;
            thetas(k) = theta;
        end
 %       phis = np.array(phis);
  %      thetas = np.array(thetas);
       dth = std(thetas);
       dph = std(phis);
     %   #print phis, thetas
      %  #print np.var(phis), np.var(thetas)
        if sum(pulse_filter) < 5
            output(q).DeltaTheta_3S2(end+1)=dth;   
            output(q).DeltaPhi_3S2(end+1)= dph;
        end
        
        if sum(pulse_filter) > 4
            output(q).DeltaTheta_3S2(end+1) = dth;   
            output(q).DeltaPhi_3S2(end+1)=dph;  
            output(q).DeltaTheta(end+1)= dth;   
            output(q).DeltaPhi(end+1)=dph;   
        end
    end
      % # pa = np.array(pulse_area_phe[pulse_filter])
      %  #DeltaE.append(np.std(pa))
    if sum(pulse_filter) <=4
        output(q).fitmask(end+1)= false;
    end
    if sum(pulse_filter) > 4
        output(q).fitmask(i)=true;
        output(q).SumAllPA(i) = sum(pulse_area_phe);
        output(q).SumfitS2PA(i) = sum(pulse_area_phe(pulse_filter));
        pa = pulse_area_phe(pulse_filter);
        output(q).firstS2PA(i) = pa(1);
        output(q).lastS2PA(i) = pa(end);
        
        output(q).DeltaE(i) = std(pa) ;
        
         pt_s2 = pulse_filter_s2;
         pa_s2 = pulse_area_phe(pt_s2);
         output(q).Delta_allS2(i) = std(pa_s2);
        
         pa_s2_mean = mean(pa_s2);
         pa_s2_std = std(pa_s2);
   %     #zscore =(pa - pa_mean / pa_std)
         zscores_pa_s2 = (pa_s2 - (ones(size(pa_s2)).*pa_s2_mean)) ./ (ones(size(pa_s2)).*pa_s2_std);
         output(q).Max_Zscore_S2_Area(i) = max(zscores_pa_s2);

        
        data_vec = [xarr(end)-xarr(1), yarr(end)-yarr(1), zarr(end)-zarr(1)];
        output(q).Ldata(i) = norm( data_vec);

        numfit = numfit + 1;
     %   #only try to fit events with 3+ S2s
        fit_this_event = true;
        xarr = x_corr_cm(pulse_filter);
        yarr = y_corr_cm(pulse_filter);
        zarr = z_cm(pulse_filter);
        
      %  #c = lambda x: np.cos(x)
       % #s = lambda x: np.sin(x)
        cos_phi =  cos(phi_corr(pulse_filter)); 
        sin_phi =  sin(phi_corr(pulse_filter));
        d_phi = sd_phiXR(pulse_filter) ./ r_corr_cm(pulse_filter);

        xerr = sigma_r(pulse_filter) .* cos_phi - r_corr_cm(pulse_filter).*sin_phi .* d_phi;
        yerr = sigma_r(pulse_filter) .* sin_phi + r_corr_cm(pulse_filter).*cos_phi .* d_phi;
        xerr = abs(xerr);
        yerr = abs(yerr);
        
        zerr = z_err(pulse_filter); % these are all one since the z position is considered known
   
        sigma_x = xerr;
        sigma_y = yerr;
        sigma_z = zerr;
        sigma_R = sigma_r(pulse_filter);
    
    else
        continue
    end
 
   % #do fitting
    if fit_this_event
        %#fitline = SVD(xarr,yarr,zarr)
   %     #point1 = fitline[0]
    %    #point2 = fitline[-1]
        luxstamp = dp.luxstamp_samples(i);
    %   #maybe use SVD to get a guess for the fit?
    %   #params = np.array([1,1,1,0,1,-2]) #initial guess
        xx=xarr;
        yy=yarr;
        zz=zarr;
        params = SVD_params(xx, yy, zz);
        params=double(params); % the result of this line is going to be our 'inital guess'
                                %for the optimization of fminsearch.
      %  #print shape(params), shape(xarr), shape(yarr), shape(zarr), shape(xerr), shape(yerr), shape(zerr)
      %  #res = scipy.optimize.minimize(lf.chi2_componentwise, params, args=(xarr,yarr,zarr,xerr,yerr,zerr,false), method='nelder-mead', options={'xtol': 1E-8, 'disp':false})
        
       % res = scipy.optimize.minimize(lf.chi2_magwise, params, args=(xarr,yarr,zarr,xerr,yerr,zerr,false), method='nelder-mead', options={'xtol': 1E-8, 'disp':false})
        %options = optimset('param1',xarr,'param2',yarr,'param3',zarr,'param4', xerr, 'param5', yerr, 'param6', zerr );
       % init = [params];
      % fh=@(params, xarr,yarr,zarr,xerr,yerr,zerr)chi2_magwise
        params_opt=fminsearch(@chi2_magwise,params);%,xarr,yarr,zarr,xerr,yerr,zerr);%, options);
        %params_opt = res.x
       % fminunc(@chi2_magwise, params)
       % #params_opt = scipy.optimize.fmin(chi2_componentwise, params, xarr,yarr,zarr,xerr,yerr,zerr, disp=false)
        %output(q).params_array(end+1) = params_opt;
        [Chi2] = chi2_magwise(params_opt);%, xarr,yarr,zarr, xerr, yerr, zerr);
        output(q).Chi2_array(end+1) = Chi2;
        
        dof_3d_line = 4 ;%#2 for direction on the unitsphere, 3 for the point - 1 b/c you can translate the point
        chi2_dof = Chi2 / (length(xarr) - dof_3d_line); 
        output(q).Chi2_dof_array(end+1)= chi2_dof;            
        
          
%% Length and Energy calculations
        top_cm = 50; % cm 
        rmax =double(24.5); % cm this is the max radius of lux assuming circular shape
        drift_to_cm = 1.51/10 ; % cm/us
        
        x0=params(1);
        y0=params(2);
        z0=params(3);
        a=params(4);
        b=params(5);
        c=params(6);
        
        X0= x0 - a./c .* z0;
        A=a./c;
        Y0=y0 - b./c .* z0;
        B=b./c;
 
        syms z; % x y
        %eqns = [q == x0 + a .* z , y == y0 + m .* z , q.^2 + y.^2 == rmax.^2];
        % for some reason matlab acts weird when you have more eqns and a single
        % unknown. So I eliminated them 
        eqn =X0.^2 + 2.*A.*X0.*z + A.^2*z.^2 + Y0.^2 + 2.*B.*Y0.*z + B.^2.*z.^2 == rmax.^2;
        z = solve(eqn, z); % z is the positions where the particle hists the boundary
        z =double (z);

        for ii=1:length(z)
            if z(ii) > 32000 
                z(ii) = 32000;
            end
            if z(ii) < 0
                z(ii) = 0;
            end
        end

        if isempty(z)
            z(1) = 0;
            z(2) = 32000;
        end

        if length(z) == 1
            S = 0;
        else
            %z is down, so take negative of A and B. This seems to be the
            %only thing that makes sense when graphed. 
            A=-A;
            B=-B;
            x = X0 + A .* z;
            y = Y0 + B .* z;
            drift_us_z = 10 .*  0.001 .* z;
            z_cm = top_cm - drift_us_z .* drift_to_cm;
    
            S = sqrt((x(1) - x(2)).^2 + (y(1) - y(2)).^2 + (z_cm(1) - z_cm(2))^2);
            S=double(S);
        end

        if S==0
            S=NaN;
        end
        
        output(q).x_int_array(end+1) = [X0]; %this is x as a function of z
        output(q).x_slope_array(end+1) = [A];
        output(q).y_int_array(end+1) = [Y0];
        output(q).y_slope_array(end+1) = [B];
        output(q).track_lengths(end+1) = S;
        
        output(q).luxstamps(end+1) = luxstamp;
        
        s1_cut = (dp.pulse_classification(:,i) == 1 | dp.s1_like_class5(:,i)) & dp.s1s2_pairing(:,i) == 1;
        s2_cut = (dp.pulse_classification(:,i) == 2 | dp.pulse_classification(:,i) == 4 | dp.s2_like_class5(:,i) == 1) & dp.x_corrected(:,i) > -99;
        %these cuts are tighter than the previous cuts in this program.
        %Might want to change the earlier cuts to be like these?
        s1_E_keV = dp.xyz_corrected_pulse_area_all_phe(:,i) .* s1_cut * 13.7/1000 * 1/0.117;
        s2_E_keV = dp.xyz_corrected_pulse_area_all_phe(:,i) .* s2_cut * 13.7/1000 * 1/12.1;

        total_E_keV = s1_E_keV + s2_E_keV;

        output(q).sum_total_E_keV(end+1) = sum(total_E_keV);
        if count<11 && output(q).Chi2_dof_array(end)<10
            hold on 
            count = count + 1;
          
            RGB = [randi(256) randi(256) randi(256)]/256 ;
            scatter(xarr, -zarr, 15, RGB)
            p=[A, X0];
            fitval=polyval (p, -zarr);
            plot(fitval, -zarr, 'color' , RGB);
             %p = polyfit(-zarr, xarr,1)
              %fitval=polyval(p, -zarr);
             %  plot(fitval, -zarr)
            %count
       %     output(i).fitval(count)
        else 
            hold off;
        end
        


    end
end
end

    
%%

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
        
        xe = xerr(j);
        ye = yerr(j);
        ze = zerr(j);
        
    
        d = distance(params, x_val,y_val,z_val);

        err = [xe,ye,ze];
        
       % #temp = d/err
        chi2 = norm(d).^2 / norm(err).^2;
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

function  params = SVD_params(xx,yy,zz)
%    """ Get the Singular Value Decomposition of the Event. 
%        This returns a line that represents the average behavior of the event, 
%        like PCA """

    evt = [xx,yy,zz];

    datamean = mean(evt);
    x_init=datamean(1);
    y_init=datamean(2);
    z_init=datamean(3);

  %  # Do an SVD on the mean-centered data.
    v = pca(evt - datamean);
   
    x_slope=v(1,1);
    y_slope=v(2,1);
    z_slope=v(3,1);
    
    params = [x_init,y_init,z_init,x_slope,y_slope,z_slope];
end

end


