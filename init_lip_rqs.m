function  rq  = init_lip_rqs (nEvents, nPulses)

     %   #RQs about Pulse Area
    
    %   rq.sumAllPulsePA = zeros(1, nEvents);
    %     #   S2 Pulse Area
              
       
       rq.sumCorrS2PA = zeros(1, nEvents); 
       rq.firstCorrS2PA = zeros(1, nEvents);
       rq.lastCorrS2PA = zeros(1, nEvents);
       rq.maxCorrS2PA = zeros(1, nEvents);
       rq.minCorrS2PA = zeros(1, nEvents);
       rq.first2CorrS2PA = zeros(1, nEvents);
       rq.last2CorrS2PA = zeros(1, nEvents);
        
       rq.sumCorrS2PA_xyz_corr = zeros(1, nEvents); 
       rq.firstCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.lastCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.maxCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.minCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.first2CorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.last2CorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.moyalTest_xyzCorrPA = zeros(1, nEvents); 
       rq.sumCorrS2PA_xyz_corr_bot = zeros(1, nEvents); 
       rq.firstCorrS2PA_xyz_corr_bot = zeros(1, nEvents);
       rq.lastCorrS2PA_xyz_corr_bot = zeros(1, nEvents);
       rq.maxCorrS2PA_xyz_corr_bot = zeros(1, nEvents);
       rq.minCorrS2PA_xyz_corr_bot = zeros(1, nEvents);
       rq.first2CorrS2PA_xyz_corr_bot = zeros(1, nEvents);
       rq.last2CorrS2PA_xyz_corr_bot = zeros(1, nEvents);
        
       rq.cumSumCorrS2PA = zeros(nPulses, nEvents);
       rq.cumSumCorrS2PA_xyz_corr = zeros(nPulses, nEvents);
       rq.cumSumCorrS2PA_xyz_corr_bot = zeros(nPulses, nEvents);
 
       rq.stdCorrS2PA = zeros(1, nEvents);
       rq.stdCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.stdCorrS2PA_xyz_corr_bot = zeros(1, nEvents);

       rq.maxZScoreCorrS2PA = zeros(1, nEvents);
       rq.maxZScoreCorrS2PA_xyz_corr = zeros(1, nEvents);
       rq.maxZScoreCorrS2PA_xyz_corr_bot = zeros(1, nEvents);
   
     %   #   S1 Pulse Area
       rq.firstS1PA = zeros(1, nEvents);
       rq.sumAllS1PA = zeros(1, nEvents);
       rq.largestS1PA = zeros(1, nEvents);
        
  %      #amount of pulse area before each S1
 %       #only keep 10
       rq.PAbeforeS1 = zeros(10, nEvents);
        
        %#RQs about numbers of pulses
       rq.nS1 = zeros(1, nEvents);
       rq.nCorrS2FitPoints = zeros(1, nEvents); 
       rq.nCorrS2PulsesBeforeCut = zeros(1, nEvents);
  
   %     #Track RQs
 %       #  from fitting
 %       #regular "chi2" fit
       rq.CorrS2Chi2 = zeros(1, nEvents);
       rq.CorrS2Chi2_dof = zeros(1, nEvents);
        

     %   #angle*distance penalty "chi2" fit
   %    rq.CorrS2Chi2_anglePen = zeros(1, nEvents);
     %  rq.CorrS2Chi2_anglePen_dof = zeros(1, nEvents);
        

  %     rq.CorrS2TrackFit_anglePen_params = zeros(6, nEvents);

       rq.CorrS2TrackFit_params = zeros(6, nEvents);

       rq.CorrS2SVD_params = zeros(6, nEvents);

%        #theta and phi of the track
       rq.CorrS2TrackTheta = zeros(1, nEvents); 
       rq.CorrS2TrackPhi = zeros(1, nEvents); 
        

  %     rq.CorrS2TrackTheta_anglePen = zeros(1, nEvents); 
  %     rq.CorrS2TrackPhi_anglePen = zeros(1, nEvents); 
        

       % #entrance and exit points of the track
       rq.CorrS2TrackEnter = zeros(3, nEvents);
       rq.CorrS2TrackExit = zeros(3, nEvents);

 %      rq.CorrS2TrackEnter_anglePen = zeros(3, nEvents);
 %      rq.CorrS2TrackExit_anglePen = zeros(3, nEvents);

    %    #Track length of the fits
       rq.CorrS2Lfit = zeros(1, nEvents); 
   %    rq.CorrS2Lfit_anglePen = zeros(1, nEvents); 

        
       % #angle between "Data direction" and fit track direction
       rq.CorrS2DataToFitAng = zeros(1, nEvents); 

     %  rq.CorrS2DataToFitAng_anglePen = zeros(1, nEvents); 

        
      %  #no fits required
       rq.avgPointDistS2Corr = zeros(1, nEvents);
       rq.stdPointDistS2Corr = zeros(1, nEvents);
        
       rq.CorrS2Ldata = zeros(1, nEvents); 
        
       rq.stdCombThetaS2Corr = zeros(1, nEvents); 
       rq.stdCombPhiS2Corr = zeros(1, nEvents); 
        
       rq.stdPairThetaS2Corr = zeros(1, nEvents); 
       rq.stdPairPhiS2Corr = zeros(1, nEvents); 
        
       rq.avgCombThetaS2Corr = zeros(1, nEvents); 
       rq.avgCombPhiS2Corr = zeros(1, nEvents); 
        
       rq.avgPairThetaS2Corr = zeros(1, nEvents); 
       rq.avgPairPhiS2Corr = zeros(1, nEvents); 
       rq.luxstamp = zeros(1, nEvents); 
      % rq.eventWithinFile = zeros(1, nEvents);
       rq.eventWithinAcq = zeros(1, nEvents);
       rq.file_number = zeros(1, nEvents);
       %rq.luxstamp = uint64(rq.luxstamp);
       
       rq.raw_S2_total_keV = zeros(1, nEvents);
       rq.raw_total_keV = zeros(1, nEvents);
               
       rq.CorrS1_keV = zeros(1, nEvents);
       rq.CorrS2_keV = zeros(1, nEvents);
       rq.total_Corr_keV = zeros(1, nEvents);
       
end
