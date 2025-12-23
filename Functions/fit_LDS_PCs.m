function [Max_eigen,idx_dir2,idx_pos2,r2_all,all_eigen,All_A]=fit_LDS_PCs(score,idx_dir,idx_pos,variance,threshold,binsl)
% FIT_LDS_PCS Fits linear dynamical systems (LDS) to neural trajectories in
% principal component space across task conditions.
%
%   [Max_eigen,idx_dir2,idx_pos2,r2_all,all_eigen,All_A] = ...
%       FIT_LDS_PCS(score,idx_dir,idx_pos,variance,threshold,binsl)
%
%   The function fits local linear dynamical systems to low-dimensional
%   neural trajectories obtained from principal component analysis (PCA).
%   LDS models are fitted in sliding temporal windows for each combination
%   of movement direction and duration. The dominant eigenvalue of the
%   fitted dynamics is used to characterise local system behaviour.Fits 
%   with poor reconstruction performance are discarded. 
%   The maximum eigenvalue is used as a summary measure of local 
%   dynamical behaviour.
%
%   INPUTS
%   ------
%
%   score : [T x D] matrix
%       Low-dimensional neural activity (e.g. PCA scores), where T is the
%       number of time bins and D is the total number of principal
%       components available.
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   variance : [D x 1] vector
%       Explained variance for each principal component.
%
%   threshold : scalar
%       Cumulative variance threshold used to determine the number of
%       principal components retained for LDS fitting. At least three
%       dimensions are always used.
%
%   binsl : scalar or empty
%       Half-length of the temporal window (in bins) used to fit each local
%       LDS. If empty, the window length is set to one quarter of the total
%       number of time points for the current condition.
%
%   OUTPUTS
%   -------
%
%   Max_eigen : [T x 1] vector
%       Maximum eigenvalue of the fitted dynamics matrix for each time bin.
%       Ill-fitted segments (r² < 0.3) or purely real eigenvalues are set to
%       NaN.
%
%   idx_dir2 : [T x 1] vector
%       Direction index corresponding to each entry of Max_eigen.
%
%   idx_pos2 : [T x 1] vector
%       Position index corresponding to each entry of Max_eigen.
%
%   r2_all : [T x 1] vector
%       Squared correlation between the
%       reconstructed and original trajectories for each LDS fit.
%
%   all_eigen : cell array
%       Cell array containing the full spectrum of eigenvalues for each
%       direction–duration condition.
%
%   All_A : [D x D x T] array
%       Estimated dynamics matrices A for each fitted LDS and time bin.
%
%
% Andrea Colins Rodriguez
% 05/05/2023
%
% 09/05/2023
% % Note: for matlab imag(nan)=0; We need to especify a
% % "complex nan" as complex(nan,nan) for specific cases. Do not take the
% % average of the entire complex number, take the separate parts.


ndim=max(find(cumsum(variance)>=threshold,1,'First'),3); % at least 3 dimensions, more if the subspace is high-dimensional

Ndir=numel(unique(idx_dir));
Nbins=numel(unique(idx_pos));

Max_eigen=nan(size(idx_dir,1),1);
idx_dir2=nan(size(idx_dir,1),1);
idx_pos2=nan(size(idx_dir,1),1);
r2_all=nan(size(idx_dir,1),1);
all_eigen=cell(Nbins,1);

counter_cond=1;

for i_dir=1:Ndir
      
    for i_dur=1:Nbins
        
        idx_this=idx_dir==i_dir & idx_pos==i_dur;
        This_score=score(idx_this,1:ndim);

         
        nTpoints=size(This_score,1);
        eigsA=nan(ndim,nTpoints);
        All_A=nan(ndim,ndim,nTpoints);
        r2=nan(nTpoints,1);

        % if the user does not choose a segment length, then  the length is
        % the total number of timebin/2

        if isempty(binsl)
        bins=round(nTpoints/4);
        else
            bins=binsl;
        end

        sLength=nan(nTpoints,1);
        sLength(bins+1:end-bins-1)=bins;       

        parfor i_time=1:nTpoints
            if ~isnan(sLength(i_time))

                [eigsA(:,i_time),r2(i_time),~,All_A(:,:,i_time)]=LDS(This_score(i_time-sLength(i_time):i_time+sLength(i_time),:));
                
           end
        end
        
        max_eig=max(eigsA,[],1);
        all_eigen{counter_cond}=eigsA;
         
         counter_cond=counter_cond+1;

         % ignore all ill fitted segments 
         % if squared correlation is below 0.3
        max_eig(r2<0.3)=complex(nan,nan);
        
        %% if the imaginary part is zero, then replace by nan
        
        % Note: for matlab imag(nan)=0;
        % that's silly! so we change it for nan+nan*i
         
        idx_imag = imag(max_eig)==0 & isnan(real(max_eig));
        max_eig(idx_imag)=complex(real(max_eig(idx_imag)),nan);

        
        % store outputs
        Max_eigen(idx_this)=max_eig;
        idx_dir2(idx_this)=idx_dir(idx_this);
        idx_pos2(idx_this)=idx_pos(idx_this);
        r2_all(idx_this)=r2;

    end
   
    
end
end