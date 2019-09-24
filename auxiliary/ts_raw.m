function [Ms, Vs] = ts_raw(M, G)% Note, M here is actually number of players.
% Demo of the TrueSkill model
%PMTKauthor Carl Rasmussen and  Joaquin Quinonero-Candela,
%PMTKurl http://mlg.eng.cam.ac.uk/teaching/4f13/1112
%PMTKmodified Kevin Murphy
% Input:
% M = number of players
% G(i,1) = id of winner for game i
% G(i,2) = id of loser for hgame i
%
% Output:
% Ms(p) = mean of skill for player p
% Ps(p) = precision of skill for player p
N = size(G,1);            % number of games

psi = inline('normpdf(x)./normcdf(x)');
lambda = inline('(normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x)');

pv = 50;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and variances
Ms = nan(M,1); 
Ps = nan(M,1);

% initialize matrices of game to skill messages - means and precisions
Mgs = zeros(N,2); 
Pgs = zeros(N,2);

% allocate matrices of skill to game messages - means and precisions
Msg = nan(N,2); 
Psg = nan(N,2);

old_Ms = Ms;
old_Ps = Ps;
exit = 0;
iter = 0;
while (exit==0 && iter<2000)
  % (1) compute marginal skills 
  iter = iter+1;
  
  old_Ms = Ms;
  old_Ps = Ps;
  
  for p=1:M
    % compute this first because it is needed for the mean update
    Ps(p) = 1/pv + sum(Pgs(G==p)); 
    Ms(p) = sum(Pgs(G==p).*Mgs(G==p))./Ps(p);
  end
  
  if (any(any(isnan(Ms))) || any(any(isnan(Ps))) )
      Ms = old_Ms;
      Ps = old_Ps;
      break;
  end
  
  % (2) compute skill to game messages
  % compute this first because it is needed for the mean update
  Psg = Ps(G) - Pgs;
  %Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Ps(G); 
  Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Psg; % KPM 
     
  % (3) compute game to performance messages
  vgt = 1 + sum(1./Psg, 2);
  mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data
   
  % (4) approximate the marginal on performance differences
  Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
  Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );
    
  % (5) compute performance to game messages
  ptg = Pt - 1./vgt;
  ptg(ptg<0.001) = 0.001;
  
  mtg = (Mt.*Pt - mgt./vgt)./ptg;   
    
  % (6) compute game to skills messages
  Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
  Mgs = [mtg, -mtg] + Msg(:,[2 1]);

end

Vs = 1 ./ Ps;
end
