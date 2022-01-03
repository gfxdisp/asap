function [init_data] = ts_solve_fast(init_data)
% Note, M here is actually number of players.
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

% RKM 2021/12: Improved vectorization to make the computations with large
% number of conditions much faster

Ng = size(init_data.G,1);            % number of games
Nc = init_data.Nc;

pv = 0.5;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and variances
Ms = init_data.Ms;
Ps = 1./init_data.Vs;

% initialize matrices of game to skill messages - means and precisions
Mgs = [init_data.Mgs]; 
Pgs = [init_data.Pgs];

% allocate matrices of skill to game messages - means and precisions
Msg = nan(Ng,2); 
Psg = nan(Ng,2);
G = init_data.G;

for iter=1:init_data.n_iter

  % (2) compute skill to game messages
  % compute this first because it is needed for the mean update
  Psg = reshape(Ps(G),size(G)) - Pgs;
  Msg = (reshape(Ps(G).*Ms(G),size(G)) - Pgs.*Mgs)./Psg; 
     
  % (3) compute game to performance messages
  vgt = 1 + sum(1./Psg, 2);
  mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data
   
  % (4) approximate the marginal on performance differences
  Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
  Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );
    
  % (5) compute performance to game messages
  ptg = Pt - 1./vgt;
  mtg = (Mt.*Pt - mgt./vgt)./ptg;   
    
  % (6) compute game to skills messages
  Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
  Mgs = [mtg, -mtg] + Msg(:,[2 1]);
  
  % (1) compute marginal skills 
  PgsMgs = Pgs.*Mgs;
  Ps = ones(Nc,1)*1/pv;
  Ms = zeros(Nc,1);
  for kk=1:numel(Pgs)
      ind = G(kk);
      Ps(ind) = Ps(ind) + Pgs(kk);
      Ms(ind) = Ms(ind) + PgsMgs(kk);    
  end
  Ms = Ms./Ps;
  
end

Vs = 1 ./ Ps;

init_data.Ms=Ms;
init_data.Vs=Vs;
init_data.Mgs=Mgs;
init_data.Pgs=Pgs;


    function [y] = psi(x)
        y = normpdf(x)./normcdf(x);
    end

    function [y]= lambda(x)
        y = (normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x);
    end
end
