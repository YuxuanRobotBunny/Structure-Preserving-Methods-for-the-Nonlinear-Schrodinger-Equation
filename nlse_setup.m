function P = nlse_setup(varargin)

    % defaults
    P = struct('lambda',1, 'L',20, 'N',256, 'T',2, 'tau',1e-3, ...
               'eta',1, 'xi',1, 'dealias',false, 'dealias_fraction',2/3);
    if mod(numel(varargin), 2) ~= 0
        error('nlse_setup:NameValuePairs', ...
              'Overrides must be supplied as name-value pairs.');
    end

    valid_names = fieldnames(P);

    % override
    for i = 1:2:numel(varargin)
        name = varargin{i};
        if ~ischar(name) && ~(isstring(name) && isscalar(name))
            error('nlse_setup:InvalidName', ...
                  'Override names must be character vectors or scalar strings.');
        end
        name = char(name);
        if ~ismember(name, valid_names)
            error('nlse_setup:UnknownParameter', ...
                  'Unknown parameter "%s". Valid names are: %s.', ...
                  name, strjoin(valid_names, ', '));
        end
        P.(name) = varargin{i+1};
    end

    if ~isscalar(P.N) || P.N ~= round(P.N) || P.N <= 0 || mod(P.N, 2) ~= 0
        error('nlse_setup:InvalidN', ...
              'P.N must be a positive even integer for the FFT wave-number layout.');
    end
    if ~isscalar(P.L) || P.L <= 0
        error('nlse_setup:InvalidL', 'P.L must be positive.');
    end
    if ~isscalar(P.T) || P.T <= 0
        error('nlse_setup:InvalidT', 'P.T must be positive.');
    end
    if ~isscalar(P.tau) || P.tau <= 0
        error('nlse_setup:InvalidTau', 'P.tau must be positive.');
    end
    if ~isscalar(P.lambda) || P.lambda <= 0
        error('nlse_setup:InvalidLambda', ...
              'P.lambda must be positive for the bright-soliton initial condition.');
    end
    if ~isscalar(P.eta) || P.eta <= 0
        error('nlse_setup:InvalidEta', 'P.eta must be positive.');
    end
    if ~isscalar(P.xi)
        error('nlse_setup:InvalidXi', 'P.xi must be scalar.');
    end
    if ~isscalar(P.dealias) || ...
       ~(islogical(P.dealias) || (isnumeric(P.dealias) && any(P.dealias == [0 1])))
        error('nlse_setup:InvalidDealias', ...
              'P.dealias must be a scalar logical or 0/1 numeric flag.');
    end
    P.dealias = logical(P.dealias);
    if ~isscalar(P.dealias_fraction) || P.dealias_fraction <= 0 || P.dealias_fraction > 1
        error('nlse_setup:InvalidDealiasFraction', ...
              'P.dealias_fraction must be in the interval (0, 1].');
    end

    % spatial grid: x_j = -L + j*h, j=0..N-1 (periodic, excludes +L)
    P.h = 2*P.L / P.N;
    P.x = (-P.L + (0:P.N-1)' * P.h);

    % wavenumbers (MATLAB FFT order: 0..N/2-1, -N/2..-1)
    k_pos = (0 : P.N/2-1)';
    k_neg = (-P.N/2 : -1)';
    P.k   = (pi/P.L) * [k_pos; k_neg];   % since 2*pi/(2L) = pi/L
    P.k2  = P.k.^2;
    mode_index = [k_pos; k_neg];
    cutoff = floor(P.dealias_fraction * (P.N/2));
    P.dealias_mask = abs(mode_index) <= cutoff;

    % bright soliton initial condition
    P.psi0 = sqrt(2*P.eta^2/P.lambda) ...
             .* sech(P.eta * P.x) .* exp(1i * P.xi * P.x);

    % adjust tau so Nsteps*tau == T exactly
    P.Nsteps = round(P.T / P.tau);
    if P.Nsteps < 1
        error('nlse_setup:InvalidNsteps', ...
              'P.tau is too large for P.T; expected round(P.T/P.tau) >= 1.');
    end
    P.tau    = P.T / P.Nsteps;
end
