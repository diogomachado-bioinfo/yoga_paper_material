function smoothed_vector = fuzzy_moving_average(input_vector, vz, varargin)
% Improved fuzzy moving average using triangular function on both sides.
%
% Input:
%   - input_vector: Input vector.
%   - vz: Half of the window size.
%   - varargin: Optional arguments [hedge, ntimes].
%       - hedge: Hedge parameter for triangular function (default is 1).
%       - ntimes: Number of smoothing iterations (default is 1).
%
% Output:
%   - smoothed_vector: Smoothed vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    input_vector = vec_in_col(input_vector)';
    
    % Default values
    hedge = 1;
    ntimes = 1;

    % Check for optional arguments
    if ~isempty(varargin)
        vin = varargin{1};
        hedge = vin(1);
        ntimes = vin(2);
    end

    % Perform smoothing
    for itm = 1:ntimes
        n = length(input_vector(1, :));
        m = length(input_vector(:, 1));
        s = 2 * vz + 1;
        t = triang(s)'.^hedge;
        smoothed_vector = zeros(m, n);

        for i = 1:n
            d = (i - vz - 1);
            t0 = 1;
            if d <= 0
                k = vz + d;
                t0 = vz - k + 1;
            end

            d = (n - (i + vz));
            tE = s;
            if d <= 0
                tE = s + d;
            end

            i0 = max(i - vz, 1);
            iE = min(i + vz, n);

            a = input_vector(:, i0:iE);
            b = repmat(t(t0:tE), m, 1);
            c = sum(t(t0:tE));

            smoothed_vector(:, i) = sum(a .* b, 2) / c;
        end

        input_vector = smoothed_vector;
    end
end
