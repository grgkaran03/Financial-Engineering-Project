function [U] = FTCS(T, K, r, sig, delta, S, Tau, h, k, isTerminal, epsilon)
	fprintf('\nRunning FTCS\n');
	m = length(S);
	n = length(Tau);
	U = zeros(n, m);

	if isTerminal
		k = -k;
	end

	U(1:end, 1) = g1(r, S(1), Tau, K);
	U(1:end, end) = g2(r, S(end), Tau, K, epsilon);
	U(1, 1:end) = f(S, K);

	for i = 2:n
		for j = 2:m-1
			aa = fa(sig, S(j));
			bb = fb(r, delta, S(j));
			cc = fc(r);
			U(i, j) = (-aa*k/h^2 + 0.5*bb*k/h)*U(i-1,j-1) + (1 + 2*aa*k/h^2 - cc*k)*U(i-1,j) + (-aa*k/h^2 - 0.5*bb*k/h)*U(i-1,j+1);
		end
	end

	if isTerminal
		U = flipud(U);
	end
end

function [y] = f(S, K)
	temp1 = zeros(size(S));
	temp2 = S - K;
	y = max([temp1; temp2]);
end

function [y] = g1(r, s, Tau, K)
	y = 0;
end

function [y] = g2(r, s, Tau, K,epsilon)
	y = Smooth(s - K*exp(-r.*Tau),epsilon);
end

function [y] = fa(sig, s)
	y = sig.^2 .* s.^2 / 2;
end

function [y] = fb(r, delta, s)
	y = (r - delta) .* s;
end

function [y] = fc(r)
	y = -r;
end

function y = Smooth(x, epsilon)
    c0 = 35 / (256 * epsilon);
    c1 = 1 / 2;
    c2 = 35 / (64 * epsilon);
    c4 = -35 / (128 * epsilon^3);
    c6 = 7 / (64 * epsilon^5);
    c8 = -5 / (256 * epsilon^7);
    
   if (x >= epsilon)
       y = x;
   elseif (x <= -epslon)
       y = 0;
   else
       y = c0 + c1*y + c2*y^2 + c4*y^4 + c6*y^6 + c8*y^8;
   end
end


