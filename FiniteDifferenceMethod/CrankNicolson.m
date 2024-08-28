function [U] = Crank(T, K, r, sig, delta, S, Tau, h, k, isTerminal, epsilon)
	fprintf('\nRunning Crank Nicolson\n');

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
		A = zeros(m, m);
		b = zeros(m, 1);

		aa = fa(sig, S);
		bb = fb(r, delta, S);
		cc = fc(r);

		A(1:m+1:end) = 2 - 2*aa*k/h^2 + cc*k;
		A(2:m+1:end) = aa(2:m)*k/h^2 - bb(2:m)*k/(2*h);
		A(m+1:m+1:end) = aa(1:m-1)*k/h^2 + bb(1:m-1)*k/(2*h);

		A(1,1) = 1;
		A(1,2) = 0;
		A(m,m) = 1;
		A(m,m-1) = 0;

		b(2:m-1) = (-aa(2:m-1)*k/h^2 + bb(2:m-1)*k/(2*h)) .* U(i-1,1:m-2) ...
					+ (2 + 2*aa(2:m-1)*k/h^2 - cc*k) .* U(i-1,2:m-1) ...
					+ (-aa(2:m-1)*k/h^2 - bb(2:m-1)*k/(2*h)) .* U(i-1,3:m);
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
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

function [y] = g2(r, s, Tau, K, epsilon)
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