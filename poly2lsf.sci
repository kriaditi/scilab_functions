// ADITI KUMARI
// FUNCTION :POLY2LSF
// input argument: a- prediction polynomial
// output argument: lsf- line spectral frequencies
//Input:a=[1 .5 .2 .4 .3]
//output:x =[0.9273 1.2838 2.0944 2.6532]
function lsf = poly2lsf(a)
								      	    // if the coefficient of the highest order is not 1, we normalize.
 	if a(1)<>1 then
 		a=a/a(1)
 	end
	if max(abs(roots(a))) >= 1.0 then
 		error('The polynomial must have all the roots inside of the unit circle.')
 	end 	
								              //Creating the sum and the difference filters
 	p=length(a)-1  							//The first element is not used
 	x=[0]
 	a1=cat(2,a,x) 							//concatenate several arrays
 	a2=flipdim(a1,2,1) 						//flipping the components of the vector
 	P1 = a1 - a2   							//Difference filter
 	Q1 = a1 + a2   							//Sum filter
 	P2=P1(1) 
 	Q2=P2(1)
 	z=poly(0,"z")
 	for i=2:length(P1) 						//extracting length
 		P2=P2+P1(i)*z^(i)
 	end
 	for i=2:length(Q1)
 		Q2=Q2+Q1(i)*z^(i)
 	end
 	divisorOdd=-1*z^2+0*z^1+1
 	divisorEven1=-1*z^1+1
 	divisorEven2=1*z^1+1
 									//we have to remove the known roots at P1 and Q1, if order is even
  									// removing both the roots from P1, if odd
 	if modulo(p,2) then                                             //Odd 
 		[r,P]= pdiv(P2,divisorOdd)
 		Q=Q1
 	else                                                            //Even 
 		[r, P] = pdiv(P2, divisorEven1)
 		[r, Q] = pdiv(P2, divisorEven2)				//dividing
 	end
 	P=coeff(P)							//extracting coefficient
 	rP = roots(coeff(P))
 	rQ = roots(coeff(Q))						// to find roots


 	aP = atan(imag(rP), real(rP))         				//2-quadrant and 4-quadrant inverse tangent
 	aQ = atan(imag(rQ), real(rQ))
 	lsf = gsort(cat(1,-aP, -aQ),'r','i')				//sorting
 	lsf = lsf(find(lsf>0))
 	return lsf

 endfunction
